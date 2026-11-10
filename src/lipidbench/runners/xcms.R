suppressWarnings(suppressMessages(library(MSnbase)))
suppressWarnings(suppressMessages(library(xcms)))
suppressMessages(library(optparse))
suppressMessages(library(BiocParallel))

option_list <- list(
  make_option(c("--stdin_config"), action = "store_true", default = FALSE),
  make_option(c("-d", "--dir"), type = "character", help = "mzML 文件夹路径"),
  make_option(c("-o", "--output"), type = "character", help = "输出文件路径"),
  make_option(c("--polarity"), type = "character", default = "positive"),
  make_option(c("--mz_tol"), type = "numeric", default = 10),
  make_option(c("--minwidth"), type = "numeric", default = 5),
  make_option(c("--maxwidth"), type = "numeric", default = 60),
  make_option(c("--noise"), type = "numeric", default = 1000),
  make_option(c("--sn"), type = "numeric", default = 3),
  make_option(c("--prefilter"), type = "numeric", default = 3),
  make_option(c("--frac"), type = "numeric", default = 0.2),
  make_option(c("--min_samples"), type = "integer", default = 1),
    make_option(c("--mzdiff"), type = "numeric", default = 0.001),
    make_option(c("--min_maxo"), type = "numeric", default = NA)
)
register(SerialParam(progressbar = TRUE))
opt <- parse_args(OptionParser(option_list = option_list))
if (isTRUE(opt$stdin_config)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install R package jsonlite for sample metadata")
    connection <- file("stdin", encoding = "UTF-8")
    request <- jsonlite::fromJSON(paste(readLines(connection, warn = FALSE), collapse = "\n"), simplifyVector = FALSE)
    close(connection)
    for (name in names(request)) opt[[name]] <- request[[name]]
}

files <- if (!is.null(opt$files)) unlist(opt$files, use.names=FALSE) else list.files(opt$dir, pattern = "\\.mzML$", full.names = TRUE, ignore.case = TRUE)
raw_data = readMSData(files,mode="onDisk",msLevel=1)
raw_data <- filterEmptySpectra(raw_data)
params <- CentWaveParam(
    ppm = opt$mz_tol,
    peakwidth = c(opt$minwidth, opt$maxwidth),
    noise = opt$noise,
    snthresh = opt$sn,
    mzdiff = opt$mzdiff,
    prefilter = c(opt$prefilter, 3000),
    mzCenterFun = "wMean",
    integrate = 1,
    fitgauss = FALSE
)
xdata <- findChromPeaks(raw_data, params, return.type="XCMSnExp")

# Retain actual local IDs, raw RTs, boundaries and areas before RT correction.
local <- chromPeaks(xdata)
native <- data.frame(Feature_ID = if (nrow(local)) rownames(local) else character(), source_file = basename(files[as.integer(local[, "sample"])]),
    mz = local[, "mz"], RT = local[, "rt"] / 60, RTmin = local[, "rtmin"] / 60,
    RTmax = local[, "rtmax"] / 60, intensity = local[, "into"], maxo = local[, "maxo"],
    RT_unit = rep("minutes", nrow(local)), check.names = FALSE)
native_path <- sub("\\.csv$", "_native.csv", opt$output, ignore.case = TRUE)
write.csv(native, native_path, row.names = FALSE, fileEncoding = "UTF-8")
if (!nrow(local)) {
    empty <- data.frame(Row.names=character(), mzmed=numeric(), rtmin=numeric(), rtmed=numeric(), rtmax=numeric(), npeaks=integer())
    for (name in basename(files)) empty[[name]] <- numeric()
    write.csv(empty, opt$output, row.names=FALSE, fileEncoding="UTF-8")
    quit(save="no", status=0)
}

groups <- rep("1", length(files))
if (!is.null(opt$samples)) {
    records <- opt$samples
    has_samples <- any(vapply(records, function(row) isTRUE(row$included) && row$sample_type == "Sample", logical(1)))
    has_qc <- any(vapply(records, function(row) isTRUE(row$included) && row$sample_type == "QC", logical(1)))
    groups <- vapply(basename(files), function(name) {
        matching <- Filter(function(row) row$source_file == name && isTRUE(row$included), records)
        if (length(matching) != 1L) stop("Missing or ambiguous sample metadata: ", name)
        row <- matching[[1]]
        if ((has_samples && row$sample_type != "Sample") || (!has_samples && has_qc && row$sample_type != "QC")) return(NA_character_)
        as.character(row$group)
    }, character(1))
}

pdp <- PeakDensityParam(sampleGroups = groups, minFraction = opt$frac,
                        minSamples = opt$min_samples, bw = 10, binSize = 0.005)
xdata <- groupChromPeaks(xdata, param = pdp)

if (length(files) > 1) {
    xdata <- adjustRtime(xdata, ObiwarpParam(binSize = 1))
    xdata <- groupChromPeaks(xdata, param = pdp)
    xdata <- fillChromPeaks(xdata)
}

peakInfo <- featureDefinitions(xdata)
pks <- chromPeaks(xdata)

peak_maxo <- sapply(peakInfo$peakidx, function(idx) {
    if (length(idx) == 0) return(NA_real_)
    vals <- pks[idx, "maxo"]
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(NA_real_)
    max(vals)
})

peak_into <- sapply(peakInfo$peakidx, function(idx) {
    if (length(idx) == 0) return(NA_real_)
    vals <- pks[idx, "into"]
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(NA_real_)
    max(vals)
})

peakInfo$maxo <- peak_maxo
peakInfo$into <- peak_into

if (length(files) == 1) {
    peak_indices <- sapply(peakInfo$peakidx, function(x) x[1])
    peakInfo$mzmin <- pks[peak_indices, "mzmin"]
    peakInfo$mzmax <- pks[peak_indices, "mzmax"]
    peakInfo$rtmin <- pks[peak_indices, "rtmin"]
    peakInfo$rtmax <- pks[peak_indices, "rtmax"]
}

intensities <- featureValues(xdata, value = "into")
snrs <- featureValues(xdata, value = "sn")
peaktable <- merge(peakInfo, intensities, by = 'row.names' , all = TRUE)

is_single_sample <- is.null(dim(snrs)) || (length(dim(snrs)) == 2 && ncol(snrs) == 1)
if (is_single_sample) {
    snrs_df <- as.data.frame(snrs)
    colnames(snrs_df) <- "sn"
    peaktable <- merge(peaktable, snrs_df, by.x = "Row.names", by.y = "row.names", all.x = TRUE)
}

drop_cols <- c("peakidx")
peaktable <- peaktable[, !(names(peaktable) %in% drop_cols)]

if (!is.na(opt$min_maxo) && "maxo" %in% names(peaktable)) {
    peaktable <- peaktable[!is.na(peaktable$maxo) & peaktable$maxo >= opt$min_maxo, , drop = FALSE]
}

write.table(peaktable, file = opt$output, sep=",", row.names= FALSE, col.names=TRUE, quote=TRUE, fileEncoding="UTF-8")
