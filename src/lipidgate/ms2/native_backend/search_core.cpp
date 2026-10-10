// Numeric search only. LipidGate's evidence, scoring and ranking stay in Python.
// Build without fast-math/FMA: tolerance endpoints must match Python float64.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#ifdef _WIN32
#define LG_EXPORT extern "C" __declspec(dllexport)
#else
#define LG_EXPORT extern "C" __attribute__((visibility("default")))
#endif

LG_EXPORT std::uint32_t lg_api_version() noexcept { return 1; }

// Template flags: 1 = glyceride chain support; 2 = precursor ion.
// Record flags: 1 = positive glyceride prefilter; 2 = precursor-only FA scorer.
// Output indexes refer to occurrences, including repeated identical fragments.
LG_EXPORT int lg_batch_match(
    const double* peak_mz, const double* peak_intensity, std::uint32_t peak_count,
    const double* template_mz, const std::uint8_t* template_flags,
    std::uint32_t template_count, const std::uint32_t* occurrences,
    std::uint32_t occurrence_count, const std::uint32_t* offsets,
    const std::uint8_t* record_flags, std::uint32_t record_count,
    const std::uint32_t* requested, std::uint32_t requested_count,
    double tolerance, std::uint8_t use_ppm, std::uint8_t prefilter,
    std::uint8_t* selected, std::int32_t* matches) noexcept {
    try {
        if (!std::isfinite(tolerance) || tolerance < 0 ||
            peak_count > static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max()) ||
            !offsets || (peak_count && (!peak_mz || !peak_intensity)) ||
            (template_count && (!template_mz || !template_flags)) ||
            (occurrence_count && (!occurrences || !matches)) ||
            (record_count && !record_flags) ||
            (requested_count && (!requested || !selected))) return 1;
        if (offsets[0] != 0 || offsets[record_count] != occurrence_count) return 1;
        for (std::uint32_t r = 0; r < record_count; ++r)
            if (offsets[r] > offsets[r + 1]) return 1;
        for (std::uint32_t i = 0; i < occurrence_count; ++i)
            if (occurrences[i] >= template_count) return 1;
        for (std::uint32_t i = 0; i < requested_count; ++i)
            if (requested[i] >= record_count) return 1;
        if (occurrence_count) std::fill(matches, matches + occurrence_count, -1);
        if (requested_count) std::fill(selected, selected + requested_count, 0);

        std::vector<double> lower(template_count), upper(template_count);
        for (std::uint32_t i = 0; i < template_count; ++i) {
            const double window = use_ppm ? std::abs(template_mz[i]) * tolerance * 1e-6 : tolerance;
            lower[i] = template_mz[i] - window;
            upper[i] = template_mz[i] + window;
        }
        // Shared template ranges are computed lazily once per query. Repeated
        // occurrences still match independently and consume separate peaks.
        std::vector<std::int32_t> left_indexes(template_count, -1), right_indexes(template_count, -1);
        const auto left_index = [&](std::uint32_t t) {
            auto& index = left_indexes[t];
            if (index < 0) index = static_cast<std::int32_t>(
                std::lower_bound(peak_mz, peak_mz + peak_count, lower[t]) - peak_mz);
            return index;
        };
        const auto right_index = [&](std::uint32_t t) {
            auto& index = right_indexes[t];
            if (index < 0) index = static_cast<std::int32_t>(
                std::upper_bound(peak_mz, peak_mz + peak_count, upper[t]) - peak_mz);
            return index;
        };
        std::vector<std::uint64_t> used(peak_count, 0);
        for (std::uint32_t q = 0; q < requested_count; ++q) {
            const auto r = requested[q];
            const auto begin = offsets[r], end = offsets[r + 1];
            bool overlap = !prefilter;
            bool precursor_only = false;
            for (auto j = begin; j < end; ++j) {
                const auto t = occurrences[j];
                if ((record_flags[r] & 2) && (template_flags[t] & 2)) precursor_only = true;
                if (overlap || !peak_count || ((record_flags[r] & 1) && !(template_flags[t] & 1))) continue;
                const auto p = left_index(t);
                if (p < static_cast<std::int32_t>(peak_count) && peak_mz[p] <= upper[t]) overlap = true;
            }
            if (!overlap) continue;
            selected[q] = 1;
            const auto generation = static_cast<std::uint64_t>(q) + 1;
            if (!peak_count) continue;
            for (auto j = begin; j < end; ++j) {
                const auto t = occurrences[j];
                if (precursor_only && !(template_flags[t] & 2)) continue;
                const auto left = left_index(t), right = right_index(t);
                std::int32_t best = -1;
                for (auto p = left; p < right; ++p) {
                    // Strict > keeps the first index on ties, exactly as Python.
                    if (used[p] != generation && (best < 0 || peak_intensity[p] > peak_intensity[best]))
                        best = static_cast<std::int32_t>(p);
                }
                if (best >= 0) { used[best] = generation; matches[j] = best; }
            }
        }
        return 0;
    } catch (...) { return 2; }
}
