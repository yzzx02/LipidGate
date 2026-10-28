// Numeric search only. LipidGate's evidence, scoring and ranking stay in Python.
// Build without fast-math/FMA: tolerance endpoints must match Python float64.
#include <algorithm>
#include <array>
#include <cstring>
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

LG_EXPORT std::uint32_t lg_plan_api_version() noexcept { return 1; }

// Immutable compact evidence columns. A program is an OR of alternatives,
// each alternative an AND of occurrence-mask population lower bounds.
// Passing a program merely admits the unchanged full Python policy scorer.
LG_EXPORT int lg_search_plan(
    const double* peak_mz, const double* peak_intensity, std::uint32_t peak_count,
    const double* masses, const std::uint8_t* fragment_flags,
    const std::uint32_t* offsets, const std::uint8_t* record_flags,
    const std::uint32_t* programs, std::uint32_t record_count,
    const std::uint32_t* plan_offsets, std::uint32_t plan_count,
    const std::uint32_t* arm_offsets, const std::uint64_t* masks, const std::uint8_t* minimums,
    std::uint32_t left, std::uint32_t right, double tolerance, std::uint8_t use_ppm,
    std::uint8_t prefilter, std::uint8_t* selected, std::int32_t* matches) noexcept {
    try {
        if (!std::isfinite(tolerance) || tolerance < 0 || left > right || right > record_count ||
            peak_count > static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max()) ||
            !offsets || !plan_offsets || !arm_offsets ||
            (record_count && (!programs || !record_flags)) ||
            (peak_count && (!peak_mz || !peak_intensity)) ||
            (right > left && !selected)) return 1;
        const auto base = offsets[left], end = offsets[right];
        if (base > end || (end > base && (!masses || !fragment_flags || !matches))) return 1;
        if (end > base) std::fill(matches, matches + end - base, -1);
        if (right > left) std::fill(selected, selected + right - left, 0);
        if (!peak_count) return 0;
        std::vector<std::uint64_t> used(peak_count, 0);
        struct Window {
            std::uint64_t key{};
            std::uint32_t left{}, right{};
            bool valid{};
        };
        // Headgroup and chain masses recur across nearby precursor candidates.
        // Cache only this invocation's exact search bounds, never intensities,
        // matches, gate decisions or results. Collisions simply replace a slot.
        std::array<Window,128> windows{};
        const bool reuse_windows=right-left>4;
        for (auto r = left; r < right; ++r) {
            const auto begin = offsets[r], stop = offsets[r + 1];
            if (begin > stop || stop > end || programs[r] >= plan_count) return 1;
            bool overlap = !prefilter;
            bool precursor_only = false;
            if (record_flags[r] & 2)
                for (auto j = begin; j < stop; ++j)
                    if (fragment_flags[j] & 2) { precursor_only = true; break; }
            std::uint64_t matched_mask = 0;
            const auto generation = static_cast<std::uint64_t>(r) + 1;
            for (auto j = begin; j < stop; ++j) {
                const double window = use_ppm ? std::abs(masses[j]) * tolerance * 1e-6 : tolerance;
                std::uint64_t bits;
                std::memcpy(&bits,masses+j,sizeof(bits));
                const auto hash=(bits^(bits>>33))*std::uint64_t(0xff51afd7ed558ccd);
                auto& cached=windows[(hash^(hash>>33)) & 127];
                const double *lo,*hi;
                if (reuse_windows && cached.valid && cached.key==bits) {
                    lo=peak_mz+cached.left;hi=peak_mz+cached.right;
                } else {
                    lo=std::lower_bound(peak_mz, peak_mz + peak_count, masses[j] - window);
                    hi=std::upper_bound(lo, peak_mz + peak_count, masses[j] + window);
                    if (reuse_windows) {
                        cached.key=bits;cached.left=static_cast<std::uint32_t>(lo-peak_mz);
                        cached.right=static_cast<std::uint32_t>(hi-peak_mz);cached.valid=true;
                    }
                }
                if (lo != hi && (!(record_flags[r] & 1) || (fragment_flags[j] & 1))) overlap = true;
                if (precursor_only && !(fragment_flags[j] & 2)) continue;
                std::int32_t best = -1;
                for (auto p = lo; p != hi; ++p) {
                    const auto index = static_cast<std::int32_t>(p - peak_mz);
                    if (used[index] != generation && (best < 0 || peak_intensity[index] > peak_intensity[best])) best = index;
                }
                if (best >= 0) {
                    used[best] = generation;
                    matches[j - base] = best;
                    if (j - begin < 64) matched_mask |= std::uint64_t(1) << (j - begin);
                }
            }
            if (!overlap) continue;
            if (record_flags[r] & 4) { selected[r-left] = 2; continue; }
            const auto program = programs[r];
            bool possible = false;
            for (auto a = plan_offsets[program]; a < plan_offsets[program + 1]; ++a) {
                bool passed = true;
                for (auto c = arm_offsets[a]; c < arm_offsets[a + 1]; ++c) {
                    if (__builtin_popcountll(matched_mask & masks[c]) < minimums[c]) { passed = false; break; }
                }
                if (passed) { possible = true; break; }
            }
            if (possible) selected[r-left] = 1;
        }
        return 0;
    } catch (...) { return 2; }
}

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
