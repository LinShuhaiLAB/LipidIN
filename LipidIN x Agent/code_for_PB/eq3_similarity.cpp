#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;

struct LibraryPeak {
    double pmz;
    double mz;
    double intensity;
};

struct PreparedSpectrum {
    std::vector<double> mz;
    std::vector<double> intensity;
};

struct SimilarityResult {
    double jcr_similarity;
    double cosine_similarity;
    double intersection_over_union;
};

static std::string strip_dtl_prefix_cpp(const std::string& label) {
    const std::string key = "_DTL_";
    std::size_t pos = label.find(key);
    if (pos == std::string::npos) {
        return label;
    }
    return label.substr(pos + key.size());
}

static std::vector<double> parse_underscore_spectrum_cpp(const std::string& text) {
    std::vector<double> values;
    std::string item;
    std::stringstream ss(text);

    while (std::getline(ss, item, '_')) {
        if (item.empty()) {
            continue;
        }

        try {
            values.push_back(std::stod(item));
        } catch (...) {
            continue;
        }
    }

    return values;
}

static PreparedSpectrum prepare_spectrum_cpp(
    const std::vector<double>& mz_input,
    const std::vector<double>& int_input
) {
    PreparedSpectrum spectrum;

    std::size_t n = std::min(mz_input.size(), int_input.size());
    if (n == 0) {
        return spectrum;
    }

    spectrum.mz.resize(n);
    spectrum.intensity.resize(n);

    double max_intensity = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        spectrum.mz[i] = mz_input[i];
        spectrum.intensity[i] = int_input[i];
        if (int_input[i] > max_intensity) {
            max_intensity = int_input[i];
        }
    }

    if (max_intensity <= 0.0) {
        std::fill(spectrum.intensity.begin(), spectrum.intensity.end(), 0.0);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            spectrum.intensity[i] /= max_intensity;
        }
    }

    bool sorted = true;
    for (std::size_t i = 1; i < n; ++i) {
        if (spectrum.mz[i] < spectrum.mz[i - 1]) {
            sorted = false;
            break;
        }
    }

    if (!sorted) {
        std::vector<std::size_t> order(n);
        for (std::size_t i = 0; i < n; ++i) {
            order[i] = i;
        }

        std::sort(
            order.begin(),
            order.end(),
            [&](std::size_t a, std::size_t b) {
                return spectrum.mz[a] < spectrum.mz[b];
            }
        );

        std::vector<double> mz_sorted(n);
        std::vector<double> int_sorted(n);

        for (std::size_t i = 0; i < n; ++i) {
            mz_sorted[i] = spectrum.mz[order[i]];
            int_sorted[i] = spectrum.intensity[order[i]];
        }

        spectrum.mz.swap(mz_sorted);
        spectrum.intensity.swap(int_sorted);
    }

    return spectrum;
}

static double calculate_jcr_only_cpp(
    const PreparedSpectrum& sample,
    const PreparedSpectrum& library,
    double ppm_err2,
    double da_err2
) {
    if (sample.mz.empty() || library.mz.empty()) {
        return 0.0;
    }

    std::size_t i = 0;
    std::size_t j = 0;

    double min_sum = 0.0;
    double max_sum = 0.0;
    double ppm_scale = ppm_err2 / 1e6;

    while (i < sample.mz.size() && j < library.mz.size()) {
        double sample_mz_i = sample.mz[i];
        double lib_mz_j = library.mz[j];

        double tol = sample_mz_i * ppm_scale;
        if (tol < da_err2) {
            tol = da_err2;
        }

        double diff = sample_mz_i - lib_mz_j;

        if (diff >= -tol && diff <= tol) {
            double s = sample.intensity[i];
            double l = library.intensity[j];

            min_sum += std::min(s, l);
            max_sum += std::max(s, l);

            ++i;
            ++j;
        } else if (diff < -tol) {
            max_sum += sample.intensity[i];
            ++i;
        } else {
            max_sum += library.intensity[j];
            ++j;
        }
    }

    while (i < sample.mz.size()) {
        max_sum += sample.intensity[i];
        ++i;
    }

    while (j < library.mz.size()) {
        max_sum += library.intensity[j];
        ++j;
    }

    if (max_sum <= 0.0) {
        return 0.0;
    }

    return min_sum / max_sum;
}

static SimilarityResult calculate_full_similarity_cpp(
    const PreparedSpectrum& sample,
    const PreparedSpectrum& library,
    double ppm_err2,
    double da_err2,
    double jcr_similarity
) {
    SimilarityResult result;
    result.jcr_similarity = jcr_similarity;
    result.cosine_similarity = 0.0;
    result.intersection_over_union = 0.0;

    if (sample.mz.empty() || library.mz.empty()) {
        result.jcr_similarity = 0.0;
        return result;
    }

    double sample_sq_sum = 0.0;
    double lib_sq_sum = 0.0;
    double lib_peak_count = 0.0;

    for (double v : sample.intensity) {
        sample_sq_sum += v * v;
    }

    for (double v : library.intensity) {
        lib_sq_sum += v * v;
        if (v > 0.0) {
            lib_peak_count += 1.0;
        }
    }

    std::size_t i = 0;
    std::size_t j = 0;

    double dot_sum = 0.0;
    double matched_lib_peak_count = 0.0;
    double ppm_scale = ppm_err2 / 1e6;

    while (i < sample.mz.size() && j < library.mz.size()) {
        double sample_mz_i = sample.mz[i];
        double lib_mz_j = library.mz[j];

        double tol = sample_mz_i * ppm_scale;
        if (tol < da_err2) {
            tol = da_err2;
        }

        double diff = sample_mz_i - lib_mz_j;

        if (diff >= -tol && diff <= tol) {
            double s = sample.intensity[i];
            double l = library.intensity[j];

            if (s > 0.0 && l > 0.0) {
                dot_sum += s * l;
                matched_lib_peak_count += 1.0;
            }

            ++i;
            ++j;
        } else if (diff < -tol) {
            ++i;
        } else {
            ++j;
        }
    }

    if (sample_sq_sum > 0.0 && lib_sq_sum > 0.0) {
        result.cosine_similarity = dot_sum / std::sqrt(sample_sq_sum * lib_sq_sum);
    }

    // 库谱中有多少个峰能在样本中找到，按峰数计，不按强度加权
    if (lib_peak_count > 0.0) {
        result.intersection_over_union = matched_lib_peak_count / lib_peak_count;
    }

    return result;
}

class FastMS2Similarity {
public:
    void load_library(
        const std::vector<double>& db_pmz,
        const std::vector<double>& db_ms2_mz,
        const std::vector<double>& db_ms2_int,
        const std::vector<std::string>& db_labels
    ) {
        label_index.clear();

        std::size_t n = db_pmz.size();
        n = std::min(n, db_ms2_mz.size());
        n = std::min(n, db_ms2_int.size());
        n = std::min(n, db_labels.size());

        for (std::size_t i = 0; i < n; ++i) {
            std::string label = strip_dtl_prefix_cpp(db_labels[i]);

            LibraryPeak peak;
            peak.pmz = db_pmz[i];
            peak.mz = db_ms2_mz[i];
            peak.intensity = db_ms2_int[i];

            label_index[label].push_back(peak);
        }
    }

    std::vector<SimilarityResult> calculate_batch(
        const std::vector<std::string>& labels,
        const std::vector<double>& sample_pmz,
        const std::vector<double>& sample_pmz_da_error,
        const std::vector<std::string>& sample_ms2_mz_text,
        const std::vector<std::string>& sample_ms2_int_text,
        double ppm_err2,
        double da_err2,
        double jcr_cutoff
    ) const {
        std::size_t n = labels.size();
        n = std::min(n, sample_pmz.size());
        n = std::min(n, sample_pmz_da_error.size());
        n = std::min(n, sample_ms2_mz_text.size());
        n = std::min(n, sample_ms2_int_text.size());

        std::vector<SimilarityResult> results;
        results.reserve(n);

        for (std::size_t row_idx = 0; row_idx < n; ++row_idx) {
            SimilarityResult zero_result;
            zero_result.jcr_similarity = 0.0;
            zero_result.cosine_similarity = 0.0;
            zero_result.intersection_over_union = 0.0;

            auto it = label_index.find(labels[row_idx]);
            if (it == label_index.end()) {
                results.push_back(zero_result);
                continue;
            }

            std::vector<double> lib_mz;
            std::vector<double> lib_int;

            const auto& peaks = it->second;
            double target_error = sample_pmz_da_error[row_idx];
            double tolerance = 1e-10 + 1e-7 * std::fabs(target_error);

            std::set<std::pair<double, double>> seen_pairs;

            for (const auto& peak : peaks) {
                double pmz_error = std::fabs(peak.pmz - sample_pmz[row_idx]);
                if (std::fabs(pmz_error - target_error) <= tolerance) {
                    std::pair<double, double> key(peak.mz, peak.intensity);
                    if (seen_pairs.insert(key).second) {
                        lib_mz.push_back(peak.mz);
                        lib_int.push_back(peak.intensity);
                    }
                }
            }

            if (lib_mz.empty() || lib_int.empty()) {
                results.push_back(zero_result);
                continue;
            }

            std::vector<double> sample_mz_values = parse_underscore_spectrum_cpp(
                sample_ms2_mz_text[row_idx]
            );
            std::vector<double> sample_int_values = parse_underscore_spectrum_cpp(
                sample_ms2_int_text[row_idx]
            );

            PreparedSpectrum sample_spectrum = prepare_spectrum_cpp(
                sample_mz_values,
                sample_int_values
            );
            PreparedSpectrum library_spectrum = prepare_spectrum_cpp(
                lib_mz,
                lib_int
            );

            if (sample_spectrum.mz.empty() || library_spectrum.mz.empty()) {
                results.push_back(zero_result);
                continue;
            }

            double jcr_similarity = calculate_jcr_only_cpp(
                sample_spectrum,
                library_spectrum,
                ppm_err2,
                da_err2
            );

            if (jcr_similarity < jcr_cutoff) {
                zero_result.jcr_similarity = jcr_similarity;
                results.push_back(zero_result);
                continue;
            }

            SimilarityResult result = calculate_full_similarity_cpp(
                sample_spectrum,
                library_spectrum,
                ppm_err2,
                da_err2,
                jcr_similarity
            );

            results.push_back(result);
        }

        return results;
    }

private:
    std::unordered_map<std::string, std::vector<LibraryPeak>> label_index;
};

PYBIND11_MODULE(eq3_similarity, m) {
    py::class_<SimilarityResult>(m, "SimilarityResult")
        .def_readonly("jcr_similarity", &SimilarityResult::jcr_similarity)
        .def_readonly("cosine_similarity", &SimilarityResult::cosine_similarity)
        .def_readonly("intersection_over_union", &SimilarityResult::intersection_over_union);

    py::class_<FastMS2Similarity>(m, "FastMS2Similarity")
        .def(py::init<>())
        .def("load_library", &FastMS2Similarity::load_library)
        .def("calculate_batch", &FastMS2Similarity::calculate_batch);
}