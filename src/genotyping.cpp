#include "genotyping.hpp"
#include "conversion.hpp"
#include "panmanUtils.hpp"
#include "panmap_utils.hpp"
#include "gap_map_utils.hpp"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <istream>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace genotyping;

void genotyping::stringSplit(const std::string& str, char delimiter, std::vector<std::string>& out) {
    out.clear();
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        if (!token.empty()) {
            out.push_back(token);
        }
    }
}

// Parse "size:prob" fields shared by the insertion (idx 4) and deletion (idx 5) matrix rows.
static void parseSizeProbFields(const std::vector<std::string>& fields,
                                std::unordered_map<int64_t, double>& out) {
    if (fields.empty()) {
        throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
    }
    for (const auto& f : fields) {
        std::vector<std::string> subFields;
        genotyping::stringSplit(f, ':', subFields);
        if (subFields.size() != 2) {
            throw std::invalid_argument("Invalid format in mutation matrix file");
        }
        int64_t size = std::stoll(subFields[0]);
        double prob = std::stod(subFields[1]);
        out[size] = prob;
    }
}

void genotyping::fillMutationMatricesFromFile(mutationMatrices& mutMat, std::ifstream& inf) {
    std::string line;
    int idx = 0;
    while (std::getline(inf, line)) {
        std::vector<double> probs;
        std::vector<std::string> fields;
        stringSplit(line, ' ', fields);

        if (fields.empty()) {
            break;
        }

        if (idx < 4) {
            if (fields.size() != 4) {
                throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
            }

            for (const auto& f : fields) {
                probs.push_back(std::stod(f));
            }
            mutMat.submat[idx] = std::move(probs);
        } else if (idx == 4) {
            parseSizeProbFields(fields, mutMat.insmat);
        } else if (idx == 5) {
            parseSizeProbFields(fields, mutMat.delmat);
        }
        idx++;
    }

    if (idx != 6) {
        throw std::invalid_argument("Received invalid mutation matrix (.mm) file");
    }

    mutMat.maxInsLogProb = 100.0;
    mutMat.maxDelLogProb = 100.0;
    if (!mutMat.insmat.empty()) {
        mutMat.maxInsLogProb =
            std::max_element(mutMat.insmat.begin(), mutMat.insmat.end(), [](const auto& a, const auto& b) {
                return a.second < b.second;
            })->second;
    }

    if (!mutMat.delmat.empty()) {
        mutMat.maxDelLogProb =
            std::max_element(mutMat.delmat.begin(), mutMat.delmat.end(), [](const auto& a, const auto& b) {
                return a.second < b.second;
            })->second;
    }

    mutMat.filled = true;
}


static constexpr int getIndexFromNucleotide(char nuc) {
    switch (nuc) {
        case 'A':
        case 'a': return 0;
        case 'C':
        case 'c': return 1;
        case 'G':
        case 'g': return 2;
        case 'T':
        case 't': return 3;
        case '*': return 4;
        default: return 5;
    }
}

std::vector<char> parse_alts(const std::string& alts_str) {
    std::vector<char> alts;
    std::istringstream alts_stream(alts_str);
    std::string alt;
    while (std::getline(alts_stream, alt, ',')) {
        if (alt.size() > 1) {
            throw std::runtime_error("Error: alt allel parsing error..");
        }
        alts.push_back(alt[0]);
    }
    return alts;
}

std::tuple<int, std::vector<double>, std::vector<int>, std::string, std::string>
parse_sample_formats(const std::string& sample_formats_str) {
    std::vector<std::string> sample_formats;
    std::istringstream sample_formats_stream(sample_formats_str);
    std::string sample_format;
    while (std::getline(sample_formats_stream, sample_format, ':')) {
        sample_formats.push_back(sample_format);
    }

    int gt = std::stoi(sample_formats[0]);
    std::vector<double> pls;
    std::istringstream pls_stream(sample_formats[1]);
    std::string pl;
    while (std::getline(pls_stream, pl, ',')) {
        pls.push_back(std::stod(pl));
    }

    std::vector<int> ads;
    std::istringstream ads_stream(sample_formats[2]);
    std::string ad;
    while (std::getline(ads_stream, ad, ',')) {
        ads.push_back(std::stoi(ad));
    }

    return std::make_tuple(gt, pls, ads, sample_formats[1], sample_formats[2]);
}

std::string genotyping::applyMutationSpectrum(const std::string& line,
                                              const std::vector<std::vector<double>>& scaled_submat) {
    std::vector<std::string> fields;
    std::istringstream line_stream(line);
    std::string field;
    while (std::getline(line_stream, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 10 || fields[0] == "#CHROM") {
        return line;
    }

    if (fields.size() != 10) {
        throw std::runtime_error("Couldn't parse VCF. Unrecognized number of fields.");
    }

    if (fields[4] == ".") {
        return "";
    }
    if (fields[7].substr(0, 2) != "DP" || getIndexFromNucleotide(fields[3][0]) > 3) {
        return fields[9][0] == '0' ? "" : line;
    }

    if (fields[3].size() > 1) {
        throw std::runtime_error("Error: reference allele parsing error.");
    }

    int ref_nuc_idx = getIndexFromNucleotide(fields[3][0]);
    std::vector<char> alts = parse_alts(fields[4]);

    auto [gt, pls, ads, pls_string, ads_string] = parse_sample_formats(fields[9]);

    std::vector<double> gls;
    if (alts.size() + 1 == pls.size()) {
        gls = pls;
    } else {
        for (int i = 0; i < pls.size(); i += 2) {
            gls.push_back(pls[i]);
            if (gls.size() == alts.size() + 1) {
                break;
            }
        }
    }

    gls[0] += scaled_submat[ref_nuc_idx][ref_nuc_idx];
    for (int i = 1; i < gls.size(); i++) {
        // scaled_submat is the 4x4 A/C/G/T substitution-spectrum prior, so only
        // apply it for base ALTs. A non-base ALT (e.g. '*' spanning deletion)
        // yields an index > 3 that would over-read the row; leave its likelihood
        // unmodified (matching how indels use the unmodified bcftools posterior).
        int alt_idx = getIndexFromNucleotide(alts[i - 1]);
        if (alt_idx <= 3) {
            gls[i] = gls[i] + scaled_submat[ref_nuc_idx][alt_idx];
        }
    }

    double min_gl = *std::min_element(gls.begin(), gls.end());
    int min_gl_index = 0;
    for (int i = 0; i < gls.size(); i++) {
        gls[i] -= min_gl;
        if (gls[i] == 0) min_gl_index = i;
    }

    if (min_gl_index == 0) return "";

    gt = min_gl_index;
    double qual = gls[0];

    std::stringstream ss;
    ss << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << fields[4] << "\t"
       << std::fixed << std::setprecision(4) << qual << "\t" << fields[6] << "\t" << fields[7] << "\t" << fields[8]
       << "\t" << gt << ":" << pls_string << ":" << ads_string;
    return ss.str();
}
