#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "Utility.hpp"

struct BinDef {
    int bin_id = -1;
    double Q2_min = -1.0;
    double Q2_max = -1.0;
    double beta_min = -1.0;
    double beta_max = -1.0;
    double xpom_min = -1.0;
    double xpom_max = -1.0;
};

inline bool StartsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

inline std::string TrimWS(const std::string& s) {
    const char* ws = " \t\r\n";
    const size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

inline std::vector<double> ParseInlineList(const std::string& line) {
    std::vector<double> values;
    const auto lbr = line.find('[');
    const auto rbr = line.find(']');
    if (lbr == std::string::npos || rbr == std::string::npos || rbr <= lbr) return values;
    std::string body = line.substr(lbr + 1, rbr - lbr - 1);
    for (char& c : body) {
        if (c == ',') c = ' ';
    }
    std::istringstream iss(body);
    double v = 0.0;
    while (iss >> v) {
        values.push_back(v);
    }
    return values;
}

inline std::vector<double> ReadInlineListFromYAML(const std::string& path, const std::string& key) {
    std::ifstream in(path);
    if (!in.is_open()) return {};
    std::string line;
    const std::string prefix = key + ":";
    while (std::getline(in, line)) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = TrimWS(line);
        if (line.empty()) continue;
        if (StartsWith(line, prefix)) {
            return ParseInlineList(line);
        }
    }
    return {};
}

inline std::vector<BinDef> ReadBinsFromYAML(const std::string& path) {
    std::vector<BinDef> bins;
    std::ifstream in(path);
    if (!in.is_open()) {
        Logger::error("cannot open YAML file " + path);
        return bins;
    }
    bool in_bins = false;
    BinDef cur;
    bool have_cur = false;
    std::string line;
    while (std::getline(in, line)) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = TrimWS(line);
        if (line.empty()) continue;
        if (StartsWith(line, "3D_bins:")) {
            in_bins = true;
            continue;
        }
        if (!in_bins) continue;
        if (StartsWith(line, "t_bins:")) continue;
        if (StartsWith(line, "W2_bins:")) continue;
        if (StartsWith(line, "-")) {
            if (have_cur) bins.push_back(cur);
            cur = BinDef();
            have_cur = true;
            const size_t pos = line.find("bin_id");
            if (pos != std::string::npos) {
                const size_t colon = line.find(':', pos);
                if (colon != std::string::npos) {
                    cur.bin_id = std::stoi(TrimWS(line.substr(colon + 1)));
                }
            }
            continue;
        }
        const size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        const std::string key = TrimWS(line.substr(0, colon));
        const std::string val_str = TrimWS(line.substr(colon + 1));
        if (val_str.empty()) continue;
        const double val = std::stod(val_str);
        if (key == "bin_id") cur.bin_id = static_cast<int>(val);
        else if (key == "Q2_min") cur.Q2_min = val;
        else if (key == "Q2_max") cur.Q2_max = val;
        else if (key == "beta_min") cur.beta_min = val;
        else if (key == "beta_max") cur.beta_max = val;
        else if (key == "x_min" || key == "x_pom_min" || key == "xpom_min") cur.xpom_min = val;
        else if (key == "x_max" || key == "x_pom_max" || key == "xpom_max") cur.xpom_max = val;
    }
    if (have_cur) bins.push_back(cur);
    return bins;
}

inline std::vector<double> UniqueSorted(std::vector<double> vals) {
    std::sort(vals.begin(), vals.end());
    const double eps = 1e-12;
    vals.erase(std::unique(vals.begin(), vals.end(), [eps](double a, double b) {
        return std::fabs(a - b) < eps;
    }), vals.end());
    return vals;
}

inline void CollectEdges(const std::vector<BinDef>& bins,
                         std::vector<double>& q2_edges,
                         std::vector<double>& beta_edges,
                         std::vector<double>& xpom_edges) {
    for (const auto& b : bins) {
        if (b.Q2_min > 0) q2_edges.push_back(b.Q2_min);
        if (b.Q2_max > 0) q2_edges.push_back(b.Q2_max);
        if (b.beta_min >= 0) beta_edges.push_back(b.beta_min);
        if (b.beta_max >= 0) beta_edges.push_back(b.beta_max);
        if (b.xpom_min > 0) xpom_edges.push_back(b.xpom_min);
        if (b.xpom_max > 0) xpom_edges.push_back(b.xpom_max);
    }
    q2_edges = UniqueSorted(q2_edges);
    beta_edges = UniqueSorted(beta_edges);
    xpom_edges = UniqueSorted(xpom_edges);
}

//==============================================================================
// Global-bin helpers and on-disk writers (shared by skim + plotter).
// Previously duplicated in DDIS_Skim_Final.cpp and DDIS_Plot_Final.cpp.
//==============================================================================

inline bool EdgesMatchWithinTolerance(double a, double b) {
    const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= 1e-12 * scale;
}

inline int FindLowerEdgeIndex(const std::vector<double>& edges, double value) {
    if (edges.size() < 2) return -1;
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        if (EdgesMatchWithinTolerance(edges[i], value)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// Convention matches the plotter: k = iBeta + iXpom*nBeta + iQ2*(nBeta*nXpom).
inline int GetGlobalBinFromBinDef(const BinDef& b,
                                  const std::vector<double>& q2_edges,
                                  const std::vector<double>& xpom_edges,
                                  const std::vector<double>& beta_edges) {
    const int iQ2 = FindLowerEdgeIndex(q2_edges, b.Q2_min);
    const int iXpom = FindLowerEdgeIndex(xpom_edges, b.xpom_min);
    const int iBeta = FindLowerEdgeIndex(beta_edges, b.beta_min);
    if (iQ2 < 0 || iXpom < 0 || iBeta < 0) return -1;

    if (iQ2 + 1 >= static_cast<int>(q2_edges.size()) ||
        iXpom + 1 >= static_cast<int>(xpom_edges.size()) ||
        iBeta + 1 >= static_cast<int>(beta_edges.size())) {
        return -1;
    }

    if (!EdgesMatchWithinTolerance(q2_edges[iQ2 + 1], b.Q2_max) ||
        !EdgesMatchWithinTolerance(xpom_edges[iXpom + 1], b.xpom_max) ||
        !EdgesMatchWithinTolerance(beta_edges[iBeta + 1], b.beta_max)) {
        return -1;
    }

    const int nBeta = static_cast<int>(beta_edges.size()) - 1;
    const int nXpom = static_cast<int>(xpom_edges.size()) - 1;
    return iBeta + iXpom * nBeta + iQ2 * (nBeta * nXpom);
}

inline bool WriteBinsTSV(const std::string& path, const std::vector<BinDef>& bins) {
    std::ofstream out(path);
    if (!out.is_open()) {
        Logger::error("cannot open " + path + " for writing bins TSV");
        return false;
    }
    out << "# bin_id\tQ2_min\tQ2_max\tbeta_min\tbeta_max\tx_pom_min\tx_pom_max\n";
    out << std::setprecision(8);
    for (const auto& b : bins) {
        out << b.bin_id << "\t"
            << b.Q2_min << "\t" << b.Q2_max << "\t"
            << b.beta_min << "\t" << b.beta_max << "\t"
            << b.xpom_min << "\t" << b.xpom_max << "\n";
    }
    return true;
}

inline bool WriteBinOccupancyYAML(const std::string& path,
                                  const std::vector<BinDef>& bins,
                                  const std::vector<int>& occupancy,
                                  const std::vector<double>& q2_edges,
                                  const std::vector<double>& xpom_edges,
                                  const std::vector<double>& beta_edges) {
    std::ofstream out(path);
    if (!out.is_open()) {
        Logger::error("cannot open " + path + " for writing bin occupancy YAML");
        return false;
    }

    long long totalAssigned = 0;
    for (const int c : occupancy) totalAssigned += c;

    out << "metadata:\n";
    out << "  description: \"Per-bin occupancy from truth global bin assignment (Q2_truth, x_pom_truth, beta_truth).\"\n";
    out << "  n_bins: " << bins.size() << "\n";
    out << "  total_assigned_events: " << totalAssigned << "\n";
    out << "  k_index_convention: \"one_based\"\n";
    out << "bins:\n";
    out << std::setprecision(8);

    for (const auto& b : bins) {
        const int kZeroBased = GetGlobalBinFromBinDef(b, q2_edges, xpom_edges, beta_edges);
        const int occ = (kZeroBased >= 0 && kZeroBased < static_cast<int>(occupancy.size()))
                            ? occupancy[kZeroBased]
                            : 0;
        out << "  - bin_id: " << b.bin_id << "\n";
        out << "    k_index: " << (kZeroBased >= 0 ? kZeroBased + 1 : -1) << "\n";
        out << "    occupancy: " << occ << "\n";
        out << "    Q2_min: " << b.Q2_min << "\n";
        out << "    Q2_max: " << b.Q2_max << "\n";
        out << "    x_pom_min: " << b.xpom_min << "\n";
        out << "    x_pom_max: " << b.xpom_max << "\n";
        out << "    beta_min: " << b.beta_min << "\n";
        out << "    beta_max: " << b.beta_max << "\n";
    }
    return true;
}
