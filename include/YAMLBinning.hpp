#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>

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
        std::cerr << "ERROR: cannot open YAML file " << path << std::endl;
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
