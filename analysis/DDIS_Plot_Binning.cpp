// DDIS binning plots: resolution, migration matrix, purity/stability/efficiency
// g++ DDIS_Plot_Binning.cpp -o DDIS_Plot_Binning $(root-config --cflags --glibs)
// ./DDIS_Plot_Binning DDIS_BinningInputs.root [binning_scheme.txt] [output.root]

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TString.h>
#include <TLegend.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <unordered_map>

#include "BinningUtility.hpp"

struct BinDef {
    int bin_id = -1;
    double Q2_min = -1.0;
    double Q2_max = -1.0;
    double beta_min = -1.0;
    double beta_max = -1.0;
    double xpom_min = -1.0;
    double xpom_max = -1.0;
};

static bool StartsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static bool EndsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string TrimWS(const std::string& s) {
    const char* ws = " \t\r\n";
    const size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

static std::vector<BinDef> ReadBinsFromYAML(const std::string& path) {
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

static std::vector<double> UniqueSorted(std::vector<double> vals) {
    std::sort(vals.begin(), vals.end());
    const double eps = 1e-12;
    vals.erase(std::unique(vals.begin(), vals.end(), [eps](double a, double b) {
        return std::fabs(a - b) < eps;
    }), vals.end());
    return vals;
}

static void CollectEdges(const std::vector<BinDef>& bins,
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

static int FindEdgeIndex(const std::vector<double>& edges, double value) {
    const double eps = 1e-12;
    for(size_t i = 0; i + 1 < edges.size(); i++) {
        if(std::fabs(edges[i] - value) < eps) return static_cast<int>(i);
    }
    return -1;
}

static int FindBinId(const std::vector<BinDef>& bins,
                     double Q2,
                     double beta,
                     double xpom) {
    const double eps = 1e-12;
    for(const auto& b : bins) {
        if(Q2 >= b.Q2_min - eps && Q2 < b.Q2_max - eps &&
           beta >= b.beta_min - eps && beta < b.beta_max - eps &&
           xpom >= b.xpom_min - eps && xpom < b.xpom_max - eps) {
            return b.bin_id;
        }
    }
    return -1;
}

static void ValidateBinIds(const BinningScheme& binning, const std::vector<BinDef>& bins) {
    int mismatches = 0;
    for(const auto& b : bins) {
        if(b.bin_id < 0) continue;
        const int iQ2 = FindEdgeIndex(binning.Q2_edges, b.Q2_min);
        const int iBeta = FindEdgeIndex(binning.beta_edges, b.beta_min);
        const int iXpom = FindEdgeIndex(binning.xpom_edges, b.xpom_min);
        if(iQ2 < 0 || iBeta < 0 || iXpom < 0) continue;
        const int k = binning.GetLinearBin(iQ2, iBeta, iXpom);
        if(k != b.bin_id) {
            if(mismatches < 10) {
                std::cerr << "WARNING: bin_id mismatch (yaml=" << b.bin_id
                          << ", computed=" << k << ") for Q2_min=" << b.Q2_min
                          << ", xpom_min=" << b.xpom_min
                          << ", beta_min=" << b.beta_min << std::endl;
            }
            mismatches++;
        }
    }
    if(mismatches > 0) {
        std::cerr << "WARNING: total bin_id mismatches: " << mismatches << std::endl;
    }
}

struct MethodHists {
    std::string tag;
    TH2D* h_res_beta;
    TH2D* h_res_xpom;
    TH2D* h_res_beta_vs_xpom;
};

struct CombinedHists {
    TH2D* h_res_Q2;
    TH2D* h_migration;
    TH2D* h_response;
    TH1D* h_purity;
    TH1D* h_stability;
    TH1D* h_efficiency;
    TH1D* h_purity_Q2;
    TH1D* h_stability_Q2;
    TH1D* h_efficiency_Q2;
    TH2D* h_purity_beta_xpom;
    TH2D* h_efficiency_beta_xpom;
    std::vector<double> truth_count;
    std::vector<double> truth_with_reco;
    std::vector<double> reco_total;
    std::vector<double> diag_count;
};

struct BinIdMap {
    std::vector<int> bin_ids;
    std::unordered_map<int, int> id_to_index;
};

static BinIdMap BuildBinIdMap(const std::vector<BinDef>& bins) {
    BinIdMap map;
    map.bin_ids.reserve(bins.size());
    for(const auto& b : bins) {
        if(b.bin_id >= 0) {
            map.bin_ids.push_back(b.bin_id);
        }
    }
    std::sort(map.bin_ids.begin(), map.bin_ids.end());
    map.bin_ids.erase(std::unique(map.bin_ids.begin(), map.bin_ids.end()), map.bin_ids.end());
    for(size_t i = 0; i < map.bin_ids.size(); i++) {
        map.id_to_index[map.bin_ids[i]] = static_cast<int>(i);
    }
    return map;
}

static int MapBinIdToIndex(const BinIdMap& map, int bin_id) {
    auto it = map.id_to_index.find(bin_id);
    if(it == map.id_to_index.end()) return -1;
    return it->second;
}

static void ApplyBinIdLabels(TH1* h, const std::vector<int>& bin_ids) {
    if(!h || bin_ids.empty()) return;
    TAxis* ax = h->GetXaxis();
    const int nBins = ax->GetNbins();
    const int nLabels = std::min(static_cast<int>(bin_ids.size()), nBins);
    for(int i = 0; i < nLabels; i++) {
        ax->SetBinLabel(i + 1, Form("%d", bin_ids[i]));
    }
    ax->LabelsOption("v");
    ax->SetLabelSize(0.025);
}

static void ApplyBinIdLabels2D(TH2* h, const std::vector<int>& bin_ids) {
    if(!h || bin_ids.empty()) return;
    TAxis* ax = h->GetXaxis();
    TAxis* ay = h->GetYaxis();
    const int nBins = ax->GetNbins();
    const int nLabels = std::min(static_cast<int>(bin_ids.size()), nBins);
    for(int i = 0; i < nLabels; i++) {
        const TString label = Form("%d", bin_ids[i]);
        ax->SetBinLabel(i + 1, label);
        ay->SetBinLabel(i + 1, label);
    }
    ax->LabelsOption("v");
    ay->LabelsOption("h");
    ax->SetLabelSize(0.025);
    ay->SetLabelSize(0.025);
}

struct ResAccum {
    std::vector<double> sum;
    std::vector<double> sumsq;
    std::vector<int> count;

    explicit ResAccum(int nBins = 0)
        : sum(nBins, 0.0), sumsq(nBins, 0.0), count(nBins, 0) {}

    void Resize(int nBins) {
        sum.assign(nBins, 0.0);
        sumsq.assign(nBins, 0.0);
        count.assign(nBins, 0);
    }

    void Fill(int k, double value) {
        if(k < 0 || k >= static_cast<int>(sum.size())) return;
        sum[k] += value;
        sumsq[k] += value * value;
        count[k] += 1;
    }

    double Mean(int k) const {
        if(k < 0 || k >= static_cast<int>(sum.size()) || count[k] == 0) return 0.0;
        return sum[k] / static_cast<double>(count[k]);
    }

    double RMS(int k) const {
        if(k < 0 || k >= static_cast<int>(sum.size()) || count[k] == 0) return 0.0;
        double mean = Mean(k);
        double var = sumsq[k] / static_cast<double>(count[k]) - mean * mean;
        return (var > 0.0) ? std::sqrt(var) : 0.0;
    }

    double RMSError(int k) const {
        if(k < 0 || k >= static_cast<int>(sum.size()) || count[k] < 2) return 0.0;
        double rms = RMS(k);
        return rms / std::sqrt(2.0 * (count[k] - 1.0));
    }
};

static TH1D* BuildResVsKHist(const std::string& name,
                             const std::string& title,
                             const ResAccum& acc) {
    const int nBins = static_cast<int>(acc.sum.size());
    TH1D* h = new TH1D(name.c_str(), title.c_str(), nBins, 0.5, nBins + 0.5);
    for(int k = 0; k < nBins; k++) {
        if(acc.count[k] <= 0) continue;
        h->SetBinContent(k + 1, acc.RMS(k));
        h->SetBinError(k + 1, acc.RMSError(k));
    }
    return h;
}

static void WriteResTable(const std::string& filename,
                          const BinningScheme& binning,
                          const ResAccum& acc,
                          const std::string& tag,
                          const std::vector<BinDef>* bins = nullptr,
                          const BinIdMap* bin_map = nullptr) {
    std::ofstream out(filename);
    out << "# tag=" << tag << "\n";
    out << "# k iQ2 iXpom iBeta Q2_min Q2_max xpom_min xpom_max beta_min beta_max count mean rms\n";
    if(bins && bin_map) {
        for(const auto& b : *bins) {
            const int idx = MapBinIdToIndex(*bin_map, b.bin_id);
            if(idx < 0 || idx >= static_cast<int>(acc.sum.size())) continue;
            const int iQ2 = FindEdgeIndex(binning.Q2_edges, b.Q2_min);
            const int iXpom = FindEdgeIndex(binning.xpom_edges, b.xpom_min);
            const int iBeta = FindEdgeIndex(binning.beta_edges, b.beta_min);
            out << b.bin_id
                << " " << iQ2
                << " " << iXpom
                << " " << iBeta
                << " " << b.Q2_min
                << " " << b.Q2_max
                << " " << b.xpom_min
                << " " << b.xpom_max
                << " " << b.beta_min
                << " " << b.beta_max
                << " " << acc.count[idx]
                << " " << acc.Mean(idx)
                << " " << acc.RMS(idx)
                << "\n";
        }
    } else {
        for(int iQ2 = 0; iQ2 < binning.nQ2(); iQ2++) {
            for(int iXpom = 0; iXpom < binning.nXpom(); iXpom++) {
                for(int iBeta = 0; iBeta < binning.nBeta(); iBeta++) {
                    const int k = binning.GetLinearBin(iQ2, iBeta, iXpom);
                    const double Q2_min = binning.Q2_edges[iQ2];
                    const double Q2_max = binning.Q2_edges[iQ2 + 1];
                    const double xpom_min = binning.xpom_edges[iXpom];
                    const double xpom_max = binning.xpom_edges[iXpom + 1];
                    const double beta_min = binning.beta_edges[iBeta];
                    const double beta_max = binning.beta_edges[iBeta + 1];
                    out << k
                        << " " << iQ2
                        << " " << iXpom
                        << " " << iBeta
                        << " " << Q2_min
                        << " " << Q2_max
                        << " " << xpom_min
                        << " " << xpom_max
                        << " " << beta_min
                        << " " << beta_max
                        << " " << acc.count[k]
                        << " " << acc.Mean(k)
                        << " " << acc.RMS(k)
                        << "\n";
                }
            }
        }
    }
}

static void SaveResVsKPlot(const std::string& title,
                           const std::vector<TH1D*>& hists,
                           const std::vector<std::string>& labels,
                           const std::string& outpath,
                           const std::vector<int>* bin_ids = nullptr) {
    if(hists.empty()) return;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas c("c_res_vs_k", "c_res_vs_k", 1100, 600);
    c.SetGridx();
    c.SetGridy();

    double ymax = 0.0;
    for(const auto* h : hists) {
        if(!h) continue;
        ymax = std::max(ymax, h->GetMaximum());
    }
    if(ymax <= 0.0) ymax = 0.1;

    for(size_t i = 0; i < hists.size(); i++) {
        TH1D* h = hists[i];
        if(!h) continue;
        h->SetTitle(title.c_str());
        h->SetMarkerStyle(20 + static_cast<int>(i));
        h->SetMarkerSize(0.9);
        h->SetLineWidth(2);
        h->GetYaxis()->SetRangeUser(0.0, 1.2 * ymax);
        if(i == 0) {
            if(bin_ids && !bin_ids->empty()) {
                ApplyBinIdLabels(h, *bin_ids);
            }
            h->Draw("E1");
        } else {
            h->Draw("E1 SAME");
        }
    }

    if(labels.size() == hists.size()) {
        TLegend leg(0.7, 0.78, 0.92, 0.92);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        for(size_t i = 0; i < hists.size(); i++) {
            if(!hists[i]) continue;
            leg.AddEntry(hists[i], labels[i].c_str(), "lep");
        }
        leg.Draw();
    }

    c.SaveAs(outpath.c_str());
}

static MethodHists BookHists(const std::string& tag, const BinningScheme& binning) {
    MethodHists h;
    h.tag = tag;

    h.h_res_beta = new TH2D(("res_beta_" + tag).c_str(),
                            ("#beta resolution vs truth (" + tag + ");#beta_{truth};(#beta_{reco}-#beta_{truth})/#beta_{truth}").c_str(),
                            binning.nBeta(), binning.beta_edges.data(),
                            100, -1.0, 1.0);

    h.h_res_xpom = new TH2D(("res_xpom_" + tag).c_str(),
                            ("x_{pom} resolution vs truth (" + tag + ");x_{pom,truth};(x_{pom,reco}-x_{pom,truth})/x_{pom,truth}").c_str(),
                            binning.nXpom(), binning.xpom_edges.data(),
                            100, -1.0, 1.0);

    h.h_res_beta_vs_xpom = new TH2D(("res_beta_vs_xpom_" + tag).c_str(),
                                    ("Resolution correlation (" + tag + ");Res x_{pom};Res #beta").c_str(),
                                    100, -1.0, 1.0, 100, -1.0, 1.0);

    return h;
}

static CombinedHists BookCombinedHists(const BinningScheme& binning, int nBinsOverride = -1) {
    CombinedHists h;
    const int nBins = (nBinsOverride > 0) ? nBinsOverride : binning.nTotalBins();

    h.h_res_Q2 = new TH2D("res_Q2_all",
                          "Q^{2} resolution vs truth;Q^{2}_{truth};(Q^{2}_{reco}-Q^{2}_{truth})/Q^{2}_{truth}",
                          binning.nQ2(), binning.Q2_edges.data(),
                          100, -1.0, 1.0);

    h.h_migration = new TH2D("migration_all",
                             "Migration matrix;k_{true};k_{reco}",
                             nBins, 0.5, nBins + 0.5, nBins, 0.5, nBins + 0.5);
    h.h_response = new TH2D("response_matrix",
                            "Response matrix (%);k_{true};k_{reco}",
                            nBins, 0.5, nBins + 0.5, nBins, 0.5, nBins + 0.5);

    h.h_purity = new TH1D("purity_all",
                          "Purity per 3D bin;k_{bin};Purity",
                          nBins, 0.5, nBins + 0.5);

    h.h_stability = new TH1D("stability_all",
                             "Stability per 3D bin;k_{bin};Stability",
                             nBins, 0.5, nBins + 0.5);

    h.h_efficiency = new TH1D("efficiency_all",
                              "Efficiency per 3D bin;k_{bin};Efficiency",
                              nBins, 0.5, nBins + 0.5);

    h.h_purity_Q2 = new TH1D("purity_vs_Q2",
                              "Purity vs Q^{2};Q^{2} [GeV^{2}];Purity",
                              binning.nQ2(), binning.Q2_edges.data());
    h.h_stability_Q2 = new TH1D("stability_vs_Q2",
                                 "Stability vs Q^{2};Q^{2} [GeV^{2}];Stability",
                                 binning.nQ2(), binning.Q2_edges.data());
    h.h_efficiency_Q2 = new TH1D("efficiency_vs_Q2",
                                  "Efficiency vs Q^{2};Q^{2} [GeV^{2}];Efficiency",
                                  binning.nQ2(), binning.Q2_edges.data());
    h.h_purity_beta_xpom = new TH2D("purity_beta_xpom",
                                    "Purity vs (#beta, x_{pom});x_{pom};#beta",
                                    binning.nXpom(), binning.xpom_edges.data(),
                                    binning.nBeta(), binning.beta_edges.data());
    h.h_efficiency_beta_xpom = new TH2D("efficiency_beta_xpom",
                                        "Efficiency vs (#beta, x_{pom});x_{pom};#beta",
                                        binning.nXpom(), binning.xpom_edges.data(),
                                        binning.nBeta(), binning.beta_edges.data());

    h.truth_count.assign(nBins, 0.0);
    h.truth_with_reco.assign(nBins, 0.0);
    h.reco_total.assign(nBins, 0.0);
    h.diag_count.assign(nBins, 0.0);

    return h;
}

static void FinalizeMetrics(CombinedHists& h) {
    const int nBins = h.h_purity->GetNbinsX();
    for(int k = 0; k < nBins; k++) {
        double truth = h.truth_count[k];
        double truth_reco = h.truth_with_reco[k];
        double reco = h.reco_total[k];
        double diag = h.diag_count[k];

        double purity = (reco > 0.0) ? (diag / reco) : 0.0;
        double stability = (truth > 0.0) ? (diag / truth) : 0.0;
        double efficiency = (truth > 0.0) ? (truth_reco / truth) : 0.0;

        h.h_purity->SetBinContent(k + 1, purity);
        h.h_stability->SetBinContent(k + 1, stability);
        h.h_efficiency->SetBinContent(k + 1, efficiency);
    }
}

static void FinalizeResponseMatrix(CombinedHists& h) {
    const int nBins = h.h_migration->GetNbinsX();
    for(int i = 1; i <= nBins; i++) {
        double truth = h.truth_count[i - 1];
        for(int j = 1; j <= nBins; j++) {
            double count = h.h_migration->GetBinContent(i, j);
            double percent = (truth > 0.0) ? (100.0 * count / truth) : 0.0;
            h.h_response->SetBinContent(i, j, percent);
        }
    }
}

static void FinalizeQ2Metrics(CombinedHists& h, const BinningScheme& binning) {
    std::vector<double> truth_Q2(binning.nQ2(), 0.0);
    std::vector<double> truth_reco_Q2(binning.nQ2(), 0.0);
    std::vector<double> reco_Q2(binning.nQ2(), 0.0);
    std::vector<double> diag_Q2(binning.nQ2(), 0.0);

    const int nBins = binning.nTotalBins();
    for(int k = 0; k < nBins; k++) {
        int iQ2 = -1, iBeta = -1, iXpom = -1;
        binning.Get3DIndices(k, iQ2, iBeta, iXpom);
        if(iQ2 < 0 || iQ2 >= binning.nQ2()) continue;
        truth_Q2[iQ2] += h.truth_count[k];
        truth_reco_Q2[iQ2] += h.truth_with_reco[k];
        reco_Q2[iQ2] += h.reco_total[k];
        diag_Q2[iQ2] += h.diag_count[k];
    }

    for(int iQ2 = 0; iQ2 < binning.nQ2(); iQ2++) {
        double purity = (reco_Q2[iQ2] > 0.0) ? (diag_Q2[iQ2] / reco_Q2[iQ2]) : 0.0;
        double stability = (truth_Q2[iQ2] > 0.0) ? (diag_Q2[iQ2] / truth_Q2[iQ2]) : 0.0;
        double efficiency = (truth_Q2[iQ2] > 0.0) ? (truth_reco_Q2[iQ2] / truth_Q2[iQ2]) : 0.0;

        h.h_purity_Q2->SetBinContent(iQ2 + 1, purity);
        h.h_stability_Q2->SetBinContent(iQ2 + 1, stability);
        h.h_efficiency_Q2->SetBinContent(iQ2 + 1, efficiency);
    }
}

static void FinalizeBetaXpomMetrics(CombinedHists& h, const BinningScheme& binning) {
    const int nBeta = binning.nBeta();
    const int nXpom = binning.nXpom();
    std::vector<double> truth(nBeta * nXpom, 0.0);
    std::vector<double> truth_reco(nBeta * nXpom, 0.0);
    std::vector<double> reco(nBeta * nXpom, 0.0);
    std::vector<double> diag(nBeta * nXpom, 0.0);

    const int nBins = binning.nTotalBins();
    for(int k = 0; k < nBins; k++) {
        int iQ2 = -1, iBeta = -1, iXpom = -1;
        binning.Get3DIndices(k, iQ2, iBeta, iXpom);
        if(iBeta < 0 || iBeta >= nBeta || iXpom < 0 || iXpom >= nXpom) continue;
        int idx = iBeta * nXpom + iXpom;
        truth[idx] += h.truth_count[k];
        truth_reco[idx] += h.truth_with_reco[k];
        reco[idx] += h.reco_total[k];
        diag[idx] += h.diag_count[k];
    }

    for(int iBeta = 0; iBeta < nBeta; iBeta++) {
        for(int iXpom = 0; iXpom < nXpom; iXpom++) {
            int idx = iBeta * nXpom + iXpom;
            double purity = (reco[idx] > 0.0) ? (diag[idx] / reco[idx]) : 0.0;
            double efficiency = (truth[idx] > 0.0) ? (truth_reco[idx] / truth[idx]) : 0.0;
            h.h_purity_beta_xpom->SetBinContent(iXpom + 1, iBeta + 1, purity * 100.0);
            h.h_efficiency_beta_xpom->SetBinContent(iXpom + 1, iBeta + 1, efficiency * 100.0);
        }
    }
}

static void SavePlots(const MethodHists& h, const std::string& outdir, const std::string& jpgdir) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas c1("c1", "c1", 900, 700);

    c1.SetLogx(false);
    h.h_res_beta->Draw("COLZ");
    c1.SaveAs((outdir + "/res_beta_" + h.tag + ".png").c_str());
    c1.SaveAs((jpgdir + "/res_beta_" + h.tag + ".jpg").c_str());

    c1.SetLogx();
    h.h_res_xpom->Draw("COLZ");
    c1.SaveAs((outdir + "/res_xpom_" + h.tag + ".png").c_str());
    c1.SaveAs((jpgdir + "/res_xpom_" + h.tag + ".jpg").c_str());

    c1.SetLogx(false);
    h.h_res_beta_vs_xpom->Draw("COLZ");
    c1.SaveAs((outdir + "/res_beta_vs_xpom_" + h.tag + ".png").c_str());
    c1.SaveAs((jpgdir + "/res_beta_vs_xpom_" + h.tag + ".jpg").c_str());
}

static void SaveCombinedPlots(const CombinedHists& h, const std::string& outdir, const std::string& jpgdir) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas c1("c1_combined", "c1_combined", 900, 700);

    c1.SetLogx();
    h.h_res_Q2->Draw("COLZ");
    c1.SaveAs((outdir + "/res_Q2.png").c_str());
    c1.SaveAs((jpgdir + "/res_Q2.jpg").c_str());

    c1.SetLogx(false);
    h.h_migration->Draw("COLZ");
    c1.SaveAs((outdir + "/migration.png").c_str());
    c1.SaveAs((jpgdir + "/migration.jpg").c_str());

    h.h_purity->SetMinimum(0.0);
    h.h_purity->SetMaximum(1.0);
    h.h_purity->Draw("HIST");
    c1.SaveAs((outdir + "/purity.png").c_str());
    c1.SaveAs((jpgdir + "/purity.jpg").c_str());

    h.h_stability->SetMinimum(0.0);
    h.h_stability->SetMaximum(1.0);
    h.h_stability->Draw("HIST");
    c1.SaveAs((outdir + "/stability.png").c_str());
    c1.SaveAs((jpgdir + "/stability.jpg").c_str());

    h.h_efficiency->SetMinimum(0.0);
    h.h_efficiency->SetMaximum(1.0);
    h.h_efficiency->Draw("HIST");
    c1.SaveAs((outdir + "/efficiency.png").c_str());
    c1.SaveAs((jpgdir + "/efficiency.jpg").c_str());

    h.h_response->SetMinimum(0.0);
    h.h_response->SetMaximum(100.0);
    h.h_response->Draw("COLZ TEXT");
    {
        double min = 0.5;
        double max = h.h_response->GetNbinsX() + 0.5;
        TLine* diag = new TLine(min, min, max, max);
        diag->SetLineColor(kBlue + 1);
        diag->SetLineStyle(2);
        diag->SetLineWidth(2);
        diag->Draw("same");
    }
    c1.SaveAs((outdir + "/response_matrix.png").c_str());
    c1.SaveAs((jpgdir + "/response_matrix.jpg").c_str());

    h.h_purity_Q2->SetMinimum(0.0);
    h.h_purity_Q2->SetMaximum(1.0);
    h.h_purity_Q2->SetMarkerStyle(20);
    h.h_purity_Q2->Draw("E1");
    c1.SaveAs((outdir + "/purity_vs_Q2.png").c_str());
    c1.SaveAs((jpgdir + "/purity_vs_Q2.jpg").c_str());

    h.h_stability_Q2->SetMinimum(0.0);
    h.h_stability_Q2->SetMaximum(1.0);
    h.h_stability_Q2->SetMarkerStyle(20);
    h.h_stability_Q2->Draw("E1");
    c1.SaveAs((outdir + "/stability_vs_Q2.png").c_str());
    c1.SaveAs((jpgdir + "/stability_vs_Q2.jpg").c_str());

    h.h_efficiency_Q2->SetMinimum(0.0);
    h.h_efficiency_Q2->SetMaximum(1.0);
    h.h_efficiency_Q2->SetMarkerStyle(20);
    h.h_efficiency_Q2->Draw("E1");
    c1.SaveAs((outdir + "/efficiency_vs_Q2.png").c_str());
    c1.SaveAs((jpgdir + "/efficiency_vs_Q2.jpg").c_str());

    c1.SetLogx();
    h.h_purity_beta_xpom->SetMinimum(0.0);
    h.h_purity_beta_xpom->SetMaximum(100.0);
    h.h_purity_beta_xpom->Draw("COLZ TEXT");
    c1.SaveAs((outdir + "/purity_beta_xpom.png").c_str());
    c1.SaveAs((jpgdir + "/purity_beta_xpom.jpg").c_str());

    h.h_efficiency_beta_xpom->SetMinimum(0.0);
    h.h_efficiency_beta_xpom->SetMaximum(100.0);
    h.h_efficiency_beta_xpom->Draw("COLZ TEXT");
    c1.SaveAs((outdir + "/efficiency_beta_xpom.png").c_str());
    c1.SaveAs((jpgdir + "/efficiency_beta_xpom.jpg").c_str());
}

static int ClampBin(int bin, int nBins) {
    if(bin < 1) return 1;
    if(bin > nBins) return nBins;
    return bin;
}

static int FindBinLow(const TAxis* axis, double edge) {
    return ClampBin(axis->FindBin(edge + 1e-9), axis->GetNbins());
}

static int FindBinHigh(const TAxis* axis, double edge) {
    return ClampBin(axis->FindBin(edge - 1e-9), axis->GetNbins());
}

static void DrawConditionalGrid(TH2D* h, const std::vector<double>& x_edges,
                                const std::vector<double>& y_edges,
                                double minCount) {
    TAxis* xaxis = h->GetXaxis();
    TAxis* yaxis = h->GetYaxis();

    for(size_t ix = 0; ix + 1 < x_edges.size(); ix++) {
        double x1 = x_edges[ix];
        double x2 = x_edges[ix + 1];
        int bx1 = FindBinLow(xaxis, x1);
        int bx2 = FindBinHigh(xaxis, x2);

        for(size_t iy = 0; iy + 1 < y_edges.size(); iy++) {
            double y1 = y_edges[iy];
            double y2 = y_edges[iy + 1];
            int by1 = FindBinLow(yaxis, y1);
            int by2 = FindBinHigh(yaxis, y2);

            double count = h->Integral(bx1, bx2, by1, by2);
            if(count < minCount) continue;

            TBox* box = new TBox(x1, y1, x2, y2);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->Draw("same");
        }
    }
}

static void WriteConditionalBins(std::ofstream& out, const std::string& header,
                                 TH2D* h, const std::vector<double>& x_edges,
                                 const std::vector<double>& y_edges,
                                 const std::vector<double>* x_edges_phys,
                                 const std::vector<double>* y_edges_phys,
                                 double minCount) {
    out << "\n# " << header << "\n";
    out << "# x_low x_high y_low y_high count\n";
    TAxis* xaxis = h->GetXaxis();
    TAxis* yaxis = h->GetYaxis();
    out << std::setprecision(8);

    for(size_t ix = 0; ix + 1 < x_edges.size(); ix++) {
        double x1 = x_edges[ix];
        double x2 = x_edges[ix + 1];
        int bx1 = FindBinLow(xaxis, x1);
        int bx2 = FindBinHigh(xaxis, x2);

        for(size_t iy = 0; iy + 1 < y_edges.size(); iy++) {
            double y1 = y_edges[iy];
            double y2 = y_edges[iy + 1];
            int by1 = FindBinLow(yaxis, y1);
            int by2 = FindBinHigh(yaxis, y2);

            double count = h->Integral(bx1, bx2, by1, by2);
            if(count < minCount) continue;

            double x_low = x_edges_phys ? (*x_edges_phys)[ix] : x1;
            double x_high = x_edges_phys ? (*x_edges_phys)[ix + 1] : x2;
            double y_low = y_edges_phys ? (*y_edges_phys)[iy] : y1;
            double y_high = y_edges_phys ? (*y_edges_phys)[iy + 1] : y2;

            out << x_low << " " << x_high << " "
                << y_low << " " << y_high << " "
                << count << "\n";
        }
    }
}

static void SaveBinningSchemePlots(TTree* tree, const BinningScheme& binning,
                                   const std::string& outdir, const std::string& jpgdir) {
    float Q2_truth = -999.0f, beta_truth = -999.0f, xpom_truth = -999.0f, y_truth = -999.0f;
    tree->SetBranchAddress("Q2_truth", &Q2_truth);
    tree->SetBranchAddress("beta_truth", &beta_truth);
    tree->SetBranchAddress("xpom_truth", &xpom_truth);
    tree->SetBranchAddress("y_truth", &y_truth);

    const int nEntries = tree->GetEntries();
    TH2D* h_Q2_xpom = new TH2D("binning_Q2_vs_xpom",
                              "Q^{2} vs x_{pom};log_{10}(x_{pom});Q^{2} [GeV^{2}]",
                              120, TMath::Log10(binning.xpom_edges.front()),
                              TMath::Log10(binning.xpom_edges.back()),
                              120, binning.Q2_edges.front(), binning.Q2_edges.back());
    TH2D* h_Q2_beta = new TH2D("binning_Q2_vs_beta",
                              "Q^{2} vs #beta;#beta;Q^{2} [GeV^{2}]",
                              120, 0.0, 1.0,
                              120, binning.Q2_edges.front(), binning.Q2_edges.back());
    TH2D* h_beta_xpom = new TH2D("binning_beta_vs_xpom",
                                "#beta vs x_{pom};log_{10}(x_{pom});#beta",
                                120, TMath::Log10(binning.xpom_edges.front()),
                                TMath::Log10(binning.xpom_edges.back()),
                                120, 0.0, 1.0);
    TH2D* h_Q2_y = new TH2D("binning_Q2_vs_y",
                           "Q^{2} vs y;y;Q^{2} [GeV^{2}]",
                           120, 0.0, 1.0,
                           120, binning.Q2_edges.front(), binning.Q2_edges.back());

    for(int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if(Q2_truth > 0 && xpom_truth > 0) {
            h_Q2_xpom->Fill(TMath::Log10(xpom_truth), Q2_truth);
        }
        if(Q2_truth > 0 && beta_truth > 0) {
            h_Q2_beta->Fill(beta_truth, Q2_truth);
        }
        if(beta_truth > 0 && xpom_truth > 0) {
            h_beta_xpom->Fill(TMath::Log10(xpom_truth), beta_truth);
        }
        if(Q2_truth > 0 && y_truth > 0) {
            h_Q2_y->Fill(y_truth, Q2_truth);
        }
    }

    const double minCount = 100.0;

    std::vector<double> xpom_log_edges;
    xpom_log_edges.reserve(binning.xpom_edges.size());
    for(double edge : binning.xpom_edges) {
        xpom_log_edges.push_back(TMath::Log10(edge));
    }

    {
        std::ofstream out("data/binning_scheme_cells.txt", std::ios::trunc);
        if(out.is_open()) {
            out << "# Bins with counts >= " << minCount << "\n";
            WriteConditionalBins(out, "Q2 vs x_pom (x_pom in linear units)",
                                 h_Q2_xpom, xpom_log_edges, binning.Q2_edges,
                                 &binning.xpom_edges, &binning.Q2_edges, minCount);
            WriteConditionalBins(out, "Q2 vs beta",
                                 h_Q2_beta, binning.beta_edges, binning.Q2_edges,
                                 &binning.beta_edges, &binning.Q2_edges, minCount);
            WriteConditionalBins(out, "beta vs x_pom (x_pom in linear units)",
                                 h_beta_xpom, xpom_log_edges, binning.beta_edges,
                                 &binning.xpom_edges, &binning.beta_edges, minCount);
        }
    }

    TCanvas c1("c1_scheme", "c1_scheme", 900, 700);
    c1.SetLogy();
    h_Q2_xpom->Draw("COLZ");
    DrawConditionalGrid(h_Q2_xpom, xpom_log_edges, binning.Q2_edges, minCount);
    c1.SaveAs((outdir + "/binning_Q2_vs_xpom.png").c_str());
    c1.SaveAs((jpgdir + "/binning_Q2_vs_xpom.jpg").c_str());

    c1.SetLogy();
    h_Q2_beta->Draw("COLZ");
    DrawConditionalGrid(h_Q2_beta, binning.beta_edges, binning.Q2_edges, minCount);
    c1.SaveAs((outdir + "/binning_Q2_vs_beta.png").c_str());
    c1.SaveAs((jpgdir + "/binning_Q2_vs_beta.jpg").c_str());

    c1.SetLogy(false);
    h_beta_xpom->Draw("COLZ");
    DrawConditionalGrid(h_beta_xpom, xpom_log_edges, binning.beta_edges, minCount);
    c1.SaveAs((outdir + "/binning_beta_vs_xpom.png").c_str());
    c1.SaveAs((jpgdir + "/binning_beta_vs_xpom.jpg").c_str());

    c1.SetLogy();
    h_Q2_y->Draw("COLZ");
    c1.SaveAs((outdir + "/binning_Q2_vs_y.png").c_str());
    c1.SaveAs((jpgdir + "/binning_Q2_vs_y.jpg").c_str());
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <DDIS_BinningInputs.root> [binning.yaml|binning_scheme.txt] [output.root]" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string binningFile = (argc > 2) ? argv[2] : "data/bins_template.yaml";
    std::string outputFile = (argc > 3) ? argv[3] : "DDIS_BinningPlots.root";

    BinningScheme binning;
    std::vector<BinDef> yaml_bins;
    const bool use_yaml = (EndsWith(binningFile, ".yaml") || EndsWith(binningFile, ".yml"));
    if (use_yaml) {
        yaml_bins = ReadBinsFromYAML(binningFile);
        if(yaml_bins.empty()) {
            std::cerr << "ERROR: No bins found in YAML: " << binningFile << std::endl;
            return 1;
        }
        CollectEdges(yaml_bins, binning.Q2_edges, binning.beta_edges, binning.xpom_edges);
        ValidateBinIds(binning, yaml_bins);
    } else {
        binning = ReadBinningScheme(binningFile);
    }
    if(binning.Q2_edges.empty() || binning.beta_edges.empty() || binning.xpom_edges.empty()) {
        std::cerr << "ERROR: Invalid binning scheme: " << binningFile << std::endl;
        return 1;
    }
    PrintBinningScheme(binning);

    TFile* infile = TFile::Open(inputFile.c_str(), "READ");
    if(!infile || infile->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file " << inputFile << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)infile->Get("binning");
    if(!tree) {
        std::cerr << "ERROR: Cannot find 'binning' tree in file" << std::endl;
        return 1;
    }

    int method = -1;
    float Q2_truth = -999.0f, Q2_reco = -999.0f;
    float x_truth = -999.0f, x_reco = -999.0f;
    float xpom_truth = -999.0f, xpom_reco = -999.0f;
    float beta_truth = -999.0f, beta_reco = -999.0f;
    float y_truth = -999.0f;

    tree->SetBranchAddress("method", &method);
    tree->SetBranchAddress("Q2_truth", &Q2_truth);
    tree->SetBranchAddress("Q2_reco", &Q2_reco);
    tree->SetBranchAddress("x_truth", &x_truth);
    tree->SetBranchAddress("x_reco", &x_reco);
    tree->SetBranchAddress("xpom_truth", &xpom_truth);
    tree->SetBranchAddress("xpom_reco", &xpom_reco);
    tree->SetBranchAddress("beta_truth", &beta_truth);
    tree->SetBranchAddress("beta_reco", &beta_reco);
    tree->SetBranchAddress("y_truth", &y_truth);

    MethodHists h_b0 = BookHists("B0", binning);
    MethodHists h_rp = BookHists("RP", binning);
    BinIdMap bin_map;
    if(use_yaml) {
        bin_map = BuildBinIdMap(yaml_bins);
        if(bin_map.bin_ids.empty()) {
            std::cerr << "ERROR: No valid bin_id entries in YAML." << std::endl;
            return 1;
        }
    }
    const int nBinsTotal = use_yaml ? static_cast<int>(bin_map.bin_ids.size()) : binning.nTotalBins();
    CombinedHists h_all = BookCombinedHists(binning, nBinsTotal);
    if(use_yaml) {
        ApplyBinIdLabels(h_all.h_purity, bin_map.bin_ids);
        ApplyBinIdLabels(h_all.h_stability, bin_map.bin_ids);
        ApplyBinIdLabels(h_all.h_efficiency, bin_map.bin_ids);
        ApplyBinIdLabels2D(h_all.h_migration, bin_map.bin_ids);
        ApplyBinIdLabels2D(h_all.h_response, bin_map.bin_ids);
    }
    ResAccum res_Q2(nBinsTotal);
    ResAccum res_beta_B0(nBinsTotal);
    ResAccum res_beta_RP(nBinsTotal);
    ResAccum res_xpom_B0(nBinsTotal);
    ResAccum res_xpom_RP(nBinsTotal);

    const Long64_t nEntries = tree->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        MethodHists* h = (method == 0) ? &h_b0 : &h_rp;

        int iQ2_true = -1;
        int iBeta_true = -1;
        int iXpom_true = -1;
        int iQ2_reco = -1;
        int iBeta_reco = -1;
        int iXpom_reco = -1;

        int k_true = -1;
        int k_reco = -1;
        if(Q2_truth > 0 && beta_truth > 0 && xpom_truth > 0) {
            if(use_yaml) {
                k_true = FindBinId(yaml_bins, Q2_truth, beta_truth, xpom_truth);
            } else {
                iQ2_true = binning.GetQ2Bin(Q2_truth);
                iBeta_true = binning.GetBetaBin(beta_truth);
                iXpom_true = binning.GetXpomBin(xpom_truth);
                if(iQ2_true >= 0 && iBeta_true >= 0 && iXpom_true >= 0) {
                    k_true = binning.GetLinearBin(iQ2_true, iBeta_true, iXpom_true);
                }
            }
        }
        if(Q2_reco > 0 && beta_reco > 0 && xpom_reco > 0) {
            if(use_yaml) {
                k_reco = FindBinId(yaml_bins, Q2_reco, beta_reco, xpom_reco);
            } else {
                iQ2_reco = binning.GetQ2Bin(Q2_reco);
                iBeta_reco = binning.GetBetaBin(beta_reco);
                iXpom_reco = binning.GetXpomBin(xpom_reco);
                if(iQ2_reco >= 0 && iBeta_reco >= 0 && iXpom_reco >= 0) {
                    k_reco = binning.GetLinearBin(iQ2_reco, iBeta_reco, iXpom_reco);
                }
            }
        }

        int idx_true = k_true;
        int idx_reco = k_reco;
        if(use_yaml) {
            idx_true = MapBinIdToIndex(bin_map, k_true);
            idx_reco = MapBinIdToIndex(bin_map, k_reco);
        }

        if(idx_true >= 0) {
            h_all.truth_count[idx_true] += 1.0;
        }
        if(idx_reco >= 0) {
            h_all.reco_total[idx_reco] += 1.0;
        }
        if(idx_true >= 0 && idx_reco >= 0) {
            h_all.truth_with_reco[idx_true] += 1.0;
            h_all.h_migration->Fill(idx_true + 1, idx_reco + 1);
            if(idx_true == idx_reco) {
                h_all.diag_count[idx_true] += 1.0;
            }
        }

        if(k_true >= 0 && k_reco >= 0) {
            if(Q2_truth != 0.0f) {
                double res_Q2_val = (Q2_reco - Q2_truth) / Q2_truth;
                h_all.h_res_Q2->Fill(Q2_truth, res_Q2_val);
                if(method == 0 && idx_true >= 0) {
                    res_Q2.Fill(idx_true, res_Q2_val);
                }
            }
            if(beta_truth != 0.0f) {
                double res_beta_val = (beta_reco - beta_truth) / beta_truth;
                h->h_res_beta->Fill(beta_truth, res_beta_val);
                if(idx_true >= 0) {
                    if(method == 0) {
                        res_beta_B0.Fill(idx_true, res_beta_val);
                    } else if(method == 1) {
                        res_beta_RP.Fill(idx_true, res_beta_val);
                    }
                }

                if(xpom_truth != 0.0f) {
                    double res_xpom_val = (xpom_reco - xpom_truth) / xpom_truth;
                    h->h_res_xpom->Fill(xpom_truth, res_xpom_val);
                    h->h_res_beta_vs_xpom->Fill(res_xpom_val, res_beta_val);
                    if(idx_true >= 0) {
                        if(method == 0) {
                            res_xpom_B0.Fill(idx_true, res_xpom_val);
                        } else if(method == 1) {
                            res_xpom_RP.Fill(idx_true, res_xpom_val);
                        }
                    }
                }
            }
        }
    }

    FinalizeMetrics(h_all);
    FinalizeQ2Metrics(h_all, binning);
    FinalizeBetaXpomMetrics(h_all, binning);
    FinalizeResponseMatrix(h_all);

    TFile* outfile = TFile::Open(outputFile.c_str(), "RECREATE");
    h_b0.h_res_beta->Write();
    h_b0.h_res_xpom->Write();
    h_b0.h_res_beta_vs_xpom->Write();
    h_all.h_res_Q2->Write();
    h_all.h_migration->Write();
    h_all.h_purity->Write();
    h_all.h_stability->Write();
    h_all.h_efficiency->Write();
    h_all.h_purity_Q2->Write();
    h_all.h_stability_Q2->Write();
    h_all.h_efficiency_Q2->Write();
    h_all.h_purity_beta_xpom->Write();
    h_all.h_efficiency_beta_xpom->Write();
    h_all.h_response->Write();

    h_rp.h_res_beta->Write();
    h_rp.h_res_xpom->Write();
    h_rp.h_res_beta_vs_xpom->Write();

    TH1D* h_res_Q2_vs_k = BuildResVsKHist(
        "res_Q2_vs_k",
        "Relative Q^{2} resolution vs global bin; k; RMS((Q^{2}_{reco}-Q^{2}_{truth})/Q^{2}_{truth})",
        res_Q2);
    h_res_Q2_vs_k->SetDirectory(nullptr);
    TH1D* h_res_beta_vs_k_B0 = BuildResVsKHist(
        "res_beta_vs_k_B0",
        "Relative #beta resolution vs global bin; k; RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_B0);
    h_res_beta_vs_k_B0->SetDirectory(nullptr);
    TH1D* h_res_beta_vs_k_RP = BuildResVsKHist(
        "res_beta_vs_k_RP",
        "Relative #beta resolution vs global bin; k; RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_RP);
    h_res_beta_vs_k_RP->SetDirectory(nullptr);
    TH1D* h_res_xpom_vs_k_B0 = BuildResVsKHist(
        "res_xpom_vs_k_B0",
        "Relative x_{pom} resolution vs global bin; k; RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_B0);
    h_res_xpom_vs_k_B0->SetDirectory(nullptr);
    TH1D* h_res_xpom_vs_k_RP = BuildResVsKHist(
        "res_xpom_vs_k_RP",
        "Relative x_{pom} resolution vs global bin; k; RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_RP);
    h_res_xpom_vs_k_RP->SetDirectory(nullptr);

    h_res_Q2_vs_k->SetMarkerColor(kBlack);
    h_res_Q2_vs_k->SetLineColor(kBlack);

    h_res_beta_vs_k_B0->SetMarkerColor(kRed + 1);
    h_res_beta_vs_k_B0->SetLineColor(kRed + 1);
    h_res_beta_vs_k_RP->SetMarkerColor(kBlue + 1);
    h_res_beta_vs_k_RP->SetLineColor(kBlue + 1);

    h_res_xpom_vs_k_B0->SetMarkerColor(kRed + 1);
    h_res_xpom_vs_k_B0->SetLineColor(kRed + 1);
    h_res_xpom_vs_k_RP->SetMarkerColor(kBlue + 1);
    h_res_xpom_vs_k_RP->SetLineColor(kBlue + 1);

    h_res_Q2_vs_k->Write();
    h_res_beta_vs_k_B0->Write();
    h_res_beta_vs_k_RP->Write();
    h_res_xpom_vs_k_B0->Write();
    h_res_xpom_vs_k_RP->Write();

    outfile->Close();

    const std::string outdir = "figs/binning";
    const std::string jpgdir = "figs/binning/jpg";
    gSystem->mkdir(outdir.c_str(), true);
    gSystem->mkdir(jpgdir.c_str(), true);
    SavePlots(h_b0, outdir, jpgdir);
    SavePlots(h_rp, outdir, jpgdir);
    SaveCombinedPlots(h_all, outdir, jpgdir);
    SaveBinningSchemePlots(tree, binning, outdir, jpgdir);

    if(use_yaml) {
        WriteResTable(outdir + "/res_vs_k_Q2.txt", binning, res_Q2, "Q2", &yaml_bins, &bin_map);
        WriteResTable(outdir + "/res_vs_k_beta_B0.txt", binning, res_beta_B0, "beta_B0", &yaml_bins, &bin_map);
        WriteResTable(outdir + "/res_vs_k_beta_RP.txt", binning, res_beta_RP, "beta_RP", &yaml_bins, &bin_map);
        WriteResTable(outdir + "/res_vs_k_xpom_B0.txt", binning, res_xpom_B0, "xpom_B0", &yaml_bins, &bin_map);
        WriteResTable(outdir + "/res_vs_k_xpom_RP.txt", binning, res_xpom_RP, "xpom_RP", &yaml_bins, &bin_map);
    } else {
        WriteResTable(outdir + "/res_vs_k_Q2.txt", binning, res_Q2, "Q2");
        WriteResTable(outdir + "/res_vs_k_beta_B0.txt", binning, res_beta_B0, "beta_B0");
        WriteResTable(outdir + "/res_vs_k_beta_RP.txt", binning, res_beta_RP, "beta_RP");
        WriteResTable(outdir + "/res_vs_k_xpom_B0.txt", binning, res_xpom_B0, "xpom_B0");
        WriteResTable(outdir + "/res_vs_k_xpom_RP.txt", binning, res_xpom_RP, "xpom_RP");
    }

    const std::vector<int>* label_ids = use_yaml ? &bin_map.bin_ids : nullptr;

    SaveResVsKPlot("Relative Q^{2} resolution vs global bin",
                   {h_res_Q2_vs_k},
                   {"EM"},
                   outdir + "/res_Q2_vs_k.png",
                   label_ids);

    SaveResVsKPlot("Relative #beta resolution vs global bin",
                   {h_res_beta_vs_k_B0, h_res_beta_vs_k_RP},
                   {"B0", "RP"},
                   outdir + "/res_beta_vs_k.png",
                   label_ids);

    SaveResVsKPlot("Relative x_{pom} resolution vs global bin",
                   {h_res_xpom_vs_k_B0, h_res_xpom_vs_k_RP},
                   {"B0", "RP"},
                   outdir + "/res_xpom_vs_k.png",
                   label_ids);

    infile->Close();
    std::cout << "Wrote plots to " << outdir << " and histograms to " << outputFile << std::endl;

    return 0;
}
