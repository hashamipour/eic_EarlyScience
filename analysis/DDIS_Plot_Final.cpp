// Inclusive DIS plotter
// g++ DDIS_Plot_Final.cpp -o DDIS_Plot_Final $(root-config --cflags --glibs)
// ./DDIS_Plot_Final <combined.root>

#include "Plotting.hpp"

#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TTree.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TString.h>
#include <TLine.h>
#include <TColor.h>
#include <TMarker.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TPad.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>

struct BinDef {
    int bin_id = -1;
    double Q2_min = -1.0;
    double Q2_max = -1.0;
    double beta_min = -1.0;
    double beta_max = -1.0;
    double xpom_min = -1.0;
    double xpom_max = -1.0;
};

static void DrawBinningGridWithCounts(TH2* hist,
                                      const std::vector<double>& xbins,
                                      const std::vector<double>& ybins,
                                      bool logX,
                                      bool logY);
static std::vector<double> BuildLinEdges(double minVal, double maxVal, int nSlices);

static bool StartsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
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

static void PlotXQ2Density(TFile* inputFile) {
    if (!inputFile) return;
    const int n_x_bins = 120;
    const int n_q2_bins = 120;
    const double q2_min = 1.0;
    const double q2_max = 1000.0;
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, n_x_bins);
    std::vector<Double_t> q2_bins = GetLogBins(q2_min, q2_max, n_q2_bins);

    auto makeHist = [&](const char* name, const char* title) {
        return new TH2D(name,
                        title,
                        x_bins.size() - 1, x_bins.data(),
                        q2_bins.size() - 1, q2_bins.data());
    };

    TH2D* h_reco_file = (TH2D*)inputFile->Get("xQ2_reco");
    TH2D* h_truth_file = (TH2D*)inputFile->Get("xQ2_truth");
    TH2D* h_reco_tmp = nullptr;
    TH2D* h_truth_tmp = nullptr;

    // Build missing reco/truth histograms from Q2_tree only when needed.
    if (!h_reco_file || !h_truth_file) {
        TTree* tree = (TTree*)inputFile->Get("Q2_tree");
        if (tree) {
            const bool hasRecoBranches = tree->GetBranch("x_EM") && tree->GetBranch("Q2_EM");
            const bool hasTruthBranches = tree->GetBranch("x_truth") && tree->GetBranch("Q2_truth");

            if (!h_reco_file && hasRecoBranches) {
                h_reco_tmp = makeHist("xQ2_reco_tmp", "Reco Event Density;x_{Bj};Q^{2} [GeV^{2}]");
            }
            if (!h_truth_file && hasTruthBranches) {
                h_truth_tmp = makeHist("xQ2_truth_tmp", "Truth Event Density;x_{Bj};Q^{2} [GeV^{2}]");
            }

            if (h_reco_tmp || h_truth_tmp) {
                float x_em = -999.0f, q2_em = -999.0f;
                float x_truth = -999.0f, q2_truth = -999.0f;
                if (h_reco_tmp) {
                    tree->SetBranchAddress("x_EM", &x_em);
                    tree->SetBranchAddress("Q2_EM", &q2_em);
                }
                if (h_truth_tmp) {
                    tree->SetBranchAddress("x_truth", &x_truth);
                    tree->SetBranchAddress("Q2_truth", &q2_truth);
                }

                const Long64_t nEntries = tree->GetEntries();
                for (Long64_t i = 0; i < nEntries; ++i) {
                    tree->GetEntry(i);
                    if (h_reco_tmp &&
                        std::isfinite(x_em) && x_em > 0.0f && x_em < 1.0f &&
                        std::isfinite(q2_em) && q2_em > 0.0f) {
                        h_reco_tmp->Fill(x_em, q2_em);
                    }
                    if (h_truth_tmp &&
                        std::isfinite(x_truth) && x_truth > 0.0f && x_truth < 1.0f &&
                        std::isfinite(q2_truth) && q2_truth > 0.0f) {
                        h_truth_tmp->Fill(x_truth, q2_truth);
                    }
                }
            }
        }
    }

    TH2* h_reco = h_reco_file ? static_cast<TH2*>(h_reco_file) : static_cast<TH2*>(h_reco_tmp);
    TH2* h_truth = h_truth_file ? static_cast<TH2*>(h_truth_file) : static_cast<TH2*>(h_truth_tmp);

    if (!h_reco && !h_truth) {
        std::cerr << "Warning: Could not find or build xQ2_reco/xQ2_truth; skipping x-Q2 density plots." << std::endl;
        delete h_reco_tmp;
        delete h_truth_tmp;
        return;
    }

    auto drawDensity = [&](TH2* h_src, const char* canvasName, const char* saveName) {
        if (!h_src) return;
        TH2* h = (TH2*)h_src->Clone(Form("%s_draw_%s", h_src->GetName(), canvasName));
        h->SetDirectory(nullptr);
        TCanvas* c = new TCanvas(canvasName, "x-Q2 Density", 1200, 1000);
        gStyle->SetOptStat(0);
        c->SetRightMargin(0.15);
        c->SetLeftMargin(0.15);
        c->SetTopMargin(0.1);
        c->SetBottomMargin(0.1);
        c->SetLogx();
        c->SetLogy();

        h->SetTitle("");
        h->GetXaxis()->SetTitle("x_{Bj}");
        h->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h->GetXaxis()->SetTitleOffset(1.3);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->Draw("COLZ");

        DrawSimLabels(inputFile);

        c->Update();
        SaveCanvas(c, saveName);

        delete h;
        delete c;
    };

    drawDensity(h_reco, "c_xq2_density_reco", "figs/inclusive/histos/xbj_q2_density_reco.png");
    drawDensity(h_truth, "c_xq2_density_truth", "figs/inclusive/histos/xbj_q2_density_truth.png");

    // Backward-compatible alias used by existing slides/docs.
    if (h_reco) {
        drawDensity(h_reco, "c_xq2_density_legacy", "figs/inclusive/histos/xbj_q2_density.png");
    } else if (h_truth) {
        drawDensity(h_truth, "c_xq2_density_legacy", "figs/inclusive/histos/xbj_q2_density.png");
    }

    delete h_reco_tmp;
    delete h_truth_tmp;
}

static void PlotDensityFromHist(TFile* inputFile,
                                const char* histName,
                                const char* xLabel,
                                const char* yLabel,
                                const char* saveName,
                                const bool logX,
                                const bool logY) {
    if (!inputFile) return;
    TH2* h = (TH2*)inputFile->Get(histName);
    if (!h) {
        if (strstr(histName, "_reco") != nullptr) {
            std::string fallback = histName;
            fallback.replace(fallback.find("_reco"), 5, "_truth");
            h = (TH2*)inputFile->Get(fallback.c_str());
            if (h) {
                std::cerr << "Warning: Histogram " << histName << " not found; using " << fallback << " instead."
                          << std::endl;
            }
        }
        if (!h) {
            std::cerr << "Warning: Histogram " << histName << " not found; skipping density plot." << std::endl;
            return;
        }
    }

    TCanvas* c = new TCanvas(Form("c_%s", histName), histName, 1200, 1000);
    gStyle->SetOptStat(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h->SetTitle("");
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->Draw("COLZ");

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
}

static void PlotDensityFromHistWithOverlay(TFile* inputFile,
                                           const char* histName,
                                           const char* xLabel,
                                           const char* yLabel,
                                           const char* saveName,
                                           const bool logX,
                                           const bool logY,
                                           const std::vector<double>& xbins,
                                           const std::vector<double>& ybins,
                                           const bool overlayLogX,
                                           const bool overlayLogY) {
    if (!inputFile) return;
    TH2* h = (TH2*)inputFile->Get(histName);
    if (!h) {
        std::cerr << "Warning: Histogram " << histName << " not found; skipping density plot." << std::endl;
        return;
    }
    if (xbins.size() < 2 || ybins.size() < 2) {
        std::cerr << "Warning: Overlay bins missing for " << histName << " (xbins="
                  << xbins.size() << ", ybins=" << ybins.size() << ")." << std::endl;
    }

    TCanvas* c = new TCanvas(Form("c_%s", histName), histName, 1200, 1000);
    gStyle->SetOptStat(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h->SetTitle("");
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->Draw("COLZ");

    DrawBinningGridWithCounts(h, xbins, ybins, overlayLogX, overlayLogY);

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
}

static void PlotGraphDensity(TFile* inputFile,
                             const char* graphName,
                             const char* title,
                             const char* xLabel,
                             const char* yLabel,
                             const char* saveName,
                             double xmin,
                             double xmax,
                             int nBins,
                             bool logX,
                             bool logY) {
    if (!inputFile) return;
    TGraph* g = (TGraph*)inputFile->Get(graphName);
    if (!g) {
        std::cerr << "Warning: Graph " << graphName << " not found; skipping density plot." << std::endl;
        return;
    }
    std::vector<Double_t> xbins = logX ? GetLogBins(xmin, xmax, nBins)
                                       : BuildLinEdges(xmin, xmax, nBins);
    std::vector<Double_t> ybins = logY ? GetLogBins(xmin, xmax, nBins)
                                       : BuildLinEdges(xmin, xmax, nBins);

    TH2D* h = new TH2D(Form("h_%s_density", graphName),
                       title,
                       xbins.size() - 1, xbins.data(),
                       ybins.size() - 1, ybins.data());

    const int n = g->GetN();
    double x = 0.0, y = 0.0;
    for (int i = 0; i < n; ++i) {
        g->GetPoint(i, x, y);
        if (!std::isfinite(x) || !std::isfinite(y)) continue;
        if (x <= 0.0 || y <= 0.0) continue;
        h->Fill(x, y);
    }

    TCanvas* c = new TCanvas(Form("c_%s_density", graphName), title, 1200, 1000);
    gStyle->SetOptStat(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h->SetTitle("");
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->Draw("COLZ");

    TLine* diag = new TLine(xmin, xmin, xmax, xmax);
    diag->SetLineColor(kBlack);
    diag->SetLineStyle(2);
    diag->SetLineWidth(2);
    diag->Draw("same");

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete diag;
    delete c;
    delete h;
}

static std::vector<double> BuildLogEdges(double minVal, double maxVal, int nSlices) {
    std::vector<double> edges;
    if (nSlices <= 0 || minVal <= 0.0 || maxVal <= 0.0 || maxVal <= minVal) {
        return edges;
    }
    edges.reserve(nSlices + 1);
    const double logMin = TMath::Log10(minVal);
    const double logMax = TMath::Log10(maxVal);
    const double step = (logMax - logMin) / nSlices;
    for (int i = 0; i <= nSlices; ++i) {
        edges.push_back(TMath::Power(10.0, logMin + step * i));
    }
    return edges;
}

static std::vector<double> BuildLinEdges(double minVal, double maxVal, int nSlices) {
    std::vector<double> edges;
    if (nSlices <= 0 || maxVal <= minVal) {
        return edges;
    }
    edges.reserve(nSlices + 1);
    const double step = (maxVal - minVal) / nSlices;
    for (int i = 0; i <= nSlices; ++i) {
        edges.push_back(minVal + step * i);
    }
    return edges;
}

static std::string FormatRange(double lo, double hi) {
    auto fmt = [](double v) {
        std::ostringstream oss;
        if (v >= 0.01 && v < 1000.0) {
            oss << std::fixed << std::setprecision(2) << v;
        } else {
            oss << std::scientific << std::setprecision(2) << v;
        }
        return oss.str();
    };
    std::ostringstream out;
    out << "[" << fmt(lo) << ", " << fmt(hi) << "]";
    return out.str();
}

static void SetGreenYellowRedPalette() {
    const Int_t nRGBs = 4;
    const Int_t nCont = 256;
    Double_t stops[nRGBs] = {0.0, 0.5, 0.75, 1.0};
    Double_t red[nRGBs]   = {0.18, 0.98, 0.98, 0.90};
    Double_t green[nRGBs] = {0.68, 0.90, 0.60, 0.20};
    Double_t blue[nRGBs]  = {0.38, 0.20, 0.20, 0.20};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
    gStyle->SetNumberContours(nCont);
}

// Adapted from analysis/Plot_BinningScheme_WithCounts.cpp
static void DrawBinningGridWithCounts(TH2* hist,
                                      const std::vector<double>& xbins,
                                      const std::vector<double>& ybins,
                                      bool logX,
                                      bool logY) {
    if (!hist || xbins.size() < 2 || ybins.size() < 2) {
        if (hist && (xbins.size() < 2 || ybins.size() < 2)) {
            std::cerr << "Warning: bin overlay skipped (empty bin edges)." << std::endl;
        }
        return;
    }

    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double ymin = hist->GetYaxis()->GetXmin();
    const double ymax = hist->GetYaxis()->GetXmax();

    const int lineColor = kRed;
    const int lineWidth = 3;
    const int lineStyle = 1;

    for (double xval : xbins) {
        if (logX && xval <= 0.0) continue;
        if (xval < xmin || xval > xmax) continue;
        TLine* line = new TLine(xval, ymin, xval, ymax);
        line->SetLineColor(lineColor);
        line->SetLineWidth(lineWidth);
        line->SetLineStyle(lineStyle);
        line->Draw("same");
    }

    for (double yval : ybins) {
        if (logY && yval <= 0.0) continue;
        if (yval < ymin || yval > ymax) continue;
        TLine* line = new TLine(xmin, yval, xmax, yval);
        line->SetLineColor(lineColor);
        line->SetLineWidth(lineWidth);
        line->SetLineStyle(lineStyle);
        line->Draw("same");
    }

    TLatex countText;
    countText.SetTextColor(kBlack);
    countText.SetTextSize(0.022);
    countText.SetTextAlign(22);
    countText.SetTextFont(42);

    TLatex idText;
    idText.SetTextColor(kBlue + 2);
    idText.SetTextSize(0.018);
    idText.SetTextAlign(22);
    idText.SetTextFont(42);

    const double tiny = 1e-12;
    const int nx = static_cast<int>(xbins.size()) - 1;
    for (size_t ix = 0; ix + 1 < xbins.size(); ++ix) {
        const double xlow = xbins[ix];
        const double xhigh = xbins[ix + 1];
        if (xhigh <= xmin || xlow >= xmax) continue;
        const double xlo = std::max(xlow, xmin);
        const double xhi = std::min(xhigh, xmax);
        const double xcenter = logX ? std::sqrt(xlo * xhi) : 0.5 * (xlo + xhi);
        const int binxlow = hist->GetXaxis()->FindBin(xlo + tiny);
        const int binxhigh = hist->GetXaxis()->FindBin(xhi - tiny);

        for (size_t iy = 0; iy + 1 < ybins.size(); ++iy) {
            const double ylow = ybins[iy];
            const double yhigh = ybins[iy + 1];
            if (yhigh <= ymin || ylow >= ymax) continue;
            const double ylo = std::max(ylow, ymin);
            const double yhi = std::min(yhigh, ymax);
            const double ycenter = logY ? std::sqrt(ylo * yhi) : 0.5 * (ylo + yhi);
            const int binylow = hist->GetYaxis()->FindBin(ylo + tiny);
            const int binyhigh = hist->GetYaxis()->FindBin(yhi - tiny);
            const int bin_id = static_cast<int>(ix) + static_cast<int>(iy) * nx;
            const double y_id = logY ? (ycenter * std::pow(yhi / ylo, 0.18))
                                     : (ycenter + 0.16 * (yhi - ylo));
            const double y_count = logY ? (ycenter / std::pow(yhi / ylo, 0.10))
                                        : (ycenter - 0.08 * (yhi - ylo));

            const double count = hist->Integral(binxlow, binxhigh, binylow, binyhigh);
            if (count > 0.0) {
                countText.DrawLatex(xcenter, y_count, Form("%d", (int)std::lround(count)));
            }
            idText.DrawLatex(xcenter, y_id, Form("%d", bin_id));
        }
    }
}

static void DrawBinsForSlice(TH2D* hist,
                             const std::vector<BinDef>& bins,
                             double beta_lo,
                             double beta_hi,
                             bool logX,
                             bool logY) {
    if (!hist || bins.empty()) return;
    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double ymin = hist->GetYaxis()->GetXmin();
    const double ymax = hist->GetYaxis()->GetXmax();
    const double eps = 1e-12;

    TLatex countText;
    countText.SetTextColor(kBlack);
    countText.SetTextSize(0.033);
    countText.SetTextAlign(22);
    countText.SetTextFont(42);

    TLatex idText;
    idText.SetTextColor(kBlue + 2);
    idText.SetTextSize(0.030);
    idText.SetTextAlign(22);
    idText.SetTextFont(42);

    for (const auto& b : bins) {
        if (b.beta_max <= beta_lo + eps || b.beta_min >= beta_hi - eps) {
            continue;
        }
        double xlow = b.xpom_min;
        double xhigh = b.xpom_max;
        double ylow = b.Q2_min;
        double yhigh = b.Q2_max;
        if (xhigh <= xmin || xlow >= xmax || yhigh <= ymin || ylow >= ymax) {
            continue;
        }
        xlow = std::max(xlow, xmin);
        xhigh = std::min(xhigh, xmax);
        ylow = std::max(ylow, ymin);
        yhigh = std::min(yhigh, ymax);
        if (xlow <= 0.0 || ylow <= 0.0) continue;

        TBox* box = new TBox(xlow, ylow, xhigh, yhigh);
        box->SetFillStyle(0);
        box->SetLineColor(kRed);
        box->SetLineWidth(2);
        box->Draw("same");

        const double xcenter = logX ? std::sqrt(xlow * xhigh) : 0.5 * (xlow + xhigh);
        const double ycenter = logY ? std::sqrt(ylow * yhigh) : 0.5 * (ylow + yhigh);
        const int binxlow = hist->GetXaxis()->FindBin(xlow + eps);
        const int binxhigh = hist->GetXaxis()->FindBin(xhigh - eps);
        const int binylow = hist->GetYaxis()->FindBin(ylow + eps);
        const int binyhigh = hist->GetYaxis()->FindBin(yhigh - eps);
        const double count = hist->Integral(binxlow, binxhigh, binylow, binyhigh);
        if (count > 0.0) {
            countText.DrawLatex(xcenter, ycenter, Form("%d", (int)std::lround(count)));
        }
        if (b.bin_id >= 0) {
            const double y_id = logY ? (ycenter * 1.25) : (ycenter + 0.06 * (yhigh - ylow));
            idText.DrawLatex(xcenter, y_id, Form("%d", b.bin_id));
        }
    }
}

static void DrawSliceGrid(TH3D* h3,
                          int sliceAxis,
                          const std::vector<double>& edges,
                          const char* projOpt,
                          const char* xTitle,
                          const char* yTitle,
                          const char* sliceLabel,
                          const char* saveName,
                          bool logX,
                          bool logY,
                          const std::vector<double>* overlayXBins = nullptr,
                          const std::vector<double>* overlayYBins = nullptr,
                          bool overlayLogX = false,
                          bool overlayLogY = false,
                          const std::vector<BinDef>* overlayBins = nullptr) {
    if (!h3 || edges.size() < 2) return;
    const int nSlices = static_cast<int>(edges.size()) - 1;
    const int nCols = 2;
    const int nRows = (nSlices + nCols - 1) / nCols;

    TCanvas* c = new TCanvas(Form("c_phase_slices_%d", sliceAxis), saveName, 2400, 2000);
    c->Divide(nCols, nRows, 0.002, 0.002);

    TAxis* axis = nullptr;
    if (sliceAxis == 1) axis = h3->GetXaxis();
    if (sliceAxis == 2) axis = h3->GetYaxis();
    if (sliceAxis == 3) axis = h3->GetZaxis();
    if (!axis) {
        delete c;
        return;
    }

    std::vector<TH2D*> projections;
    projections.reserve(nSlices);

    for (int i = 0; i < nSlices; ++i) {
        const double lo = edges[i];
        const double hi = edges[i + 1];
        const double tiny = 1e-12;
        const int binLo = axis->FindBin(lo + tiny);
        const int binHi = axis->FindBin(hi - tiny);
        axis->SetRange(binLo, binHi);

        TH2D* h2raw = (TH2D*)h3->Project3D(projOpt);
        if (!h2raw) continue;
        TH2D* h2 = (TH2D*)h2raw->Clone(Form("slice_%s_%d", projOpt, i));
        h2->SetDirectory(nullptr);
        delete h2raw;
        h2->SetTitle("");
        h2->GetXaxis()->SetTitle(xTitle);
        h2->GetYaxis()->SetTitle(yTitle);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.2);
        if (strcmp(xTitle, "#beta") == 0) {
            h2->GetXaxis()->SetRangeUser(0.0, 1.0);
        }
        if (strcmp(yTitle, "#beta") == 0) {
            h2->GetYaxis()->SetRangeUser(0.0, 1.0);
        }

        TPad* pad = (TPad*)c->cd(i + 1);
        pad->SetRightMargin(0.14);
        pad->SetLeftMargin(0.14);
        pad->SetTopMargin(0.12);
        pad->SetBottomMargin(0.12);
        if (logX) pad->SetLogx();
        if (logY) pad->SetLogy();
        h2->Draw("COLZ");

        if (overlayBins) {
            DrawBinsForSlice(h2, *overlayBins, lo, hi, overlayLogX, overlayLogY);
        } else if (overlayXBins && overlayYBins) {
            DrawBinningGridWithCounts(h2, *overlayXBins, *overlayYBins, overlayLogX, overlayLogY);
        }

        TLatex latex;
        latex.SetTextSize(0.06);
        latex.SetNDC();
        latex.SetTextColor(kBlack);
        const std::string label = std::string(sliceLabel) + " " + FormatRange(lo, hi);
        latex.DrawLatex(0.15, 0.88, label.c_str());
        projections.push_back(h2);
    }

    axis->SetRange(0, 0);
    c->Update();
    SaveCanvas(c, saveName);
    for (TH2D* proj : projections) {
        delete proj;
    }
    delete c;
}

static void PlotPhaseSpaceSlices(TFile* inputFile, const std::string& yamlPath) {
    if (!inputFile) return;
    TH3D* h3 = (TH3D*)inputFile->Get("phase3D_reco");
    if (!h3) {
        std::cerr << "Warning: phase3D_reco not found; skipping phase-space slices." << std::endl;
        return;
    }
    if (!yamlPath.empty()) {
        TH3D* h3_yaml = (TH3D*)inputFile->Get("phase3D_reco_yaml");
        if (h3_yaml) {
            h3 = h3_yaml;
        }
    }

    const int nSlices = 4;
    std::vector<double> q2_edges = BuildLogEdges(h3->GetXaxis()->GetXmin(), h3->GetXaxis()->GetXmax(), nSlices);
    std::vector<double> xpom_edges = BuildLogEdges(h3->GetYaxis()->GetXmin(), h3->GetYaxis()->GetXmax(), nSlices);
    std::vector<double> beta_edges = BuildLinEdges(h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetXmax(), nSlices);

    // Editable bin edges for (x_pom, Q^2) overlays in beta slices
    // Modify these vectors to move/add/remove bins; counts update on recompile+plot.
    const int n_xpom_bins_overlay = 4;
    const int n_q2_bins_overlay = 10;
    std::vector<double> xpom_overlay_bins = BuildLogEdges(h3->GetYaxis()->GetXmin(),
                                                         h3->GetYaxis()->GetXmax(),
                                                         n_xpom_bins_overlay);
    std::vector<double> q2_overlay_bins = BuildLogEdges(h3->GetXaxis()->GetXmin(),
                                                       h3->GetXaxis()->GetXmax(),
                                                       n_q2_bins_overlay);
    std::vector<BinDef> yaml_bins;
    const bool yaml_requested = !yamlPath.empty();
    if (!yamlPath.empty()) {
        yaml_bins = ReadBinsFromYAML(yamlPath);
        std::vector<double> y_q2, y_beta, y_xpom;
        CollectEdges(yaml_bins, y_q2, y_beta, y_xpom);
        if (y_q2.size() >= 2 && y_beta.size() >= 2 && y_xpom.size() >= 2) {
            q2_edges = y_q2;
            beta_edges = y_beta;
            xpom_edges = y_xpom;
            xpom_overlay_bins = y_xpom;
            q2_overlay_bins = y_q2;
            std::cout << "Using YAML bins for phase slices: " << yamlPath << std::endl;
        } else {
            std::cerr << "Warning: YAML bins invalid; using default slice edges." << std::endl;
            yaml_bins.clear();
        }
    }

    DrawSliceGrid(h3, 3, beta_edges,
                  "xy",
                  "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                  "figs/diffractive/histos/phase_slices_beta.png",
                  true, true,
                  (yaml_requested ? nullptr : &xpom_overlay_bins),
                  (yaml_requested ? nullptr : &q2_overlay_bins),
                  true, true,
                  yaml_requested ? &yaml_bins : nullptr);
}

static void PlotRecoSetComparison(TFile* inputFile,
                                  const char* effTruthHistNameSetA,
                                  const char* effRecoHistNameSetA,
                                  const char* recoHistNameSetA,
                                  const char* recoHistNameSetB,
                                  const char* xLabel,
                                  const char* title,
                                  const char* saveName,
                                  const bool logX,
                                  const bool logY,
                                  const char* setALabel,
                                  const char* setBLabel) {
    if (!inputFile) return;
    TH1* h_truth_setA = (TH1*)inputFile->Get(effTruthHistNameSetA);
    TH1* h_reco_setA_for_eff = (TH1*)inputFile->Get(effRecoHistNameSetA);
    TH1* h_setA = (TH1*)inputFile->Get(recoHistNameSetA);
    TH1* h_setB = (TH1*)inputFile->Get(recoHistNameSetB);
    if (!h_truth_setA || !h_reco_setA_for_eff || !h_setA || !h_setB) {
        std::cerr << "Warning: Missing hist(s) for common-eff reco comparison: "
                  << effTruthHistNameSetA << ", " << effRecoHistNameSetA << ", "
                  << recoHistNameSetA << ", " << recoHistNameSetB << std::endl;
        return;
    }

    TH1D* h_eff = (TH1D*)h_reco_setA_for_eff->Clone(Form("eff_ref_%s", effRecoHistNameSetA));
    h_eff->SetDirectory(nullptr);
    h_eff->Divide(h_reco_setA_for_eff, h_truth_setA, 1.0, 1.0, "B");

    TH1D* hA_corr = (TH1D*)h_setA->Clone(Form("setA_corr_%s", recoHistNameSetA));
    TH1D* hB_corr = (TH1D*)h_setB->Clone(Form("setB_corr_%s", recoHistNameSetB));
    hA_corr->SetDirectory(nullptr);
    hB_corr->SetDirectory(nullptr);

    const int nbins = hA_corr->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        const double eff = h_eff->GetBinContent(i);
        const double a = hA_corr->GetBinContent(i);
        const double b = hB_corr->GetBinContent(i);
        const double a_err = hA_corr->GetBinError(i);
        const double b_err = hB_corr->GetBinError(i);
        if (eff > 0.0) {
            hA_corr->SetBinContent(i, a / eff);
            hA_corr->SetBinError(i, a_err / eff);
            hB_corr->SetBinContent(i, b / eff);
            hB_corr->SetBinError(i, b_err / eff);
        } else {
            hA_corr->SetBinContent(i, 0.0);
            hA_corr->SetBinError(i, 0.0);
            hB_corr->SetBinContent(i, 0.0);
            hB_corr->SetBinError(i, 0.0);
        }
    }

    TCanvas* c = new TCanvas(Form("c_cmp_%s", recoHistNameSetA), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    hA_corr->SetStats(false);
    hA_corr->SetTitle("");
    hA_corr->GetXaxis()->SetTitle(xLabel);
    hA_corr->GetYaxis()->SetTitle("Efficiency-corrected counts");
    hA_corr->GetXaxis()->SetTitleOffset(1.1);
    hA_corr->GetYaxis()->SetTitleOffset(1.2);

    hA_corr->SetLineWidth(2);
    hA_corr->SetLineColor(kBlack);
    hB_corr->SetMarkerStyle(20);
    hB_corr->SetMarkerSize(1.1);
    hB_corr->SetMarkerColor(kBlack);
    hB_corr->SetLineColor(kBlack);

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= hA_corr->GetNbinsX(); ++i) {
        const double a = hA_corr->GetBinContent(i);
        const double b = hB_corr->GetBinContent(i);
        max_val = std::max(max_val, a + hA_corr->GetBinError(i));
        max_val = std::max(max_val, b + hB_corr->GetBinError(i));
        if (a > 0.0 && a < min_pos) min_pos = a;
        if (b > 0.0 && b < min_pos) min_pos = b;
    }
    if (max_val > 0.0) {
        hA_corr->SetMaximum(max_val * 1.25);
        if (logY) {
            const double min_y = (min_pos < std::numeric_limits<double>::max()) ? (min_pos * 0.5) : 1e-6;
            hA_corr->SetMinimum(std::max(min_y, 1e-6));
        } else {
            hA_corr->SetMinimum(0.0);
        }
    }

    hA_corr->Draw("HIST");
    hB_corr->Draw("PE SAME");

    TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hA_corr, setALabel, "l");
    legend->AddEntry(hB_corr, setBLabel, "p");
    legend->AddEntry((TObject*)nullptr, "Common efficiency from Set A", "");
    legend->Draw();

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
    delete h_eff;
    delete hA_corr;
    delete hB_corr;
}

static void PlotRecoSetComparisonBR(TFile* inputFile,
                                    const char* effTruthHistNameSetA_B0,
                                    const char* effRecoHistNameSetA_B0,
                                    const char* recoHistNameSetA_B0,
                                    const char* recoHistNameSetB_B0,
                                    const char* effTruthHistNameSetA_RP,
                                    const char* effRecoHistNameSetA_RP,
                                    const char* recoHistNameSetA_RP,
                                    const char* recoHistNameSetB_RP,
                                    const char* xLabel,
                                    const char* title,
                                    const char* saveName,
                                    const bool logX,
                                    const bool logY,
                                    const char* setALabel,
                                    const char* setBLabel,
                                    const bool includeSum,
                                    const bool drawSetAAsTruth,
                                    const bool drawSingleSetALine) {
    if (!inputFile) return;
    TH1* h_truthA_b0_src = (TH1*)inputFile->Get(effTruthHistNameSetA_B0);
    TH1* h_recoA_b0_for_eff_src = (TH1*)inputFile->Get(effRecoHistNameSetA_B0);
    TH1* h_truthA_rp_src = (TH1*)inputFile->Get(effTruthHistNameSetA_RP);
    TH1* h_recoA_rp_for_eff_src = (TH1*)inputFile->Get(effRecoHistNameSetA_RP);
    TH1* hA_b0_src = (TH1*)inputFile->Get(recoHistNameSetA_B0);
    TH1* hB_b0_src = (TH1*)inputFile->Get(recoHistNameSetB_B0);
    TH1* hA_rp_src = (TH1*)inputFile->Get(recoHistNameSetA_RP);
    TH1* hB_rp_src = (TH1*)inputFile->Get(recoHistNameSetB_RP);
    if (!h_truthA_b0_src || !h_recoA_b0_for_eff_src || !h_truthA_rp_src || !h_recoA_rp_for_eff_src ||
        !hA_b0_src || !hB_b0_src || !hA_rp_src || !hB_rp_src) {
        std::cerr << "Warning: Missing hist(s) for common-eff B0/RP reco comparison: "
                  << effTruthHistNameSetA_B0 << ", " << effRecoHistNameSetA_B0 << ", "
                  << recoHistNameSetA_B0 << ", " << recoHistNameSetB_B0 << ", "
                  << effTruthHistNameSetA_RP << ", " << effRecoHistNameSetA_RP << ", "
                  << recoHistNameSetA_RP << ", " << recoHistNameSetB_RP << std::endl;
        return;
    }

    TH1D* h_eff_b0 = (TH1D*)h_recoA_b0_for_eff_src->Clone(Form("eff_ref_b0_%s", effRecoHistNameSetA_B0));
    TH1D* h_eff_rp = (TH1D*)h_recoA_rp_for_eff_src->Clone(Form("eff_ref_rp_%s", effRecoHistNameSetA_RP));
    h_eff_b0->SetDirectory(nullptr);
    h_eff_rp->SetDirectory(nullptr);
    h_eff_b0->Divide(h_recoA_b0_for_eff_src, h_truthA_b0_src, 1.0, 1.0, "B");
    h_eff_rp->Divide(h_recoA_rp_for_eff_src, h_truthA_rp_src, 1.0, 1.0, "B");

    TH1D* hA_b0_corr = (TH1D*)hA_b0_src->Clone(Form("setA_b0_corr_%s", recoHistNameSetA_B0));
    TH1D* hB_b0_corr = (TH1D*)hB_b0_src->Clone(Form("setB_b0_corr_%s", recoHistNameSetB_B0));
    TH1D* hA_rp_corr = (TH1D*)hA_rp_src->Clone(Form("setA_rp_corr_%s", recoHistNameSetA_RP));
    TH1D* hB_rp_corr = (TH1D*)hB_rp_src->Clone(Form("setB_rp_corr_%s", recoHistNameSetB_RP));
    hA_b0_corr->SetDirectory(nullptr);
    hB_b0_corr->SetDirectory(nullptr);
    hA_rp_corr->SetDirectory(nullptr);
    hB_rp_corr->SetDirectory(nullptr);

    const int nbins_b0 = hA_b0_corr->GetNbinsX();
    for (int i = 1; i <= nbins_b0; ++i) {
        const double eff = h_eff_b0->GetBinContent(i);
        const double a = hA_b0_corr->GetBinContent(i);
        const double b = hB_b0_corr->GetBinContent(i);
        const double a_err = hA_b0_corr->GetBinError(i);
        const double b_err = hB_b0_corr->GetBinError(i);
        if (eff > 0.0) {
            hA_b0_corr->SetBinContent(i, a / eff);
            hA_b0_corr->SetBinError(i, a_err / eff);
            hB_b0_corr->SetBinContent(i, b / eff);
            hB_b0_corr->SetBinError(i, b_err / eff);
        } else {
            hA_b0_corr->SetBinContent(i, 0.0);
            hA_b0_corr->SetBinError(i, 0.0);
            hB_b0_corr->SetBinContent(i, 0.0);
            hB_b0_corr->SetBinError(i, 0.0);
        }
    }
    const int nbins_rp = hA_rp_corr->GetNbinsX();
    for (int i = 1; i <= nbins_rp; ++i) {
        const double eff = h_eff_rp->GetBinContent(i);
        const double a = hA_rp_corr->GetBinContent(i);
        const double b = hB_rp_corr->GetBinContent(i);
        const double a_err = hA_rp_corr->GetBinError(i);
        const double b_err = hB_rp_corr->GetBinError(i);
        if (eff > 0.0) {
            hA_rp_corr->SetBinContent(i, a / eff);
            hA_rp_corr->SetBinError(i, a_err / eff);
            hB_rp_corr->SetBinContent(i, b / eff);
            hB_rp_corr->SetBinError(i, b_err / eff);
        } else {
            hA_rp_corr->SetBinContent(i, 0.0);
            hA_rp_corr->SetBinError(i, 0.0);
            hB_rp_corr->SetBinContent(i, 0.0);
            hB_rp_corr->SetBinError(i, 0.0);
        }
    }

    TH1D* hA_b0_truth = nullptr;
    TH1D* hA_rp_truth = nullptr;
    TH1D* hA_b0_draw = hA_b0_corr;
    TH1D* hA_rp_draw = hA_rp_corr;
    if (drawSetAAsTruth) {
        hA_b0_truth = (TH1D*)h_truthA_b0_src->Clone(Form("setA_b0_truth_%s", recoHistNameSetA_B0));
        hA_rp_truth = (TH1D*)h_truthA_rp_src->Clone(Form("setA_rp_truth_%s", recoHistNameSetA_RP));
        hA_b0_truth->SetDirectory(nullptr);
        hA_rp_truth->SetDirectory(nullptr);
        hA_b0_draw = hA_b0_truth;
        hA_rp_draw = hA_rp_truth;
    }

    TH1D* h_sum_A = nullptr;
    TH1D* h_sum_B = nullptr;
    if (includeSum) {
        h_sum_A = (TH1D*)hA_b0_draw->Clone(Form("sum_setA_%s", recoHistNameSetA_B0));
        h_sum_A->SetDirectory(nullptr);
        h_sum_A->Add(hA_rp_draw);
        h_sum_B = (TH1D*)hB_b0_corr->Clone(Form("sum_setB_%s", recoHistNameSetB_B0));
        h_sum_B->SetDirectory(nullptr);
        h_sum_B->Add(hB_rp_corr);
    }
    TH1D* hA_single = nullptr;
    if (drawSingleSetALine) {
        hA_single = (TH1D*)hA_b0_draw->Clone(Form("single_setA_%s", recoHistNameSetA_B0));
        hA_single->SetDirectory(nullptr);
        hA_single->Add(hA_rp_draw);
    }
    TH1D* h_ref_draw = drawSingleSetALine ? hA_single : hA_b0_draw;

    TCanvas* c = new TCanvas(Form("c_cmp_%s", recoHistNameSetA_B0), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h_ref_draw->SetStats(false);
    h_ref_draw->SetTitle("");
    h_ref_draw->GetXaxis()->SetTitle(xLabel);
    h_ref_draw->GetYaxis()->SetTitle("Efficiency-corrected counts");
    h_ref_draw->GetXaxis()->SetTitleOffset(1.1);
    h_ref_draw->GetYaxis()->SetTitleOffset(1.2);

    const int color_b0 = kRed + 1;
    const int color_rp = kBlue + 1;
    const int color_setA = drawSingleSetALine ? kBlack : color_b0;

    hA_b0_draw->SetLineWidth(2);
    hA_b0_draw->SetLineColor(color_setA);
    hA_b0_draw->SetLineStyle(1);
    if (drawSingleSetALine) {
        h_ref_draw->SetLineWidth(2);
        h_ref_draw->SetLineColor(kBlack);
        h_ref_draw->SetLineStyle(1);
    }
    hB_b0_corr->SetMarkerStyle(20);
    hB_b0_corr->SetMarkerSize(1.0);
    hB_b0_corr->SetMarkerColor(color_b0);
    hB_b0_corr->SetLineColor(color_b0);
    hA_rp_draw->SetLineWidth(2);
    hA_rp_draw->SetLineColor(color_rp);
    hA_rp_draw->SetLineStyle(2);
    hB_rp_corr->SetMarkerStyle(20);
    hB_rp_corr->SetMarkerSize(1.0);
    hB_rp_corr->SetMarkerColor(color_rp);
    hB_rp_corr->SetLineColor(color_rp);

    if (includeSum) {
        if (!drawSingleSetALine) {
            h_sum_A->SetLineWidth(2);
            h_sum_A->SetLineColor(kGreen + 2);
            h_sum_A->SetLineStyle(3);
        }
        h_sum_B->SetMarkerStyle(21);
        h_sum_B->SetMarkerSize(1.0);
        h_sum_B->SetMarkerColor(kGreen + 2);
        h_sum_B->SetLineColor(kGreen + 2);
    }

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= h_ref_draw->GetNbinsX(); ++i) {
        max_val = std::max(max_val, h_ref_draw->GetBinContent(i) + h_ref_draw->GetBinError(i));
        max_val = std::max(max_val, hB_b0_corr->GetBinContent(i) + hB_b0_corr->GetBinError(i));
        max_val = std::max(max_val, hB_rp_corr->GetBinContent(i) + hB_rp_corr->GetBinError(i));
        if (!drawSingleSetALine) {
            max_val = std::max(max_val, hA_rp_draw->GetBinContent(i) + hA_rp_draw->GetBinError(i));
        }
        if (includeSum) {
            if (!drawSingleSetALine) {
                max_val = std::max(max_val, h_sum_A->GetBinContent(i) + h_sum_A->GetBinError(i));
            }
            max_val = std::max(max_val, h_sum_B->GetBinContent(i) + h_sum_B->GetBinError(i));
        }
        const double vA0 = h_ref_draw->GetBinContent(i);
        const double vB0 = hB_b0_corr->GetBinContent(i);
        const double vBR = hB_rp_corr->GetBinContent(i);
        if (vA0 > 0.0 && vA0 < min_pos) min_pos = vA0;
        if (vB0 > 0.0 && vB0 < min_pos) min_pos = vB0;
        if (vBR > 0.0 && vBR < min_pos) min_pos = vBR;
        if (!drawSingleSetALine) {
            const double vAR = hA_rp_draw->GetBinContent(i);
            if (vAR > 0.0 && vAR < min_pos) min_pos = vAR;
        }
        if (includeSum) {
            if (!drawSingleSetALine) {
                const double vAS = h_sum_A->GetBinContent(i);
                if (vAS > 0.0 && vAS < min_pos) min_pos = vAS;
            }
            const double vBS = h_sum_B->GetBinContent(i);
            if (vBS > 0.0 && vBS < min_pos) min_pos = vBS;
        }
    }
    if (max_val > 0.0) {
        h_ref_draw->SetMaximum(max_val * 1.25);
        if (logY) {
            const double min_y = (min_pos < std::numeric_limits<double>::max()) ? (min_pos * 0.5) : 1e-6;
            h_ref_draw->SetMinimum(std::max(min_y, 1e-6));
        } else {
            h_ref_draw->SetMinimum(0.0);
        }
    }

    h_ref_draw->Draw("HIST");
    if (!drawSingleSetALine) {
        hA_rp_draw->Draw("HIST SAME");
    }
    hB_b0_corr->Draw("PE SAME");
    hB_rp_corr->Draw("PE SAME");
    if (includeSum) {
        if (!drawSingleSetALine) {
            h_sum_A->Draw("HIST SAME");
        }
        h_sum_B->Draw("PE SAME");
    }

    const bool legend_top_left = (strstr(saveName, "xpom_effcorr") != nullptr) ||
                                 (strstr(saveName, "t_effcorr") != nullptr);
    TLegend* legend = legend_top_left
        ? new TLegend(0.15, 0.7, 0.45, 0.9)
        : new TLegend(0.55, 0.65, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    if (drawSingleSetALine) {
        legend->AddEntry(h_ref_draw, setALabel, "l");
    } else {
        legend->AddEntry(hA_b0_draw, Form("B0 %s", setALabel), "l");
        legend->AddEntry(hA_rp_draw, Form("RP %s", setALabel), "l");
    }
    legend->AddEntry(hB_b0_corr, Form("B0 %s", setBLabel), "p");
    legend->AddEntry(hB_rp_corr, Form("RP %s", setBLabel), "p");
    if (includeSum) {
        if (!drawSingleSetALine) {
            legend->AddEntry(h_sum_A, Form("B0+RP %s", setALabel), "l");
        }
        legend->AddEntry(h_sum_B, Form("B0+RP %s", setBLabel), "p");
    }
    legend->AddEntry((TObject*)nullptr, "Common efficiency from Set A", "");
    legend->Draw();

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
    delete h_eff_b0;
    delete h_eff_rp;
    delete hA_b0_corr;
    delete hB_b0_corr;
    delete hA_rp_corr;
    delete hB_rp_corr;
    delete hA_b0_truth;
    delete hA_rp_truth;
    delete h_sum_A;
    delete h_sum_B;
    delete hA_single;
}

static void PlotRecoSetComparisonBRWithSetBUncorrected(TFile* inputFile,
                                                       const char* effTruthHistNameSetA_B0,
                                                       const char* effRecoHistNameSetA_B0,
                                                       const char* recoHistNameSetA_B0,
                                                       const char* recoHistNameSetB_B0,
                                                       const char* effTruthHistNameSetA_RP,
                                                       const char* effRecoHistNameSetA_RP,
                                                       const char* recoHistNameSetA_RP,
                                                       const char* recoHistNameSetB_RP,
                                                       const char* xLabel,
                                                       const char* title,
                                                       const char* saveName,
                                                       const bool logX,
                                                       const bool logY) {
    if (!inputFile) return;

    TH1* h_truthA_b0_src = (TH1*)inputFile->Get(effTruthHistNameSetA_B0);
    TH1* h_recoA_b0_for_eff_src = (TH1*)inputFile->Get(effRecoHistNameSetA_B0);
    TH1* h_truthA_rp_src = (TH1*)inputFile->Get(effTruthHistNameSetA_RP);
    TH1* h_recoA_rp_for_eff_src = (TH1*)inputFile->Get(effRecoHistNameSetA_RP);
    TH1* hA_b0_src = (TH1*)inputFile->Get(recoHistNameSetA_B0);
    TH1* hB_b0_src = (TH1*)inputFile->Get(recoHistNameSetB_B0);
    TH1* hA_rp_src = (TH1*)inputFile->Get(recoHistNameSetA_RP);
    TH1* hB_rp_src = (TH1*)inputFile->Get(recoHistNameSetB_RP);
    if (!h_truthA_b0_src || !h_recoA_b0_for_eff_src || !h_truthA_rp_src || !h_recoA_rp_for_eff_src ||
        !hA_b0_src || !hB_b0_src || !hA_rp_src || !hB_rp_src) {
        std::cerr << "Warning: Missing hist(s) for B0/RP set-B corrected/uncorrected comparison: "
                  << effTruthHistNameSetA_B0 << ", " << effRecoHistNameSetA_B0 << ", "
                  << recoHistNameSetA_B0 << ", " << recoHistNameSetB_B0 << ", "
                  << effTruthHistNameSetA_RP << ", " << effRecoHistNameSetA_RP << ", "
                  << recoHistNameSetA_RP << ", " << recoHistNameSetB_RP << std::endl;
        return;
    }

    TH1D* h_eff_b0 = (TH1D*)h_recoA_b0_for_eff_src->Clone(Form("eff_ref_b0_var_%s", effRecoHistNameSetA_B0));
    TH1D* h_eff_rp = (TH1D*)h_recoA_rp_for_eff_src->Clone(Form("eff_ref_rp_var_%s", effRecoHistNameSetA_RP));
    h_eff_b0->SetDirectory(nullptr);
    h_eff_rp->SetDirectory(nullptr);
    h_eff_b0->Divide(h_recoA_b0_for_eff_src, h_truthA_b0_src, 1.0, 1.0, "B");
    h_eff_rp->Divide(h_recoA_rp_for_eff_src, h_truthA_rp_src, 1.0, 1.0, "B");

    TH1D* hA_sum_reco = (TH1D*)hA_b0_src->Clone(Form("sum_setA_reco_%s", recoHistNameSetA_B0));
    hA_sum_reco->SetDirectory(nullptr);
    hA_sum_reco->Add(hA_rp_src);

    TH1D* hB_b0_corr = (TH1D*)hB_b0_src->Clone(Form("setB_b0_corr_var_%s", recoHistNameSetB_B0));
    TH1D* hB_rp_corr = (TH1D*)hB_rp_src->Clone(Form("setB_rp_corr_var_%s", recoHistNameSetB_RP));
    TH1D* hB_b0_unc = (TH1D*)hB_b0_src->Clone(Form("setB_b0_unc_var_%s", recoHistNameSetB_B0));
    TH1D* hB_rp_unc = (TH1D*)hB_rp_src->Clone(Form("setB_rp_unc_var_%s", recoHistNameSetB_RP));
    hB_b0_corr->SetDirectory(nullptr);
    hB_rp_corr->SetDirectory(nullptr);
    hB_b0_unc->SetDirectory(nullptr);
    hB_rp_unc->SetDirectory(nullptr);

    for (int i = 1; i <= hB_b0_corr->GetNbinsX(); ++i) {
        const double eff = h_eff_b0->GetBinContent(i);
        const double v = hB_b0_corr->GetBinContent(i);
        const double e = hB_b0_corr->GetBinError(i);
        if (eff > 0.0) {
            hB_b0_corr->SetBinContent(i, v / eff);
            hB_b0_corr->SetBinError(i, e / eff);
        } else {
            hB_b0_corr->SetBinContent(i, 0.0);
            hB_b0_corr->SetBinError(i, 0.0);
        }
    }
    for (int i = 1; i <= hB_rp_corr->GetNbinsX(); ++i) {
        const double eff = h_eff_rp->GetBinContent(i);
        const double v = hB_rp_corr->GetBinContent(i);
        const double e = hB_rp_corr->GetBinError(i);
        if (eff > 0.0) {
            hB_rp_corr->SetBinContent(i, v / eff);
            hB_rp_corr->SetBinError(i, e / eff);
        } else {
            hB_rp_corr->SetBinContent(i, 0.0);
            hB_rp_corr->SetBinError(i, 0.0);
        }
    }

    TCanvas* c = new TCanvas(Form("c_cmp_unc_%s", recoHistNameSetA_B0), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    hA_sum_reco->SetStats(false);
    hA_sum_reco->SetTitle("");
    hA_sum_reco->GetXaxis()->SetTitle(xLabel);
    hA_sum_reco->GetYaxis()->SetTitle("Counts");
    hA_sum_reco->GetXaxis()->SetTitleOffset(1.1);
    hA_sum_reco->GetYaxis()->SetTitleOffset(1.2);
    hA_sum_reco->SetLineColor(kBlack);
    hA_sum_reco->SetLineWidth(2);
    hA_sum_reco->SetLineStyle(1);

    const int color_b0 = kRed + 1;
    const int color_rp = kBlue + 1;
    hB_b0_corr->SetMarkerStyle(20);
    hB_b0_corr->SetMarkerSize(1.0);
    hB_b0_corr->SetMarkerColor(color_b0);
    hB_b0_corr->SetLineColor(color_b0);
    hB_b0_unc->SetMarkerStyle(24);
    hB_b0_unc->SetMarkerSize(1.0);
    hB_b0_unc->SetMarkerColor(color_b0);
    hB_b0_unc->SetLineColor(color_b0);

    hB_rp_corr->SetMarkerStyle(20);
    hB_rp_corr->SetMarkerSize(1.0);
    hB_rp_corr->SetMarkerColor(color_rp);
    hB_rp_corr->SetLineColor(color_rp);
    hB_rp_unc->SetMarkerStyle(24);
    hB_rp_unc->SetMarkerSize(1.0);
    hB_rp_unc->SetMarkerColor(color_rp);
    hB_rp_unc->SetLineColor(color_rp);

    auto updateRange = [](TH1* h, double& maxVal, double& minPos) {
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            const double c = h->GetBinContent(i);
            const double e = h->GetBinError(i);
            maxVal = std::max(maxVal, c + e);
            if (c > 0.0 && c < minPos) minPos = c;
        }
    };
    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    updateRange(hA_sum_reco, max_val, min_pos);
    updateRange(hB_b0_corr, max_val, min_pos);
    updateRange(hB_b0_unc, max_val, min_pos);
    updateRange(hB_rp_corr, max_val, min_pos);
    updateRange(hB_rp_unc, max_val, min_pos);
    if (max_val > 0.0) {
        hA_sum_reco->SetMaximum(max_val * 1.25);
        if (logY) {
            const double min_y = (min_pos < std::numeric_limits<double>::max()) ? (min_pos * 0.5) : 1e-6;
            hA_sum_reco->SetMinimum(std::max(min_y, 1e-6));
        } else {
            hA_sum_reco->SetMinimum(0.0);
        }
    }

    hA_sum_reco->Draw("HIST");
    hB_b0_unc->Draw("PE SAME");
    hB_b0_corr->Draw("PE SAME");
    hB_rp_unc->Draw("PE SAME");
    hB_rp_corr->Draw("PE SAME");

    TLegend* legend = new TLegend(0.15, 0.65, 0.48, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hA_sum_reco, "MC (Set-A reco sum)", "l");
    legend->AddEntry(hB_b0_corr, "B0 pseudo-data corrected", "p");
    legend->AddEntry(hB_b0_unc, "B0 pseudo-data uncorrected", "p");
    legend->AddEntry(hB_rp_corr, "RP pseudo-data corrected", "p");
    legend->AddEntry(hB_rp_unc, "RP pseudo-data uncorrected", "p");
    legend->AddEntry((TObject*)nullptr, "Efficiency from Set A", "");
    legend->Draw();

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
    delete h_eff_b0;
    delete h_eff_rp;
    delete hA_sum_reco;
    delete hB_b0_corr;
    delete hB_rp_corr;
    delete hB_b0_unc;
    delete hB_rp_unc;
}

static void PlotRecoSetComparisonBRSumOnly(TFile* inputFile,
                                           const char* effTruthHistNameSetA_B0,
                                           const char* effRecoHistNameSetA_B0,
                                           const char* recoHistNameSetA_B0,
                                           const char* recoHistNameSetB_B0,
                                           const char* effTruthHistNameSetA_RP,
                                           const char* effRecoHistNameSetA_RP,
                                           const char* recoHistNameSetA_RP,
                                           const char* recoHistNameSetB_RP,
                                           const char* xLabel,
                                           const char* title,
                                           const char* saveName,
                                           const bool logX,
                                           const bool logY,
                                           const char* setALabel,
                                           const char* setBLabel) {
    if (!inputFile) return;
    TH1* h_truthA_b0_src = (TH1*)inputFile->Get(effTruthHistNameSetA_B0);
    TH1* h_recoA_b0_for_eff_src = (TH1*)inputFile->Get(effRecoHistNameSetA_B0);
    TH1* h_truthA_rp_src = (TH1*)inputFile->Get(effTruthHistNameSetA_RP);
    TH1* h_recoA_rp_for_eff_src = (TH1*)inputFile->Get(effRecoHistNameSetA_RP);
    TH1* hA_b0_src = (TH1*)inputFile->Get(recoHistNameSetA_B0);
    TH1* hB_b0_src = (TH1*)inputFile->Get(recoHistNameSetB_B0);
    TH1* hA_rp_src = (TH1*)inputFile->Get(recoHistNameSetA_RP);
    TH1* hB_rp_src = (TH1*)inputFile->Get(recoHistNameSetB_RP);
    if (!h_truthA_b0_src || !h_recoA_b0_for_eff_src || !h_truthA_rp_src || !h_recoA_rp_for_eff_src ||
        !hA_b0_src || !hB_b0_src || !hA_rp_src || !hB_rp_src) {
        std::cerr << "Warning: Missing hist(s) for sum-only B0/RP reco comparison." << std::endl;
        return;
    }

    TH1D* h_eff_b0 = (TH1D*)h_recoA_b0_for_eff_src->Clone(Form("eff_ref_b0_sumonly_%s", effRecoHistNameSetA_B0));
    TH1D* h_eff_rp = (TH1D*)h_recoA_rp_for_eff_src->Clone(Form("eff_ref_rp_sumonly_%s", effRecoHistNameSetA_RP));
    h_eff_b0->SetDirectory(nullptr);
    h_eff_rp->SetDirectory(nullptr);
    h_eff_b0->Divide(h_recoA_b0_for_eff_src, h_truthA_b0_src, 1.0, 1.0, "B");
    h_eff_rp->Divide(h_recoA_rp_for_eff_src, h_truthA_rp_src, 1.0, 1.0, "B");

    TH1D* hA_b0_corr = (TH1D*)hA_b0_src->Clone(Form("setA_b0_corr_sumonly_%s", recoHistNameSetA_B0));
    TH1D* hB_b0_corr = (TH1D*)hB_b0_src->Clone(Form("setB_b0_corr_sumonly_%s", recoHistNameSetB_B0));
    TH1D* hA_rp_corr = (TH1D*)hA_rp_src->Clone(Form("setA_rp_corr_sumonly_%s", recoHistNameSetA_RP));
    TH1D* hB_rp_corr = (TH1D*)hB_rp_src->Clone(Form("setB_rp_corr_sumonly_%s", recoHistNameSetB_RP));
    hA_b0_corr->SetDirectory(nullptr);
    hB_b0_corr->SetDirectory(nullptr);
    hA_rp_corr->SetDirectory(nullptr);
    hB_rp_corr->SetDirectory(nullptr);

    for (int i = 1; i <= hA_b0_corr->GetNbinsX(); ++i) {
        const double eff = h_eff_b0->GetBinContent(i);
        if (eff > 0.0) {
            hA_b0_corr->SetBinContent(i, hA_b0_corr->GetBinContent(i) / eff);
            hA_b0_corr->SetBinError(i, hA_b0_corr->GetBinError(i) / eff);
            hB_b0_corr->SetBinContent(i, hB_b0_corr->GetBinContent(i) / eff);
            hB_b0_corr->SetBinError(i, hB_b0_corr->GetBinError(i) / eff);
        } else {
            hA_b0_corr->SetBinContent(i, 0.0); hA_b0_corr->SetBinError(i, 0.0);
            hB_b0_corr->SetBinContent(i, 0.0); hB_b0_corr->SetBinError(i, 0.0);
        }
    }
    for (int i = 1; i <= hA_rp_corr->GetNbinsX(); ++i) {
        const double eff = h_eff_rp->GetBinContent(i);
        if (eff > 0.0) {
            hA_rp_corr->SetBinContent(i, hA_rp_corr->GetBinContent(i) / eff);
            hA_rp_corr->SetBinError(i, hA_rp_corr->GetBinError(i) / eff);
            hB_rp_corr->SetBinContent(i, hB_rp_corr->GetBinContent(i) / eff);
            hB_rp_corr->SetBinError(i, hB_rp_corr->GetBinError(i) / eff);
        } else {
            hA_rp_corr->SetBinContent(i, 0.0); hA_rp_corr->SetBinError(i, 0.0);
            hB_rp_corr->SetBinContent(i, 0.0); hB_rp_corr->SetBinError(i, 0.0);
        }
    }

    TH1D* h_sum_B_uncorr = (TH1D*)hB_b0_src->Clone(Form("sumB_uncorr_%s", recoHistNameSetB_B0));
    TH1D* h_sum_A_corr = (TH1D*)hA_b0_corr->Clone(Form("sumA_corr_%s", recoHistNameSetA_B0));
    TH1D* h_sum_B_corr = (TH1D*)hB_b0_corr->Clone(Form("sumB_corr_%s", recoHistNameSetB_B0));
    h_sum_B_uncorr->SetDirectory(nullptr);
    h_sum_A_corr->SetDirectory(nullptr);
    h_sum_B_corr->SetDirectory(nullptr);
    h_sum_B_uncorr->Add(hB_rp_src);
    h_sum_A_corr->Add(hA_rp_corr);
    h_sum_B_corr->Add(hB_rp_corr);

    TCanvas* c = new TCanvas(Form("c_cmp_sumonly_%s", recoHistNameSetA_B0), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h_sum_A_corr->SetStats(false);
    h_sum_A_corr->SetTitle("");
    h_sum_A_corr->GetXaxis()->SetTitle(xLabel);
    h_sum_A_corr->GetYaxis()->SetTitle("Counts");
    h_sum_A_corr->GetXaxis()->SetTitleOffset(1.1);
    h_sum_A_corr->GetYaxis()->SetTitleOffset(1.2);

    h_sum_A_corr->SetLineColor(kBlack);
    h_sum_A_corr->SetLineStyle(1);
    h_sum_A_corr->SetLineWidth(2);

    h_sum_B_uncorr->SetMarkerStyle(24);
    h_sum_B_uncorr->SetMarkerSize(1.0);
    h_sum_B_uncorr->SetMarkerColor(kBlue + 1);
    h_sum_B_uncorr->SetLineColor(kBlue + 1);
    h_sum_B_corr->SetMarkerStyle(20);
    h_sum_B_corr->SetMarkerSize(1.0);
    h_sum_B_corr->SetMarkerColor(kRed + 1);
    h_sum_B_corr->SetLineColor(kRed + 1);

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= h_sum_A_corr->GetNbinsX(); ++i) {
        const double vals[] = {
            h_sum_A_corr->GetBinContent(i),
            h_sum_B_uncorr->GetBinContent(i),
            h_sum_B_corr->GetBinContent(i)
        };
        const double errs[] = {
            h_sum_A_corr->GetBinError(i),
            h_sum_B_uncorr->GetBinError(i),
            h_sum_B_corr->GetBinError(i)
        };
        for (int j = 0; j < 3; ++j) {
            max_val = std::max(max_val, vals[j] + errs[j]);
            if (vals[j] > 0.0 && vals[j] < min_pos) min_pos = vals[j];
        }
    }
    if (max_val > 0.0) {
        h_sum_A_corr->SetMaximum(max_val * 1.25);
        if (logY) {
            const double min_y = (min_pos < std::numeric_limits<double>::max()) ? (min_pos * 0.5) : 1e-6;
            h_sum_A_corr->SetMinimum(std::max(min_y, 1e-6));
        } else {
            h_sum_A_corr->SetMinimum(0.0);
        }
    }

    h_sum_A_corr->Draw("HIST");
    h_sum_B_corr->Draw("PE SAME");
    h_sum_B_uncorr->Draw("PE SAME");

    TLegend* legend = new TLegend(0.52, 0.64, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h_sum_A_corr, setALabel, "l");
    legend->AddEntry(h_sum_B_corr, Form("%s corrected (B0+RP)", setBLabel), "p");
    legend->AddEntry(h_sum_B_uncorr, Form("%s uncorrected (B0+RP)", setBLabel), "p");
    legend->AddEntry((TObject*)nullptr, "Common efficiency from Set A", "");
    legend->Draw();

    DrawSimLabels(inputFile);
    c->Update();
    SaveCanvas(c, saveName);

    delete c;
    delete h_eff_b0;
    delete h_eff_rp;
    delete hA_b0_corr;
    delete hB_b0_corr;
    delete hA_rp_corr;
    delete hB_rp_corr;
    delete h_sum_B_uncorr;
    delete h_sum_A_corr;
    delete h_sum_B_corr;
}

static void PlotRelResVsK(TFile* inputFile,
                          const std::vector<std::string>& histNames,
                          const std::vector<std::string>& labels,
                          const std::string& title,
                          const std::string& outpath) {
    if (!inputFile) return;
    std::vector<TH1*> hists;
    hists.reserve(histNames.size());
    for (const auto& name : histNames) {
        TH1* h = dynamic_cast<TH1*>(inputFile->Get(name.c_str()));
        if (h) hists.push_back(h);
    }
    if (hists.empty()) return;

    TCanvas c("c_relres_vs_k", "c_relres_vs_k", 1100, 600);
    gStyle->SetOptTitle(0);
    c.SetGridx();
    c.SetGridy();

    double ymax = 0.0;
    for (const auto* h : hists) {
        if (!h) continue;
        ymax = std::max(ymax, h->GetMaximum());
    }
    if (ymax <= 0.0) ymax = 0.1;

    const int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
    for (size_t i = 0; i < hists.size(); i++) {
        TH1* h = hists[i];
        if (!h) continue;
        h->SetTitle("");
        h->SetMarkerStyle(20 + static_cast<int>(i));
        h->SetMarkerSize(0.9);
        h->SetLineWidth(2);
        h->SetMarkerColor(colors[i % 4]);
        h->SetLineColor(colors[i % 4]);
        h->GetYaxis()->SetRangeUser(0.0, 1.2 * ymax);
        if (i == 0) {
            h->Draw("E1");
        } else {
            h->Draw("E1 SAME");
        }
    }

    if (labels.size() == hists.size()) {
        TLegend leg(0.7, 0.78, 0.92, 0.92);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        for (size_t i = 0; i < hists.size(); i++) {
            leg.AddEntry(hists[i], labels[i].c_str(), "lep");
        }
        leg.Draw();
    }

    DrawSimLabels(inputFile);
    SaveCanvas(&c, outpath.c_str());
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <combined.root> <bins.yaml>" << std::endl;
        return 1;
    }

    gErrorIgnoreLevel = kWarning;

    TString inputFileName = argv[1];

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetNumberContours(120);
    gStyle->SetPalette(kBlueRedYellow);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.048, "XYZ");
    gStyle->SetLabelSize(0.038, "XYZ");
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetOptFit(111);
    SetGreenYellowRedPalette();

    TFile* inputFile = TFile::Open(inputFileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return 1;
    }

    std::string yamlPath = argv[2];
    std::vector<BinDef> yaml_bins = ReadBinsFromYAML(yamlPath);
    std::vector<double> q2_edges;
    std::vector<double> beta_edges;
    std::vector<double> xpom_edges;
    if (!yaml_bins.empty()) {
        CollectEdges(yaml_bins, q2_edges, beta_edges, xpom_edges);
    } else {
        std::cerr << "Error: No 3D_bins loaded from " << yamlPath << std::endl;
        return 1;
    }

    std::vector<PlotOptions*> plots;
    PlotOptions1D* plot_ptr = nullptr;
    PlotOptionsBinnedRelRes* binned_plot_ptr = nullptr;

    // Plot outputs are routed to figs/organized/... in SaveCanvas.
    gSystem->mkdir("figs/organized", kTRUE);

    // Acceptance/purity plots (if tracking histograms exist)
    TH1D* h_gen_Q2 = (TH1D*)inputFile->Get("h_gen_Q2");
    TH1D* h_gen_and_reco_after_cuts_Q2_EM = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_EM");
    TH1D* h_gen_and_reco_after_cuts_Q2_DA = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_DA");
    TH1D* h_gen_and_reco_after_cuts_Q2_Sigma = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_Sigma");
    TH1D* h_reco_Q2_EM = (TH1D*)inputFile->Get("h_reco_Q2_EM");
    TH1D* h_reco_Q2_DA = (TH1D*)inputFile->Get("h_reco_Q2_DA");
    TH1D* h_reco_Q2_Sigma = (TH1D*)inputFile->Get("h_reco_Q2_Sigma");
    const char* controlSetALabel = "MC";
    const char* controlSetBLabel = "pseudo-data";

    PlotXQ2Density(inputFile);
    PlotRecoSetComparisonBRSumOnly(inputFile,
                                   "t_truth_mc_B0",
                                   "t_reco_mc_B0",
                                   "t_reco_mc_B0",
                                   "t_reco_pdata_B0",
                                   "t_truth_mc_RP",
                                   "t_reco_mc_RP",
                                   "t_reco_mc_RP",
                                   "t_reco_pdata_RP",
                                   "|t| [GeV^{2}]",
                                   "|t| (B0+RP) Uncorrected/Corrected",
                                   "figs/diffractive/efficiency/t_effcorr.png",
                                   true,
                                   false,
                                   controlSetALabel,
                                   controlSetBLabel);
    PlotRecoSetComparisonBRSumOnly(inputFile,
                                   "beta_truth_mc_B0",
                                   "beta_reco_mc_B0",
                                   "beta_reco_mc_B0",
                                   "beta_reco_pdata_B0",
                                   "beta_truth_mc_RP",
                                   "beta_reco_mc_RP",
                                   "beta_reco_mc_RP",
                                   "beta_reco_pdata_RP",
                                   "#beta",
                                   "#beta (B0+RP) Uncorrected/Corrected",
                                   "figs/diffractive/efficiency/beta_effcorr.png",
                                   false,
                                   false,
                                   controlSetALabel,
                                   controlSetBLabel);
    PlotRecoSetComparisonBR(inputFile,
                            "xpom_truth_mc_B0",
                            "xpom_reco_mc_B0",
                            "xpom_reco_mc_B0",
                            "xpom_reco_pdata_B0",
                            "xpom_truth_mc_RP",
                            "xpom_reco_mc_RP",
                            "xpom_reco_mc_RP",
                            "xpom_reco_pdata_RP",
                            "x_{pom}",
                            "x_{pom} Reco Comparison (MC vs pseudo-data)",
                            "figs/diffractive/efficiency/xpom_effcorr.png",
                            true,
                            false,
                            controlSetALabel,
                            controlSetBLabel,
                            true,
                            false,
                            true);
    PlotRecoSetComparisonBR(inputFile,
                            "theta_truth_mc_B0",
                            "theta_reco_mc_B0",
                            "theta_reco_mc_B0",
                            "theta_reco_pdata_B0",
                            "theta_truth_mc_RP",
                            "theta_reco_mc_RP",
                            "theta_reco_mc_RP",
                            "theta_reco_pdata_RP",
                            "#theta [mrad]",
                            "#theta Reco Comparison (MC vs pseudo-data)",
                            "figs/diffractive/efficiency/theta_effcorr.png",
                            false,
                            true,
                            controlSetALabel,
                            controlSetBLabel,
                            false,
                            true,
                            true);

    PlotRecoSetComparison(inputFile,
                          "Q2_truth_mc",
                          "Q2_reco_mc",
                          "Q2_reco_mc",
                          "Q2_reco_pdata",
                          "Q^{2} [GeV^{2}]",
                          "Q^{2} Reco Comparison (MC vs pseudo-data)",
                          "figs/inclusive/efficiency/q2_effcorr.png",
                          true,
                          false,
                          controlSetALabel,
                          controlSetBLabel);
    PlotRecoSetComparison(inputFile,
                          "x_truth_mc",
                          "x_reco_mc",
                          "x_reco_mc",
                          "x_reco_pdata",
                          "x_{Bj}",
                          "x_{Bj} Reco Comparison (MC vs pseudo-data)",
                          "figs/inclusive/efficiency/xbj_effcorr.png",
                          true,
                          false,
                          controlSetALabel,
                          controlSetBLabel);
    PlotRecoSetComparison(inputFile,
                          "y_truth_mc",
                          "y_reco_mc",
                          "y_reco_mc",
                          "y_reco_pdata",
                          "y",
                          "y Reco Comparison (MC vs pseudo-data)",
                          "figs/inclusive/efficiency/y_effcorr.png",
                          false,
                          false,
                          controlSetALabel,
                          controlSetBLabel);
    PlotRecoSetComparison(inputFile,
                          "W2_truth_mc",
                          "W2_reco_mc",
                          "W2_reco_mc",
                          "W2_reco_pdata",
                          "W^{2} [GeV^{2}]",
                          "W^{2} Reco Comparison (MC vs pseudo-data)",
                          "figs/inclusive/efficiency/w2_effcorr.png",
                          true,
                          false,
                          controlSetALabel,
                          controlSetBLabel);
    PlotRecoSetComparison(inputFile,
                          "MX2_truth_mc",
                          "MX2_reco_mc",
                          "MX2_reco_mc",
                          "MX2_reco_pdata",
                          "M_{X}^{2} [GeV^{2}]",
                          "M_{X}^{2} Reco Comparison (MC vs pseudo-data)",
                          "figs/diffractive/efficiency/mx2_effcorr.png",
                          true,
                          false,
                          controlSetALabel,
                          controlSetBLabel);
    PlotDensityFromHistWithOverlay(inputFile, "beta_Q2_reco", "#beta", "Q^{2} [GeV^{2}]",
                                   "figs/diffractive/histos/beta_q2_density.png",
                                   false, true, beta_edges, q2_edges, false, true);
    PlotDensityFromHist(inputFile, "t_Q2_reco", "|t| [GeV^{2}]", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/t_q2_density.png", true, true);
    PlotDensityFromHistWithOverlay(inputFile, "xpom_Q2_reco", "x_{pom}", "Q^{2} [GeV^{2}]",
                                   "figs/diffractive/histos/xpom_q2_density.png",
                                   true, true, xpom_edges, q2_edges, true, true);
    PlotDensityFromHist(inputFile, "beta_t_reco", "#beta", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/beta_t_density.png", false, true);
    PlotDensityFromHist(inputFile, "xbj_t_reco", "x_{Bj}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xbj_t_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_t_reco", "x_{pom}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xpom_t_density.png", true, true);
    PlotDensityFromHistWithOverlay(inputFile, "xpom_beta_reco", "x_{pom}", "#beta",
                                   "figs/diffractive/histos/xpom_beta_density.png",
                                   true, false, xpom_edges, beta_edges, true, false);
    PlotDensityFromHist(inputFile, "xbj_beta_reco", "x_{Bj}", "#beta",
                        "figs/diffractive/histos/xbj_beta_density.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_xpom_reco", "x_{Bj}", "x_{pom}",
                        "figs/diffractive/histos/xbj_xpom_density.png", true, true);
    // MC truth density versions
    PlotDensityFromHist(inputFile, "beta_Q2_truth", "#beta", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/beta_q2_density_truth.png", false, true);
    PlotDensityFromHist(inputFile, "t_Q2_truth", "|t| [GeV^{2}]", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/t_q2_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_Q2_truth", "x_{pom}", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/xpom_q2_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "beta_t_truth", "#beta", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/beta_t_density_truth.png", false, true);
    PlotDensityFromHist(inputFile, "xbj_t_truth", "x_{Bj}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xbj_t_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_t_truth", "x_{pom}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xpom_t_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_beta_truth", "x_{pom}", "#beta",
                        "figs/diffractive/histos/xpom_beta_density_truth.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_beta_truth", "x_{Bj}", "#beta",
                        "figs/diffractive/histos/xbj_beta_density_truth.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_xpom_truth", "x_{Bj}", "x_{pom}",
                        "figs/diffractive/histos/xbj_xpom_density_truth.png", true, true);
    PlotPhaseSpaceSlices(inputFile, yamlPath);

    PlotGraphDensity(inputFile,
                     "g_W2_EM",
                     "W^{2} Correlation (EM, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/histos/w2_corr_unbinned_em_density.png",
                     10.0, 1.0e4, 140, true, true);
    PlotGraphDensity(inputFile,
                     "g_W2_DA",
                     "W^{2} Correlation (DA, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/histos/w2_corr_unbinned_da_density.png",
                     10.0, 1.0e4, 140, true, true);
    PlotGraphDensity(inputFile,
                     "g_W2_Sigma",
                     "W^{2} Correlation (#Sigma, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/histos/w2_corr_unbinned_sigma_density.png",
                     10.0, 1.0e4, 140, true, true);

    if (h_gen_Q2 && h_gen_and_reco_after_cuts_Q2_EM && h_reco_Q2_EM) {
        TH1D* h_acceptance_Q2_EM = (TH1D*)h_gen_and_reco_after_cuts_Q2_EM->Clone("h_acceptance_Q2_EM");
        h_acceptance_Q2_EM->Divide(h_gen_Q2);
        h_acceptance_Q2_EM->SetTitle("Acceptance vs Q^{2} (EM);Q^{2} [GeV^{2}];Acceptance");

        TH1D* h_acceptance_Q2_DA = (TH1D*)h_gen_and_reco_after_cuts_Q2_DA->Clone("h_acceptance_Q2_DA");
        h_acceptance_Q2_DA->Divide(h_gen_Q2);

        TH1D* h_acceptance_Q2_Sigma = (TH1D*)h_gen_and_reco_after_cuts_Q2_Sigma->Clone("h_acceptance_Q2_Sigma");
        h_acceptance_Q2_Sigma->Divide(h_gen_Q2);

        TH1D* h_purity_Q2_EM = (TH1D*)h_gen_and_reco_after_cuts_Q2_EM->Clone("h_purity_Q2_EM");
        h_purity_Q2_EM->Divide(h_reco_Q2_EM);
        h_purity_Q2_EM->SetTitle("Purity vs Q^{2} (EM);Q^{2} [GeV^{2}];Purity");

        TH1D* h_purity_Q2_DA = (TH1D*)h_gen_and_reco_after_cuts_Q2_DA->Clone("h_purity_Q2_DA");
        h_purity_Q2_DA->Divide(h_reco_Q2_DA);

        TH1D* h_purity_Q2_Sigma = (TH1D*)h_gen_and_reco_after_cuts_Q2_Sigma->Clone("h_purity_Q2_Sigma");
        h_purity_Q2_Sigma->Divide(h_reco_Q2_Sigma);

        plot_ptr = new PlotOptions1D(
            {"h_acceptance_Q2_EM", "h_acceptance_Q2_DA", "h_acceptance_Q2_Sigma"},
            {"EM Method", "DA Method", "Sigma Method"},
            {"pe", "pe", "pe"},
            "Acceptance vs Q^{2}",
            "Q^{2} [GeV^{2}]",
            "Acceptance",
            "figs/performance/acceptance_vs_Q2.png",
            true,
            false
        );
        plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
        plots.push_back(plot_ptr);

        plot_ptr = new PlotOptions1D(
            {"h_purity_Q2_EM", "h_purity_Q2_DA", "h_purity_Q2_Sigma"},
            {"EM Method", "DA Method", "Sigma Method"},
            {"pe", "pe", "pe"},
            "Purity vs Q^{2}",
            "Q^{2} [GeV^{2}]",
            "Purity",
            "figs/performance/purity_vs_Q2.png",
            true,
            false
        );
        plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
        plots.push_back(plot_ptr);
    }

    // =================================================================
    // Q2/xy distributions and PDFs
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA", "h_Q2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "Q^{2} Reconstruction Methods",
        "Q^{2}",
        "# of events",
        "figs/inclusive/histos/q2_methods_hist.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA", "h_Q2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "Q^{2} PDF Comparison",
        "Q^{2}",
        "PDF",
        "figs/inclusive/histos/q2_methods_pdf.png",
        true,
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"x_truth", "x_EM", "x_DA", "x_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{Bj} Reconstruction Methods",
        "x_{Bj}",
        "# of events",
        "figs/inclusive/histos/xbj_methods_hist.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"x_truth", "x_EM", "x_DA", "x_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{Bj} PDF Comparison",
        "x_{Bj}",
        "PDF",
        "figs/inclusive/histos/xbj_methods_pdf.png",
        true,
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"y_truth", "y_EM", "y_DA", "y_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "y (inelasticity) Reconstruction Methods",
        "y",
        "# of events",
        "figs/inclusive/histos/y_methods_hist.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"y_truth", "y_EM", "y_DA", "y_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "y (inelasticity) PDF Comparison",
        "y",
        "PDF",
        "figs/inclusive/histos/y_methods_pdf.png",
        false,
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // =================================================================
    // NEW: W² distributions
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"W2_truth", "W2_EM", "W2_DA", "W2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "W^{2} Reconstruction Methods",
        "W^{2} [GeV^{2}]",
        "# of events",
        "figs/inclusive/histos/w2_methods_hist.png",
        true,
        true
    );
    plot_ptr->SetRangeX(10.0, 1.0e4);
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"W2_truth", "W2_EM", "W2_DA", "W2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "W^{2} PDF Comparison",
        "W^{2} [GeV^{2}]",
        "PDF",
        "figs/inclusive/histos/w2_methods_pdf.png",
        true,
        true,
        true
    );
    plot_ptr->SetRangeX(10.0, 1.0e4);
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // =================================================================
    // NEW: Scattered electron leptonic quantities
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"Ep_e_truth", "Ep_e"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Scattered Electron Energy",
        "E'_{e} [GeV]",
        "Counts",
        "figs/inclusive/histos/electron_energy.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"Ep_e_truth", "Ep_e"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Scattered Electron Energy",
        "E'_{e} [GeV]",
        "Counts",
        "figs/inclusive/histos/electron_energy_logy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"phi_e_truth", "phi_e"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Scattered Electron Azimuthal Angle",
        "#phi_{e} [rad]",
        "Counts",
        "figs/inclusive/histos/electron_phi.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.4, 0.2, 0.6, 0.4);
    plots.push_back(plot_ptr);

    // =================================================================
    // NEW: Scattered electron p_{T} distributions
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"pT_e_truth", "pT_e"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Scattered Electron p_{T}",
        "p_{T}^{e} [GeV]",
        "Counts",
        "figs/inclusive/histos/electron_pt.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"pT_e_truth", "pT_e"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Scattered Electron p_{T}",
        "p_{T}^{e} [GeV]",
        "Counts",
        "figs/inclusive/histos/electron_pt_logy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // =================================================================
    // Unbinned correlation plots (Reco vs MC) - one per method
    // =================================================================
#if 0
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_Q2_EM"},
        {"EM"},
        {kRed},
        {20},
        "Q^{2} Correlation (EM, Unbinned)",
        "Q^{2}_{truth} [GeV^{2}]",
        "Q^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/q2_corr_unbinned_em.png",
        {1.0, 300.0},
        {1.0, 300.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_Q2_DA"},
        {"DA"},
        {kBlue},
        {21},
        "Q^{2} Correlation (DA, Unbinned)",
        "Q^{2}_{truth} [GeV^{2}]",
        "Q^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/q2_corr_unbinned_da.png",
        {1.0, 300.0},
        {1.0, 300.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_Q2_Sigma"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "Q^{2} Correlation (#Sigma, Unbinned)",
        "Q^{2}_{truth} [GeV^{2}]",
        "Q^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/q2_corr_unbinned_sigma.png",
        {1.0, 300.0},
        {1.0, 300.0},
        true,
        true
    ));

#if 0
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_x_EM"},
        {"EM"},
        {kRed},
        {20},
        "x_{Bj} Correlation (EM, Unbinned)",
        "x_{truth}",
        "x_{reco}",
        "figs/inclusive/histos/xbj_corr_unbinned_em.png",
        {1e-4, 1.0},
        {1e-4, 1.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_x_DA"},
        {"DA"},
        {kBlue},
        {21},
        "x_{Bj} Correlation (DA, Unbinned)",
        "x_{truth}",
        "x_{reco}",
        "figs/inclusive/histos/xbj_corr_unbinned_da.png",
        {1e-4, 1.0},
        {1e-4, 1.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_x_Sigma"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "x_{Bj} Correlation (#Sigma, Unbinned)",
        "x_{truth}",
        "x_{reco}",
        "figs/inclusive/histos/xbj_corr_unbinned_sigma.png",
        {1e-4, 1.0},
        {1e-4, 1.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_y_EM"},
        {"EM"},
        {kRed},
        {20},
        "y Correlation (EM, Unbinned)",
        "y_{truth}",
        "y_{reco}",
        "figs/inclusive/histos/y_corr_unbinned_em.png",
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_y_DA"},
        {"DA"},
        {kBlue},
        {21},
        "y Correlation (DA, Unbinned)",
        "y_{truth}",
        "y_{reco}",
        "figs/inclusive/histos/y_corr_unbinned_da.png",
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_y_Sigma"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "y Correlation (#Sigma, Unbinned)",
        "y_{truth}",
        "y_{reco}",
        "figs/inclusive/histos/y_corr_unbinned_sigma.png",
        {0.0, 1.0},
        {0.0, 1.0}
    ));
#endif

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_W2_EM"},
        {"EM"},
        {kRed},
        {20},
        "W^{2} Correlation (EM, Unbinned)",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/w2_corr_unbinned_em.png",
        {10.0, 1.0e4},
        {10.0, 1.0e4},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_W2_DA"},
        {"DA"},
        {kBlue},
        {21},
        "W^{2} Correlation (DA, Unbinned)",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/w2_corr_unbinned_da.png",
        {10.0, 1.0e4},
        {10.0, 1.0e4},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_W2_Sigma"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "W^{2} Correlation (#Sigma, Unbinned)",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/w2_corr_unbinned_sigma.png",
        {10.0, 1.0e4},
        {10.0, 1.0e4},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_W2_EM", "g_W2_DA", "g_W2_Sigma"},
        {"EM", "DA", "Sigma"},
        {kRed, kBlue, kGreen + 2},
        {20, 21, 22},
        "W^{2} Correlation (EM/DA/#Sigma, Unbinned)",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/histos/w2_corr_unbinned_all.png",
        {10.0, 1.0e4},
        {10.0, 1.0e4},
        true,
        true
    ));

#if 0
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_Ep_e"},
        {"Electron"},
        {kBlack},
        {20},
        "E'_{e} Correlation (Unbinned)",
        "E'_{e,truth} [GeV]",
        "E'_{e,reco} [GeV]",
        "figs/inclusive/histos/electron_energy_corr_unbinned.png",
        {0.0, 20.0},
        {0.0, 20.0}
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_phi_e"},
        {"Electron"},
        {kBlack},
        {20},
        "#phi_{e} Correlation (Unbinned)",
        "#phi_{e,truth} [rad]",
        "#phi_{e,reco} [rad]",
        "figs/inclusive/histos/electron_phi_corr_unbinned.png",
        {-3.2, 3.2},
        {-3.2, 3.2}
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_pT_e"},
        {"Electron"},
        {kBlack},
        {20},
        "p_{T}^{e} Correlation (Unbinned)",
        "p_{T,truth}^{e} [GeV]",
        "p_{T,reco}^{e} [GeV]",
        "figs/inclusive/histos/electron_pt_corr_unbinned.png",
        {0.0, 10.0},
        {0.0, 10.0}
    ));
#endif
#endif

    // =================================================================
    // NEW: Unfolding response matrices (raw and normalized)
    // =================================================================
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_Q2",
        "Q^{2}_{truth} [GeV^{2}]",
        "Q^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_q2_electron_raw.png",
        true,
        true,
        {1.0, 300},
        {1.0, 300}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_x",
        "x_{truth}",
        "x_{reco}",
        "figs/inclusive/response/response_xbj_electron_raw.png",
        true,
        true,
        {1e-4, 1.0},
        {1e-4, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_y",
        "y_{truth}",
        "y_{reco}",
        "figs/inclusive/response/response_y_electron_raw.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_Q2_rowNorm",
        "Q^{2}_{truth} [GeV^{2}]",
        "Q^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_q2_electron_rowNorm.png",
        true,
        true,
        {1.0, 300},
        {1.0, 300}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_Q2_colNorm",
        "Q^{2}_{truth} [GeV^{2}]",
        "Q^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_q2_electron_colNorm.png",
        true,
        true,
        {1.0, 300},
        {1.0, 300}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_x_rowNorm",
        "x_{truth}",
        "x_{reco}",
        "figs/inclusive/response/response_xbj_electron_rowNorm.png",
        true,
        true,
        {1e-4, 1.0},
        {1e-4, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_x_colNorm",
        "x_{truth}",
        "x_{reco}",
        "figs/inclusive/response/response_xbj_electron_colNorm.png",
        true,
        true,
        {1e-4, 1.0},
        {1e-4, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_y_rowNorm",
        "y_{truth}",
        "y_{reco}",
        "figs/inclusive/response/response_y_electron_rowNorm.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_y_colNorm",
        "y_{truth}",
        "y_{reco}",
        "figs/inclusive/response/response_y_electron_colNorm.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_EM",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_em.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_DA",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_da.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_Sigma",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_sigma.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_EM_rowNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_em_rowNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_EM_colNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_em_colNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_DA_rowNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_da_rowNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_DA_colNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_da_colNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_Sigma_rowNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_sigma_rowNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_Sigma_colNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_sigma_colNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));

    // =================================================================
    // Overall resolution plots (Q2/x/y)
    // =================================================================
    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_EM",
        "#frac{Q^{2}_{EM} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Counts",
        -0.01, 0.01,
        "figs/inclusive/resolution/q2_relres_em.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_DA",
        "#frac{Q^{2}_{DA} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Counts",
        -0.005, 0.02,
        "figs/inclusive/resolution/q2_relres_da.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_Sigma",
        "#frac{Q^{2}_{#Sigma} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Counts",
        -0.01, 0.01,
        "figs/inclusive/resolution/q2_relres_sigma.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_EM",
        "#frac{x_{EM} - x_{MC}}{ x_{MC}}",
        "Counts",
        -0.025, 0.02,
        "figs/inclusive/resolution/xbj_relres_em.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_DA",
        "#frac{x_{DA} - x_{MC}}{X_{MC}}",
        "Counts",
        0., 0.,
        "figs/inclusive/resolution/xbj_relres_da.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_Sigma",
        "#frac{x_{#Sigma} - x_{MC}}{X_{MC}}",
        "Counts",
        0., 0.,
        "figs/inclusive/resolution/xbj_relres_sigma.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_EM",
        "#frac{y_{EM} - y_{MC}}{ y_{MC}}",
        "Counts",
        -0.009, 0.009,
        "figs/inclusive/resolution/y_relres_em.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_DA",
        "#frac{y_{DA} - y_{MC}}{y_{MC}}",
        "Counts",
        0., 0.,
        "figs/inclusive/resolution/y_relres_da.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_Sigma",
        "#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "Counts",
        0., 0.,
        "figs/inclusive/resolution/y_relres_sigma.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_EM",
        "#frac{W^{2}_{EM} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Counts",
        -0.5, 0.5,
        "figs/inclusive/resolution/w2_relres_em.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_DA",
        "#frac{W^{2}_{DA} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Counts",
        -0.5, 0.5,
        "figs/inclusive/resolution/w2_relres_da.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_Sigma",
        "#frac{W^{2}_{#Sigma} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Counts",
        -0.5, 0.5,
        "figs/inclusive/resolution/w2_relres_sigma.png"
    ));

    // =================================================================
    // Binned resolution plots (Q2/x/y)
    // =================================================================
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_EM",
        "Relative bin by bin resolution (EM);Q^{2}_{MC};#frac{Q^{2}_{EM} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{EM}",
        "",
        {
         {-0.015, 0.015}, {-0.015, 0.015}, {-0.015, 0.015}, {-0.015, 0.015}, {-0.015, 0.015},
         {-0.015, 0.015}, {-0.015, 0.015}, {-0.014, 0.015}, {-0.025, 0.025}, {-0.01, 0.012},
         {-0.027, 0.028}, {-0.018, 0.02}, {-0.022, 0.02}, {-0.02, 0.015}, {-0.018, 0.02},
         {-0.02, 0.015}, {-0.02, 0.017}, {-0.017, 0.02}, {-0.02, 0.02}, {-0.04, 0.04},
         {-0.025, 0.03}, {-0.015, 0.025}, {-0.05, 0.06}
        },
        "figs/inclusive/resolution/q2_relres_binned_em.png",
        "figs/inclusive/resolution/profile/q2_relres_binned_em",
        std::make_pair(5.0, 200),
        true
    ));

    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_DA",
        "Relative bin by bin resolution (DA);Q^{2}_{MC};#frac{Q^{2}_{DA} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{DA}",
        "",
        {
          {-0.008, 0.015}, {-0.008, 0.015}, {-0.008, 0.015}, {-0.008, 0.02}, {-0.01, 0.025},
          {-0.0, 0.0}, {-0, 0}, {-0., 0.}, {-0., 0.}, {-0., 0.},
          {-0.015, 0.03}, {-0.009, 0.02}, {-0.01, 0.02}, {-0.01, 0.02}, {-0.01, 0.02},
          {-0.01, 0.02}, {-0.004, 0.02}, {-0.017, 0.027}, {-0.025, 0.03}, {-0.08, 0.08},
          {-0.05, 0.06}, {-0.05, 0.065}, {-0.05, 0.06}
        },
        "figs/inclusive/resolution/q2_relres_binned_da.png",
        "figs/inclusive/resolution/profile/q2_relres_binned_da",
        std::make_pair(5.0, 200),
        true
    ));

    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_Sigma",
        ";Q^{2}_{MC};#frac{Q^{2}_{#Sigma} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{#Sigma}",
        "",
        {
         {-0.07, 0.02}, {-0.07, 0.02}, {-0.07, 0.01}, {-0.07, 0.01}, {-0.07, 0.01},
         {-0.0, 0.0}, {-0.015, 0.015}, {-0.014, 0.015}, {-0.025, 0.025}, {-0.01, 0.012},
         {-0.027, 0.028}, {-0.018, 0.02}, {-0.022, 0.02}, {-0.02, 0.015}, {-0.018, 0.02},
         {-0.02, 0.015}, {-0.02, 0.017}, {-0.017, 0.02}, {-0.02, 0.02}, {-0.04, 0.04},
         {-0.025, 0.03}, {-0.015, 0.025}, {-0.05, 0.06}
        },
        "figs/inclusive/resolution/q2_relres_binned_sigma.png",
        "figs/inclusive/resolution/profile/q2_relres_binned_sigma",
        std::make_pair(5.0, 200),
        true
    ));

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_EM",
        ";x_{MC};#frac{x_{EM} - x_{MC}}{x_{MC}}",
        "x_{EM}",
        "",
        {
         {-0.0, 0.0}, {-0.022, 0.02}, {-0.02, 0.02}, {-0.02, 0.02}, {-0.02, 0.02},
         {-0.015, 0.015}, {-0.015, 0.015}, {-0.01, 0.01}, {-0.01, 0.01}, {-0.015, 0.015},
         {-0.015, 0.015}, {-0.015, 0.015}, {-0.015, 0.015}, {-0.02, 0.02}, {-0.015, 0.015},
         {-0.015, 0.015}, {-0.015, 0.015}, {0, 0}, {0, 0}, {0, 0},
         {0, 0}, {0, 0}, {0, 0}
        },
        "figs/inclusive/resolution/xbj_relres_binned_em.png",
        "figs/inclusive/resolution/profile/xbj_relres_binned_em",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_DA",
        ";x_{MC};#frac{x_{DA} - x_{MC}}{x_{MC}}",
        "x_{DA}",
        "",
        {
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0}
        },
        "figs/inclusive/resolution/xbj_relres_binned_da.png",
        "figs/inclusive/resolution/profile/xbj_relres_binned_da",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_Sigma",
        ";x_{MC};#frac{x_{#Sigma} - x_{MC}}{x_{MC}}",
        "x_{#Sigma}",
        "",
        {
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0}
        },
        "figs/inclusive/resolution/xbj_relres_binned_sigma.png",
        "figs/inclusive/resolution/profile/xbj_relres_binned_sigma",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_EM",
        ";y_{MC};#frac{y_{EM} - y_{MC}}{y_{MC}}",
        "y_{EM}",
        "",
        {
            {0,0},{-0.09,0.08},{-0.05,0.05},{-0.04,0.04},{-0.02,0.02},
            {-0.01,0.01},{-0.01,0.01},{-0.01,0.01},{-0.01,0.01},{0,0}
        },
        "figs/inclusive/resolution/y_relres_binned_em.png",
        "figs/inclusive/resolution/profile/y_relres_binned_em",
        std::make_pair(0.0, 1.0),
        false
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_DA",
        ";y_{MC};#frac{y_{DA} - y_{MC}}{y_{MC}}",
        "y_{DA}",
        "",
        {
            {-0.2,0.05},{-0.2,0.05},{-0.2,0.05},{-0.18,0.1},{-0.15,0.05},
            {-0.12,0.1},{-0.1,0.06},{-0.06,0.04},{-0.06,0.04},{-0.05,0.03}
        },
        "figs/inclusive/resolution/y_relres_binned_da.png",
        "figs/inclusive/resolution/profile/y_relres_binned_da",
        std::make_pair(0.0, 1.0),
        false
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_Sigma",
        ";y_{MC};#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "y_{#Sigma}",
        "",
        {
            {-0.2,0.05},{-0.15,0.05},{-0.1,0.05},{-0.05,0.02},{-0.02,0.02},
            {0,0},{0,0},{0,0},{0,0},{0,0}
        },
        "figs/inclusive/resolution/y_relres_binned_sigma.png",
        "figs/inclusive/resolution/profile/y_relres_binned_sigma",
        std::make_pair(0.0, 1.0),
        false
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "W2_RelRes_binned_EM",
        "W^{2} relative bin by bin resolution (EM);W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{EM} - W^{2}_{MC}}{W^{2}_{MC}}",
        "W^{2}_{EM}",
        "",
        {
            {0,0}, {-0.4,0.4}, {-0.2,0.2}, {-0.15,0.15}, {-0.1,0.1},
            {-0.05,0.05}, {-0.04,0.04}
        },
        "figs/inclusive/resolution/w2_relres_binned_em.png",
        "figs/inclusive/resolution/profile/w2_relres_binned_em",
        std::make_pair(10.0, 1.0e4),
        true
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "W2_RelRes_binned_DA",
        "W^{2} relative bin by bin resolution (DA);W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{DA} - W^{2}_{MC}}{W^{2}_{MC}}",
        "W^{2}_{DA}",
        "",
        {
            {-0.3,0.2}, {-0.25,0.15}, {-0.25,0.15}, {-0.3,0.1}, {-0.25,0.1},
            {-0.2,0.1}, {-0.1,0.1}
        },
        "figs/inclusive/resolution/w2_relres_binned_da.png",
        "figs/inclusive/resolution/profile/w2_relres_binned_da",
        std::make_pair(10.0, 1.0e4),
        true
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "W2_RelRes_binned_Sigma",
        "W^{2} relative bin by bin resolution (#Sigma);W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{#Sigma} - W^{2}_{MC}}{W^{2}_{MC}}",
        "W^{2}_{#Sigma}",
        "",
        {
            {-0,0}, {-0.25,0.06}, {-0.18,0.02}, {-0.,0.}, {-0.25,0.02},
            {-0.15,0.02}, {0,0}
        },
        "figs/inclusive/resolution/w2_relres_binned_sigma.png",
        "figs/inclusive/resolution/profile/w2_relres_binned_sigma",
        std::make_pair(10.0, 1.0e4),
        true
    );
    plots.push_back(binned_plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xL_MC", "xL_B0", "xL_RP"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "pe", "pe"},
        "x_{L} Distributions",
        "x_{L}",
        "Counts",
        "figs/diffractive/histos/xL_distributions.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

#if 0
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xL_B0", "g_xL_RP"},
        {"B0", "RP"},
        {kRed, kBlue},
        {20, 24},
        "x_{L} Correlation (Unbinned)",
        "Truth x_{L}",
        "Reco x_{L}",
        "figs/diffractive/histos/xL_corr_unbinned.png",
        {0.75, 1.05},
        {0.75, 1.05}
    ));
#endif

    // =================================================================
    // Diffractive: x_{pom} (definition) - split by B0/RP
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"xpom_truth_all", "xpom_reco_EM_all", "xpom_reco_DA_all", "xpom_reco_Sigma_all"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{pom} Distributions (All)",
        "x_{pom}",
        "Counts",
        "figs/diffractive/histos/xpom_def_comparison_all.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-4, 0.3);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_truth_B0", "xpom_reco_EM_B0", "xpom_reco_DA_B0", "xpom_reco_Sigma_B0"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{pom} Distributions (B0)",
        "x_{pom}",
        "Counts",
        "figs/diffractive/histos/xpom_def_comparison_b0.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-4, 0.3);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_truth_RP", "xpom_reco_EM_RP", "xpom_reco_DA_RP", "xpom_reco_Sigma_RP"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{pom} Distributions (RP)",
        "x_{pom}",
        "Counts",
        "figs/diffractive/histos/xpom_def_comparison_rp.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-4, 0.3);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_EM_B0"},
        {"EM"},
        {kRed},
        {20},
        "x_{pom} Correlation (EM, B0)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/histos/xpom_corr_unbinned_em_b0.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_DA_B0"},
        {"DA"},
        {kBlue},
        {21},
        "x_{pom} Correlation (DA, B0)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/histos/xpom_corr_unbinned_da_b0.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_Sigma_B0"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "x_{pom} Correlation (#Sigma, B0)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/histos/xpom_corr_unbinned_sigma_b0.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_EM_RP"},
        {"EM"},
        {kRed},
        {20},
        "x_{pom} Correlation (EM, RP)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/histos/xpom_corr_unbinned_em_rp.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_DA_RP"},
        {"DA"},
        {kBlue},
        {21},
        "x_{pom} Correlation (DA, RP)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/histos/xpom_corr_unbinned_da_rp.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_Sigma_RP"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "x_{pom} Correlation (#Sigma, RP)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/histos/xpom_corr_unbinned_sigma_rp.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_xpom_EM_B0",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/response/response_xpom_em_b0.png",
        true,
        true,
        {1e-4, 0.3},
        {1e-4, 0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_xpom_EM_RP",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/response/response_xpom_em_rp.png",
        true,
        true,
        {1e-4, 0.3},
        {1e-4, 0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_xpom_DA_B0",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/response/response_xpom_da_b0.png",
        true,
        true,
        {1e-4, 0.3},
        {1e-4, 0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_xpom_DA_RP",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/response/response_xpom_da_rp.png",
        true,
        true,
        {1e-4, 0.3},
        {1e-4, 0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_xpom_Sigma_B0",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/response/response_xpom_sigma_b0.png",
        true,
        true,
        {1e-4, 0.3},
        {1e-4, 0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_xpom_Sigma_RP",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/response/response_xpom_sigma_rp.png",
        true,
        true,
        {1e-4, 0.3},
        {1e-4, 0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_beta_EM_B0",
        "Truth #beta",
        "Reco #beta",
        "figs/diffractive/response/response_beta_em_b0.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_beta_EM_RP",
        "Truth #beta",
        "Reco #beta",
        "figs/diffractive/response/response_beta_em_rp.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_beta_DA_B0",
        "Truth #beta",
        "Reco #beta",
        "figs/diffractive/response/response_beta_da_b0.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_beta_DA_RP",
        "Truth #beta",
        "Reco #beta",
        "figs/diffractive/response/response_beta_da_rp.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_beta_Sigma_B0",
        "Truth #beta",
        "Reco #beta",
        "figs/diffractive/response/response_beta_sigma_b0.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_beta_Sigma_RP",
        "Truth #beta",
        "Reco #beta",
        "figs/diffractive/response/response_beta_sigma_rp.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    // =================================================================
    // Diffractive: beta = x_{Bj}/x_{pom} (B0/RP)
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"beta_truth_all", "beta_reco_EM_all", "beta_reco_DA_all", "beta_reco_Sigma_all"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "#beta Distributions (All)",
        "#beta",
        "Counts",
        "figs/diffractive/histos/beta_comparison_all.png",
        false,
        false
    );
    plot_ptr->SetRangeX(0.0, 1.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"beta_truth_B0", "beta_reco_EM_B0", "beta_reco_DA_B0", "beta_reco_Sigma_B0"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "#beta Distributions (B0)",
        "#beta",
        "Counts",
        "figs/diffractive/histos/beta_comparison_b0.png",
        false,
        false
    );
    plot_ptr->SetRangeX(0.0, 1.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"beta_truth_RP", "beta_reco_EM_RP", "beta_reco_DA_RP", "beta_reco_Sigma_RP"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "#beta Distributions (RP)",
        "#beta",
        "Counts",
        "figs/diffractive/histos/beta_comparison_rp.png",
        false,
        false
    );
    plot_ptr->SetRangeX(0.0, 1.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_EM_B0"},
        {"EM"},
        {kRed},
        {20},
        "#beta Correlation (EM, B0)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/histos/beta_corr_unbinned_em_b0.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_DA_B0"},
        {"DA"},
        {kBlue},
        {21},
        "#beta Correlation (DA, B0)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/histos/beta_corr_unbinned_da_b0.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_Sigma_B0"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "#beta Correlation (#Sigma, B0)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/histos/beta_corr_unbinned_sigma_b0.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_EM_RP"},
        {"EM"},
        {kRed},
        {20},
        "#beta Correlation (EM, RP)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/histos/beta_corr_unbinned_em_rp.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_DA_RP"},
        {"DA"},
        {kBlue},
        {21},
        "#beta Correlation (DA, RP)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/histos/beta_corr_unbinned_da_rp.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_Sigma_RP"},
        {"Sigma"},
        {kGreen + 2},
        {22},
        "#beta Correlation (#Sigma, RP)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/histos/beta_corr_unbinned_sigma_rp.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    // =================================================================
    // Diffractive: M_X^2 plots
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"MX2_truth", "MX2_reco"},
        {"MC Truth", "Reco"},
        {"hist", "pe"},
        "M_{X}^{2} Distributions",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/diffractive/histos/MX2_comparison.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-3, 1000.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"MX2_truth", "MX2_reco"},
        {"MC Truth", "Reco"},
        {"hist", "pe"},
        "M_{X}^{2} Distributions",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/diffractive/histos/MX2_comparison_logy.png",
        true,
        true
    );
    plot_ptr->SetRangeX(1e-3, 1000.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_MX2"},
        {"Reco vs Truth"},
        {kBlack},
        {20},
        "M_{X}^{2} Correlation (Unbinned)",
        "M_{X,truth}^{2} [GeV^{2}]",
        "M_{X,reco}^{2} [GeV^{2}]",
        "figs/diffractive/histos/MX2_corr_unbinned.png",
        {1e-3, 1000.0},
        {1e-3, 1000.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsRelRes(
        "MX2_RelRes",
        "#frac{M_{X,reco}^{2} - M_{X,truth}^{2}}{M_{X,truth}^{2}}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/MX2_resolution.png"
    ));

#if 0
    plots.push_back(new PlotOptions2D(
        "MX2_corr",
        "M_{X,truth}^{2} [GeV^{2}]",
        "M_{X,reco}^{2} [GeV^{2}]",
        "figs/diffractive/histos/MX2_truth_vs_reco.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 1000.0}
    ));
#endif

#if 0
    plots.push_back(new PlotOptions2D(
        "MX2_t_truth",
        "M_{X}^{2} [GeV^{2}]",
        "|t| [GeV^{2}]",
        "figs/diffractive/histos/MX2_vs_t_truth.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptions2D(
        "MX2_t_B0",
        "M_{X}^{2} [GeV^{2}]",
        "|t| [GeV^{2}]",
        "figs/diffractive/histos/MX2_vs_t_b0.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptions2D(
        "MX2_t_RP",
        "M_{X}^{2} [GeV^{2}]",
        "|t| [GeV^{2}]",
        "figs/diffractive/histos/MX2_vs_t_rp.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 2.0}
    ));
#endif

    // =================================================================
    // Diffractive: Mandelstam t plots
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| Distributions",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/diffractive/histos/t_distributions.png",
        true,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| Distributions",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/diffractive/histos/t_distributions_logy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"dsigma_dt_MC", "dsigma_dt_B0", "dsigma_dt_RP", "dsigma_dt_Sum"},
        {"MC Truth", "B0 Reco", "RP Reco", "B0+RP Sum"},
        {"hist", "pe", "pe", "pe"},
        "d#sigma/dt",
        "|t| [GeV^{2}]",
        "d#sigma/dt [nb/GeV^{2}]",
        "figs/diffractive/histos/dsigma_dt.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"theta_MC", "theta_B0", "theta_RP"},
        {"MC Truth", "B0 Reco", "RP Reco (#theta #leq 5 mrad)"},
        {"hist", "pe", "pe"},
        "Proton Scattering Angles",
        "#theta [mrad]",
        "Counts",
        "figs/diffractive/histos/theta_distributions.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

#if 0
    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_t_B0", "g_t_RP"},
        {"B0", "RP"},
        {kRed, kBlue},
        {20, 24},
        "|t| Correlation (Unbinned)",
        "Truth |t| [GeV^{2}]",
        "Reco |t| [GeV^{2}]",
        "figs/diffractive/histos/t_corr_unbinned.png",
        {1e-3, 2.0},
        {1e-3, 2.0},
        true,
        true
    ));
#endif

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}]",
        "figs/diffractive/response/response_t_b0.png",
        true,
        true,
        {1e-3, 2.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_RP",
        "Truth |t| [GeV^{2}]",
        "RP Reco |t| [GeV^{2}]",
        "figs/diffractive/response/response_t_rp.png",
        true,
        true,
        {1e-3, 2.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptionsRelRes(
        "t_res_B0",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -999., -999.,
        "figs/diffractive/resolution/t_res_b0.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "t_res_RP",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -999., -999.,
        "figs/diffractive/resolution/t_res_rp.png"
    ));

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "t_RelRes_binned_B0",
        ";|t|_{truth} [GeV^{2}];#frac{|t|_{reco}-|t|_{truth}}{|t|_{truth}}",
        "|t|_{B0}",
        "",
        {},
        "figs/diffractive/resolution/t_relres_binned_b0.png",
        "diffractive/resolution/profile/t_relres_binned_b0",
        std::make_pair(1e-3, 2.0),
        true
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "t_RelRes_binned_RP",
        ";|t|_{truth} [GeV^{2}];#frac{|t|_{reco}-|t|_{truth}}{|t|_{truth}}",
        "|t|_{RP}",
        "",
        {},
        "figs/diffractive/resolution/t_relres_binned_rp.png",
        "diffractive/resolution/profile/t_relres_binned_rp",
        std::make_pair(1e-3, 2.0),
        true
    );
    plots.push_back(binned_plot_ptr);

    // =================================================================
    // Legacy plot set from combined/t/xy plotters
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "E-p_{z} Distribution",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/distributions/EPz_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "E-p_{z} Distribution",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/distributions/EPz_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Counts",
        "figs/distributions/eta_max_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Counts",
        "figs/distributions/eta_max_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"MX2_truth", "MX2_reco"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "M_{X}^{2} Distribution",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/distributions/MX2_distribution.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-3, 1000.0);
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"MX2_truth", "MX2_reco"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "M_{X}^{2} Distribution",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/distributions/MX2_distribution_logY.png",
        true,
        true
    );
    plot_ptr->SetRangeX(1e-3, 1000.0);
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| PDF Comparison",
        "|t| [GeV^{2}]",
        "PDF",
        "figs/distributions/t_pdf_comparison.png",
        true,
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsRelRes(
        "xL_res_B0",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/xL_resolution_B0.png"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "xL_res_RP",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/xL_resolution_RP.png"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "xpom_res_B0",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/xpom_resolution_B0.png"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "xpom_res_RP",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/xpom_resolution_RP.png"
    ));
    plots.push_back(new PlotOptions1D(
        {"beta_res_B0"},
        {"B0 Reco"},
        {"hist"},
        "B0 #beta Resolution",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Counts",
        "figs/resolutions/simple/beta_resolution_B0.png",
        false,
        false
    ));
    plots.push_back(new PlotOptions1D(
        {"beta_res_RP"},
        {"RP Reco"},
        {"hist"},
        "RP #beta Resolution",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Counts",
        "figs/resolutions/simple/beta_resolution_RP.png",
        false,
        false
    ));

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xL_RelRes_binned_B0",
        "B0 x_{L} Resolution vs Truth x_{L}",
        "x_{L,truth}",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/xL_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/xL_B0",
        std::make_pair(0.75, 1.05),
        false
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xL_RelRes_binned_RP",
        "RP x_{L} Resolution vs Truth x_{L}",
        "x_{L,truth}",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/xL_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/xL_RP",
        std::make_pair(0.75, 1.05),
        false
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_B0",
        "B0 x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/xpom_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/xpom_B0",
        std::make_pair(1e-4, 0.4),
        true
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_RP",
        "RP x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/xpom_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/xpom_RP",
        std::make_pair(1e-4, 0.4),
        true
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_B0",
        "B0 #beta Resolution vs Truth #beta",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/beta_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/beta_B0",
        std::make_pair(0.0, 1.0),
        false
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_RP",
        "RP #beta Resolution vs Truth #beta",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/beta_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/beta_RP",
        std::make_pair(0.0, 1.0),
        false
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "MX2_RelRes_binned",
        "M_{X}^{2} Resolution vs Truth M_{X}^{2}",
        "M_{X,truth}^{2} [GeV^{2}]",
        "(M_{X,reco}^{2} - M_{X,truth}^{2})/M_{X,truth}^{2}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/MX2_resolution_binned.png",
        "figs/resolutions/binned/bins/MX2",
        std::make_pair(1e-3, 1000.0),
        true
    );
    plots.push_back(binned_plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_EM",
        "Q^{2} (true) [GeV^{2}]",
        "Q^{2} (EM) [GeV^{2}]",
        "figs/response_matrices/response_matrix_Q2_EM.png",
        true,
        true,
        {1.0, 300.0},
        {1.0, 300.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_DA",
        "Q^{2} (true) [GeV^{2}]",
        "Q^{2} (DA) [GeV^{2}]",
        "figs/response_matrices/response_matrix_Q2_DA.png",
        true,
        true,
        {1.0, 300.0},
        {1.0, 300.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_Sigma",
        "Q^{2} (true) [GeV^{2}]",
        "Q^{2} (Sigma) [GeV^{2}]",
        "figs/response_matrices/response_matrix_Q2_Esigma.png",
        true,
        true,
        {1.0, 300.0},
        {1.0, 300.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_EM",
        "x_{Bj} (true)",
        "x_{Bj} (EM)",
        "figs/response_matrices/response_matrix_x_EM.png",
        true,
        true,
        {1e-3, 0.3},
        {1e-3, 0.3}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_DA",
        "x_{Bj} (true)",
        "x_{Bj} (DA)",
        "figs/response_matrices/response_matrix_x_DA.png",
        true,
        true,
        {1e-3, 0.3},
        {1e-3, 0.3}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_Sigma",
        "x_{Bj} (true)",
        "x_{Bj} (Sigma)",
        "figs/response_matrices/response_matrix_x_Sigma.png",
        true,
        true,
        {1e-3, 0.3},
        {1e-3, 0.3}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_EM",
        "y (true)",
        "y (EM)",
        "figs/response_matrices/response_matrix_y_EM.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_DA",
        "y (true)",
        "y (DA)",
        "figs/response_matrices/response_matrix_y_DA.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_Sigma",
        "y (true)",
        "y (Sigma)",
        "figs/response_matrices/response_matrix_y_Sigma.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}]",
        "figs/response_matrices/response_matrix_t_B0.png",
        true,
        true,
        {1e-3, 2.0},
        {1e-3, 2.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_RP",
        "Truth |t| [GeV^{2}]",
        "RP Reco |t| [GeV^{2}]",
        "figs/response_matrices/response_matrix_t_RP.png",
        true,
        true,
        {1e-3, 2.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_B0",
        "Truth x_{L}",
        "B0 Reco x_{L}",
        "figs/response_matrices/response_matrix_xL_B0.png",
        false,
        false,
        {0.75, 1.05},
        {0.75, 1.05}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_RP",
        "Truth x_{L}",
        "RP Reco x_{L}",
        "figs/response_matrices/response_matrix_xL_RP.png",
        false,
        false,
        {0.75, 1.05},
        {0.75, 1.05}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0",
        "Truth x_{pom} (1-x_{L})",
        "B0 Reco x_{pom} (1-x_{L})",
        "figs/response_matrices/response_matrix_xpom_B0.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_RP",
        "Truth x_{pom} (1-x_{L})",
        "RP Reco x_{pom} (1-x_{L})",
        "figs/response_matrices/response_matrix_xpom_RP.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_B0",
        "Truth #beta",
        "B0 Reco #beta",
        "figs/response_matrices/response_matrix_beta_B0.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_RP",
        "Truth #beta",
        "RP Reco #beta",
        "figs/response_matrices/response_matrix_beta_RP.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "MX2_corr",
        "M_{X,truth}^{2} [GeV^{2}]",
        "M_{X,reco}^{2} [GeV^{2}]",
        "figs/response_matrices/response_matrix_MX2.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 1000.0}
    ));

    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_def_MC"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "MC Truth x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_MC_logxy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_B0", "xpom_def_B0"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "B0 Reco x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_B0_logxy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_RP", "xpom_def_RP"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "RP Reco x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_RP_logxy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_B0", "xpom_RP"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "hist", "hist"},
        "x_{pom} Comparison (All)",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_all_logxy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.45, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_MC",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_MC.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_B0",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_B0.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_RP",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_RP.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plot_ptr = new PlotOptions1D(
        {"theta_all_TS", "theta_B0"},
        {"All Truth-Seeded Protons", "B0 Accepted Protons"},
        {"hist", "hist"},
        "Proton Scattering Angle Distribution",
        "#theta [mrad]",
        "Counts",
        "figs/distributions/theta_comparison_B0_acceptance.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"theta_all_TS", "theta_B0"},
        {"All Truth-Seeded Protons", "B0 Accepted Protons"},
        {"hist", "hist"},
        "Proton Scattering Angle Distribution",
        "#theta [mrad]",
        "Counts",
        "figs/distributions/theta_comparison_B0_acceptance_logxy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"beta_MC", "beta_B0", "beta_RP"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "E1", "E1"},
        "#beta = x_{Bj} / x_{pom} Distributions",
        "#beta",
        "Counts",
        "figs/distributions/beta_distributions_logy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_Q2",
        "Q^{2} [GeV^{2}]",
        "#beta",
        "figs/distributions/beta_vs_Q2.png",
        true,
        false,
        {1.0, 300.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_xpom",
        "x_{pom}",
        "#beta",
        "figs/distributions/beta_vs_xpom.png",
        true,
        false,
        {1e-4, 0.4},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_t",
        "|t| [GeV^{2}]",
        "#beta",
        "figs/distributions/beta_vs_t.png",
        false,
        false,
        {1e-3, 2.0},
        {0.0, 1.0}
    ));

    // =================================================================
    // =================================================================
    // Execute all standard plots
    // =================================================================
    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    for (const auto& plot : plots) {
        delete plot;
    }

    // =================================================================
    // Relative resolution vs global bin index (k)
    // =================================================================
    PlotRelResVsK(inputFile,
                  {"Q2_RelRes_vs_k_EM", "Q2_RelRes_vs_k_DA", "Q2_RelRes_vs_k_Sigma"},
                  {"EM", "DA", "Sigma"},
                  "Q^{2} relative resolution vs k",
                  "figs/resolutions/binned/q2_relres_vs_k.png");

    PlotRelResVsK(inputFile,
                  {"xpom_RelRes_vs_k_EM_B0", "xpom_RelRes_vs_k_DA_B0", "xpom_RelRes_vs_k_Sigma_B0"},
                  {"EM", "DA", "Sigma"},
                  "x_{pom} relative resolution vs k (B0)",
                  "figs/resolutions/binned/xpom_relres_vs_k_b0.png");
    PlotRelResVsK(inputFile,
                  {"xpom_RelRes_vs_k_EM_RP", "xpom_RelRes_vs_k_DA_RP", "xpom_RelRes_vs_k_Sigma_RP"},
                  {"EM", "DA", "Sigma"},
                  "x_{pom} relative resolution vs k (RP)",
                  "figs/resolutions/binned/xpom_relres_vs_k_rp.png");

    PlotRelResVsK(inputFile,
                  {"beta_RelRes_vs_k_EM_B0", "beta_RelRes_vs_k_DA_B0", "beta_RelRes_vs_k_Sigma_B0"},
                  {"EM", "DA", "Sigma"},
                  "#beta relative resolution vs k (B0)",
                  "figs/resolutions/binned/beta_relres_vs_k_b0.png");
    PlotRelResVsK(inputFile,
                  {"beta_RelRes_vs_k_EM_RP", "beta_RelRes_vs_k_DA_RP", "beta_RelRes_vs_k_Sigma_RP"},
                  {"EM", "DA", "Sigma"},
                  "#beta relative resolution vs k (RP)",
                  "figs/resolutions/binned/beta_relres_vs_k_rp.png");

    // =================================================================
    // Legacy 2D maps (circle style), EPz correlation, and cross sections
    // =================================================================
    auto createCirclePlot = [&](const char* histName, const char* saveName,
                                const char* xTitle, const char* yTitle,
                                bool logX = true, bool logY = true) {
        TProfile2D* prof = (TProfile2D*)inputFile->Get(histName);
        if (!prof) return;

        TCanvas* c = new TCanvas(Form("c_circle_%s", histName), "Circle Plot", 1200, 900);
        c->SetLogx(logX);
        c->SetLogy(logY);
        c->SetGrid();
        c->SetRightMargin(0.15);

        TH2D* frame = new TH2D(Form("frame_%s", histName), Form(";%s;%s", xTitle, yTitle),
                               prof->GetNbinsX(), prof->GetXaxis()->GetXmin(), prof->GetXaxis()->GetXmax(),
                               prof->GetNbinsY(), prof->GetYaxis()->GetXmin(), prof->GetYaxis()->GetXmax());
        frame->Draw();

        for (int ix = 1; ix <= prof->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= prof->GetNbinsY(); iy++) {
                double entries = prof->GetBinEntries(prof->GetBin(ix, iy));
                if (entries < 10) continue;
                double x = prof->GetXaxis()->GetBinCenter(ix);
                double y = prof->GetYaxis()->GetBinCenter(iy);
                double rms = prof->GetBinError(ix, iy);
                double markerSize = TMath::Min(rms * 1000.0, 3.0);
                TMarker* marker = new TMarker(x, y, 20);
                marker->SetMarkerSize(markerSize);
                if (entries < 100) {
                    marker->SetMarkerStyle(24);
                    marker->SetMarkerColor(kBlue);
                } else {
                    marker->SetMarkerStyle(20);
                    marker->SetMarkerColor(kRed);
                }
                marker->Draw();
            }
        }
        SaveCanvas(c, saveName);
        delete frame;
        delete c;
    };

    auto createBestMethodPlot = [&](const char* histNameEM, const char* histNameDA, const char* histNameSigma,
                                    const char* saveName, const char* xTitle, const char* yTitle,
                                    bool logX = true, bool logY = true) {
        TProfile2D* profEM = (TProfile2D*)inputFile->Get(histNameEM);
        TProfile2D* profDA = (TProfile2D*)inputFile->Get(histNameDA);
        TProfile2D* profSigma = (TProfile2D*)inputFile->Get(histNameSigma);
        if (!profEM || !profDA || !profSigma) return;

        TCanvas* c = new TCanvas(Form("c_best_%s", histNameEM), "Best Method", 1200, 900);
        c->SetLogx(logX);
        c->SetLogy(logY);
        c->SetGrid();
        c->SetRightMargin(0.15);

        TH2D* frame = new TH2D(Form("frame_best_%s", histNameEM), Form(";%s;%s", xTitle, yTitle),
                               profEM->GetNbinsX(), profEM->GetXaxis()->GetXmin(), profEM->GetXaxis()->GetXmax(),
                               profEM->GetNbinsY(), profEM->GetYaxis()->GetXmin(), profEM->GetYaxis()->GetXmax());
        frame->Draw();

        for (int ix = 1; ix <= profEM->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= profEM->GetNbinsY(); iy++) {
                const double rmsEM = profEM->GetBinError(ix, iy);
                const double rmsDA = profDA->GetBinError(ix, iy);
                const double rmsSigma = profSigma->GetBinError(ix, iy);
                const int entriesEM = static_cast<int>(profEM->GetBinEntries(profEM->GetBin(ix, iy)));
                const int entriesDA = static_cast<int>(profDA->GetBinEntries(profDA->GetBin(ix, iy)));
                const int entriesSigma = static_cast<int>(profSigma->GetBinEntries(profSigma->GetBin(ix, iy)));

                int bestMethod = -1;
                double bestRms = std::numeric_limits<double>::infinity();
                int bestEntries = 0;
                if (std::isfinite(rmsEM) && rmsEM > 0.0 && entriesEM >= 5 && rmsEM < bestRms) {
                    bestRms = rmsEM;
                    bestMethod = 0;
                    bestEntries = entriesEM;
                }
                if (std::isfinite(rmsDA) && rmsDA > 0.0 && entriesDA >= 5 && rmsDA < bestRms) {
                    bestRms = rmsDA;
                    bestMethod = 1;
                    bestEntries = entriesDA;
                }
                if (std::isfinite(rmsSigma) && rmsSigma > 0.0 && entriesSigma >= 5 && rmsSigma < bestRms) {
                    bestRms = rmsSigma;
                    bestMethod = 2;
                    bestEntries = entriesSigma;
                }
                if (bestMethod < 0) continue;

                const int markerStyle = (bestEntries < 100) ? 24 : 20; // hollow for low stats, filled for high stats
                TMarker* marker = new TMarker(profEM->GetXaxis()->GetBinCenter(ix),
                                              profEM->GetYaxis()->GetBinCenter(iy), markerStyle);
                double markerSize = 0.3 + 25.0 * bestRms;
                if (markerSize > 4.0) markerSize = 4.0;
                marker->SetMarkerSize(markerSize);
                if (bestMethod == 0) marker->SetMarkerColor(kRed);
                else if (bestMethod == 1) marker->SetMarkerColor(kBlue);
                else marker->SetMarkerColor(kGreen + 2);
                marker->Draw();
            }
        }

        TLegend* leg = new TLegend(0.7, 0.82, 0.9, 0.95);
        TMarker mEM(0, 0, 20); mEM.SetMarkerColor(kRed);
        TMarker mDA(0, 0, 20); mDA.SetMarkerColor(kBlue);
        TMarker mSigma(0, 0, 20); mSigma.SetMarkerColor(kGreen + 2);
        leg->AddEntry(&mEM, "EM Best", "p");
        leg->AddEntry(&mDA, "DA Best", "p");
        leg->AddEntry(&mSigma, "Sigma Best", "p");
        leg->AddEntry((TObject*)nullptr, "Size = 0.3 + 25#timesRMS", "");
        leg->AddEntry((TObject*)nullptr, "Hollow: N<100, Filled: N#geq100", "");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        SaveCanvas(c, saveName);
        delete leg;
        delete frame;
        delete c;
    };

    createCirclePlot("Q2_RelRes_vs_xy_EM", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_EM.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_DA", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_DA.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_Sigma", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_Sigma.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_EM", "figs/resolutions/2d_maps/x_RelRes_xQ2_EM.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_DA", "figs/resolutions/2d_maps/x_RelRes_xQ2_DA.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_Sigma", "figs/resolutions/2d_maps/x_RelRes_xQ2_Sigma.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_EM", "figs/resolutions/2d_maps/y_RelRes_xQ2_EM.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_DA", "figs/resolutions/2d_maps/y_RelRes_xQ2_DA.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_Sigma", "figs/resolutions/2d_maps/y_RelRes_xQ2_Sigma.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("t_RelRes_vs_xpomQ2_B0", "figs/resolutions/2d_maps/t_RelRes_xpomQ2_B0.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("t_RelRes_vs_xpomQ2_RP", "figs/resolutions/2d_maps/t_RelRes_xpomQ2_RP.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("xpom_RelRes_vs_xpomQ2_B0", "figs/resolutions/2d_maps/xpom_RelRes_xpomQ2_B0.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("xpom_RelRes_vs_xpomQ2_RP", "figs/resolutions/2d_maps/xpom_RelRes_xpomQ2_RP.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("beta_RelRes_vs_betaQ2_B0", "figs/resolutions/2d_maps/beta_RelRes_betaQ2_B0.png", "#beta", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("beta_RelRes_vs_betaQ2_RP", "figs/resolutions/2d_maps/beta_RelRes_betaQ2_RP.png", "#beta", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("xL_RelRes_vs_xLQ2_B0", "figs/resolutions/2d_maps/xL_RelRes_xLQ2_B0.png", "x_{L}", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("xL_RelRes_vs_xLQ2_RP", "figs/resolutions/2d_maps/xL_RelRes_xLQ2_RP.png", "x_{L}", "Q^{2} [GeV^{2}]", false, true);

    createBestMethodPlot("Q2_RelRes_vs_xy_EM", "Q2_RelRes_vs_xy_DA", "Q2_RelRes_vs_xy_Sigma",
                         "figs/resolutions/2d_maps/Q2_RelRes_Q2x_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("x_RelRes_vs_xQ2_EM", "x_RelRes_vs_xQ2_DA", "x_RelRes_vs_xQ2_Sigma",
                         "figs/resolutions/2d_maps/x_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("y_RelRes_vs_xQ2_EM", "y_RelRes_vs_xQ2_DA", "y_RelRes_vs_xQ2_Sigma",
                         "figs/resolutions/2d_maps/y_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    if (auto* profMX2 = (TProfile2D*)inputFile->Get("MX2_RelRes_vs_MX2Q2")) {
        TCanvas* cMX2 = new TCanvas("c_MX2_relres_map", "MX2 Resolution", 1200, 900);
        cMX2->SetRightMargin(0.15);
        cMX2->SetLogx();
        cMX2->SetLogy();
        profMX2->SetTitle("M_{X}^{2} Resolution vs (M_{X}^{2}, Q^{2});M_{X}^{2} [GeV^{2}];Q^{2} [GeV^{2}]");
        profMX2->Draw("COLZ TEXT");
        SaveCanvas(cMX2, "figs/resolutions/2d_maps/MX2_RelRes_MX2Q2.png");
        delete cMX2;
    }

    if (auto* hEPz2D = (TH2D*)inputFile->Get("h_EPz_2D")) {
        TCanvas* cEPz = new TCanvas("c_EPz_2D", "E-pz Correlation", 1200, 900);
        cEPz->SetRightMargin(0.15);
        cEPz->SetGrid();
        hEPz2D->Draw("COLZ");
        SaveCanvas(cEPz, "figs/distributions/EPz_2D.png");
        delete cEPz;
    }

    TH1D* hDSigB0 = (TH1D*)inputFile->Get("dsigma_dt_B0");
    TH1D* hDSigRP = (TH1D*)inputFile->Get("dsigma_dt_RP");
    TH1D* hDSigMC = (TH1D*)inputFile->Get("dsigma_dt_MC");
    TH1D* hDSigSum = (TH1D*)inputFile->Get("dsigma_dt_Sum");
    if (hDSigSum) {
        hDSigSum = (TH1D*)hDSigSum->Clone("dsigma_dt_Sum_local");
        hDSigSum->SetDirectory(nullptr);
    } else if (hDSigB0 && hDSigRP) {
        hDSigSum = (TH1D*)hDSigB0->Clone("dsigma_dt_Sum_local");
        hDSigSum->SetDirectory(nullptr);
        hDSigSum->Add(hDSigRP);
    }
    TF1* fitExp = nullptr;
    double bValue = 0.0;
    double bError = 0.0;
    if (hDSigMC) {
        fitExp = new TF1("fit_exp_plot_final", "[0]*TMath::Exp(-[1]*x)", 0.01, 2.0);
        fitExp->SetParameters(1000.0, 5.0);
        fitExp->SetLineColor(kMagenta + 2);
        fitExp->SetLineWidth(3);
        fitExp->SetLineStyle(2);
        hDSigMC->Fit(fitExp, "RSQ");
        bValue = fitExp->GetParameter(1);
        bError = fitExp->GetParError(1);
    }
    if (hDSigMC && hDSigB0 && hDSigRP && hDSigSum && fitExp) {
        TCanvas* cLin = new TCanvas("c_dsigma_with_fit_linear", "dSigma/dt", 1200, 900);
        cLin->SetLogx();
        cLin->SetGrid();
        hDSigMC->SetLineColor(kBlack);
        hDSigMC->SetLineWidth(2);
        hDSigMC->Draw("HIST");
        hDSigB0->SetMarkerStyle(20); hDSigB0->SetMarkerColor(kRed); hDSigB0->SetLineColor(kRed); hDSigB0->Draw("PE SAME");
        hDSigRP->SetMarkerStyle(20); hDSigRP->SetMarkerColor(kBlue); hDSigRP->SetLineColor(kBlue); hDSigRP->Draw("PE SAME");
        hDSigSum->SetMarkerStyle(20); hDSigSum->SetMarkerColor(kOrange + 7); hDSigSum->SetLineColor(kOrange + 7); hDSigSum->Draw("PE SAME");
        fitExp->Draw("SAME");
        TLegend leg1(0.56, 0.55, 0.87, 0.9);
        leg1.SetBorderSize(0); leg1.SetFillStyle(0);
        leg1.AddEntry(hDSigMC, "MC Truth", "l");
        leg1.AddEntry(hDSigB0, "B0 Reco", "pe");
        leg1.AddEntry(hDSigRP, "RP Reco", "pe");
        leg1.AddEntry(hDSigSum, "B0+RP Sum", "pe");
        leg1.AddEntry(fitExp, Form("Fit: e^{-b|t|}, b = %.2f #pm %.2f GeV^{-2}", bValue, bError), "l");
        leg1.Draw();
        SaveCanvas(cLin, "figs/cross_sections/dsigma_dt_with_fit.png");
        delete cLin;

        TCanvas* cLog = new TCanvas("c_dsigma_with_fit_logy", "dSigma/dt", 1200, 900);
        cLog->SetLogx();
        cLog->SetLogy();
        cLog->SetGrid();
        hDSigMC->Draw("HIST");
        hDSigB0->Draw("PE SAME");
        hDSigRP->Draw("PE SAME");
        hDSigSum->Draw("PE SAME");
        fitExp->Draw("SAME");
        TLegend leg2(0.56, 0.55, 0.87, 0.9);
        leg2.SetBorderSize(0); leg2.SetFillStyle(0);
        leg2.AddEntry(hDSigMC, "MC Truth", "l");
        leg2.AddEntry(hDSigB0, "B0 Reco", "pe");
        leg2.AddEntry(hDSigRP, "RP Reco", "pe");
        leg2.AddEntry(hDSigSum, "B0+RP Sum", "pe");
        leg2.AddEntry(fitExp, Form("Fit: e^{-b|t|}, b = %.2f #pm %.2f GeV^{-2}", bValue, bError), "l");
        leg2.Draw();
        SaveCanvas(cLog, "figs/cross_sections/dsigma_dt_logy_with_fit.png");
        delete cLog;
    }
    delete fitExp;
    delete hDSigSum;

    TH3D* hD3MC = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_MC");
    TH3D* hD3B0 = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_B0");
    TH3D* hD3RP = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_RP");
    TH3D* hD3Sum = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_Sum");
    TH3D* hD3SumLocal = nullptr;
    if (!hD3Sum && hD3B0 && hD3RP) {
        hD3SumLocal = (TH3D*)hD3B0->Clone("d3sigma_dQ2dbeta_dxpom_Sum_local");
        hD3SumLocal->SetDirectory(nullptr);
        hD3SumLocal->Add(hD3RP);
        hD3Sum = hD3SumLocal;
    }
    if (hD3MC && hD3B0 && hD3RP) {
        const int nQ2 = hD3MC->GetNbinsX();
        const int nBeta = hD3MC->GetNbinsY();
        const int nXpom = hD3MC->GetNbinsZ();

        TCanvas* cBeta = new TCanvas("c_d3sigma_beta", "d3sigma vs beta", 2400, 1600);
        cBeta->Divide(nQ2, nXpom);
        for (int iQ2 = 1; iQ2 <= nQ2; ++iQ2) {
            for (int iXpom = 1; iXpom <= nXpom; ++iXpom) {
                int pad = (iXpom - 1) * nQ2 + iQ2;
                cBeta->cd(pad);
                gPad->SetLogy();
                TH1D* pMC = hD3MC->ProjectionY(Form("proj_beta_MC_%d_%d_pf", iQ2, iXpom), iQ2, iQ2, iXpom, iXpom);
                TH1D* pB0 = hD3B0->ProjectionY(Form("proj_beta_B0_%d_%d_pf", iQ2, iXpom), iQ2, iQ2, iXpom, iXpom);
                TH1D* pRP = hD3RP->ProjectionY(Form("proj_beta_RP_%d_%d_pf", iQ2, iXpom), iQ2, iQ2, iXpom, iXpom);
                TH1D* pSum = hD3Sum ? hD3Sum->ProjectionY(Form("proj_beta_SUM_%d_%d_pf", iQ2, iXpom), iQ2, iQ2, iXpom, iXpom) : nullptr;
                pMC->SetLineColor(kBlack);
                pB0->SetMarkerColor(kRed); pB0->SetMarkerStyle(20);
                pRP->SetMarkerColor(kBlue); pRP->SetMarkerStyle(20);
                if (pSum) { pSum->SetMarkerColor(kGreen + 2); pSum->SetMarkerStyle(21); }
                pMC->SetTitle(Form("Q^{2} bin %d, x_{pom} bin %d;#beta;d^{3}#sigma", iQ2, iXpom));
                pMC->Draw("HIST");
                pB0->Draw("PE SAME");
                pRP->Draw("PE SAME");
                if (pSum) pSum->Draw("PE SAME");
                if (iQ2 == 1 && iXpom == 1) {
                    TLegend leg(0.55, 0.63, 0.88, 0.89);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.AddEntry(pMC, "MC Truth", "l");
                    leg.AddEntry(pB0, "B0 Reco", "pe");
                    leg.AddEntry(pRP, "RP Reco", "pe");
                    if (pSum) leg.AddEntry(pSum, "B0+RP Sum", "pe");
                    leg.Draw();
                }
            }
        }
        SaveCanvas(cBeta, "figs/cross_sections/d3sigma_vs_beta.png");
        delete cBeta;

        TCanvas* cXpom = new TCanvas("c_d3sigma_xpom", "d3sigma vs xpom", 2400, 1000);
        cXpom->Divide(nQ2, nBeta);
        for (int iQ2 = 1; iQ2 <= nQ2; ++iQ2) {
            for (int iBeta = 1; iBeta <= nBeta; ++iBeta) {
                int pad = (iBeta - 1) * nQ2 + iQ2;
                cXpom->cd(pad);
                gPad->SetLogx();
                gPad->SetLogy();
                TH1D* pMC = hD3MC->ProjectionZ(Form("proj_xpom_MC_%d_%d_pf", iQ2, iBeta), iQ2, iQ2, iBeta, iBeta);
                TH1D* pB0 = hD3B0->ProjectionZ(Form("proj_xpom_B0_%d_%d_pf", iQ2, iBeta), iQ2, iQ2, iBeta, iBeta);
                TH1D* pRP = hD3RP->ProjectionZ(Form("proj_xpom_RP_%d_%d_pf", iQ2, iBeta), iQ2, iQ2, iBeta, iBeta);
                TH1D* pSum = hD3Sum ? hD3Sum->ProjectionZ(Form("proj_xpom_SUM_%d_%d_pf", iQ2, iBeta), iQ2, iQ2, iBeta, iBeta) : nullptr;
                pMC->SetLineColor(kBlack);
                pB0->SetMarkerColor(kRed); pB0->SetMarkerStyle(20);
                pRP->SetMarkerColor(kBlue); pRP->SetMarkerStyle(20);
                if (pSum) { pSum->SetMarkerColor(kGreen + 2); pSum->SetMarkerStyle(21); }
                pMC->SetTitle(Form("Q^{2} bin %d, #beta bin %d;x_{pom};d^{3}#sigma", iQ2, iBeta));
                pMC->Draw("HIST");
                pB0->Draw("PE SAME");
                pRP->Draw("PE SAME");
                if (pSum) pSum->Draw("PE SAME");
                if (iQ2 == 1 && iBeta == 1) {
                    TLegend leg(0.55, 0.63, 0.88, 0.89);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.AddEntry(pMC, "MC Truth", "l");
                    leg.AddEntry(pB0, "B0 Reco", "pe");
                    leg.AddEntry(pRP, "RP Reco", "pe");
                    if (pSum) leg.AddEntry(pSum, "B0+RP Sum", "pe");
                    leg.Draw();
                }
            }
        }
        SaveCanvas(cXpom, "figs/cross_sections/d3sigma_vs_xpom.png");
        delete cXpom;
    }
    delete hD3SumLocal;

    // Legacy path compatibility: mirror the same plots to older output names.
    auto copyPlotAlias = [](const char* src, const char* dst) {
        if (!src || !dst) return;
        if (gSystem->AccessPathName(src)) return;
        std::string dstPath(dst);
        const size_t slash = dstPath.find_last_of('/');
        if (slash != std::string::npos) {
            gSystem->mkdir(dstPath.substr(0, slash).c_str(), true);
        }
        gSystem->CopyFile(src, dst, true);
    };

    const std::vector<std::pair<const char*, const char*>> legacyPlotAliases = {
        {"figs/inclusive/histos/q2_methods_hist.png", "figs/distributions/Q2_hist.png"},
        {"figs/inclusive/histos/q2_methods_pdf.png", "figs/distributions/Q2_pdf.png"},
        {"figs/inclusive/histos/xbj_methods_hist.png", "figs/distributions/x_hist.png"},
        {"figs/inclusive/histos/xbj_methods_pdf.png", "figs/distributions/x_pdf.png"},
        {"figs/inclusive/histos/y_methods_hist.png", "figs/distributions/y_hist.png"},
        {"figs/inclusive/histos/y_methods_pdf.png", "figs/distributions/y_pdf.png"},
        {"figs/diffractive/histos/t_distributions.png", "figs/distributions/t_distributions.png"},
        {"figs/diffractive/histos/t_distributions_logy.png", "figs/distributions/t_distributions_logy.png"},
        {"figs/diffractive/histos/theta_distributions.png", "figs/distributions/theta_distributions.png"},
        {"figs/inclusive/resolution/q2_relres_em.png", "figs/resolutions/simple/DDIS_Q2RelRes_EM.png"},
        {"figs/inclusive/resolution/q2_relres_da.png", "figs/resolutions/simple/DDIS_Q2RelRes_DA.png"},
        {"figs/inclusive/resolution/q2_relres_sigma.png", "figs/resolutions/simple/DDIS_Q2RelRes_Sigma.png"},
        {"figs/inclusive/resolution/xbj_relres_em.png", "figs/resolutions/simple/DDIS_RelRes_xBj_EM.png"},
        {"figs/inclusive/resolution/xbj_relres_da.png", "figs/resolutions/simple/DDIS_RelRes_xBj_DA.png"},
        {"figs/inclusive/resolution/xbj_relres_sigma.png", "figs/resolutions/simple/DDIS_RelRes_x_Sigma.png"},
        {"figs/inclusive/resolution/y_relres_em.png", "figs/resolutions/simple/DDIS_RelRes_y_EM.png"},
        {"figs/inclusive/resolution/y_relres_da.png", "figs/resolutions/simple/DDIS_RelRes_y_DA.png"},
        {"figs/inclusive/resolution/y_relres_sigma.png", "figs/resolutions/simple/DDIS_RelRes_y_Sigma.png"},
        {"figs/diffractive/resolution/t_res_b0.png", "figs/resolutions/simple/t_resolution_B0.png"},
        {"figs/diffractive/resolution/t_res_rp.png", "figs/resolutions/simple/t_resolution_RP.png"},
        {"figs/inclusive/resolution/q2_relres_binned_em.png", "figs/resolutions/binned/DDIS_Q2RelRes_binned_EM.png"},
        {"figs/inclusive/resolution/q2_relres_binned_da.png", "figs/resolutions/binned/DDIS_Q2RelRes_binned_DA.png"},
        {"figs/inclusive/resolution/q2_relres_binned_sigma.png", "figs/resolutions/binned/DDIS_Q2RelRes_binned_Sigma.png"},
        {"figs/inclusive/resolution/xbj_relres_binned_em.png", "figs/resolutions/binned/DDIS_RelRes_binned_x_EM.png"},
        {"figs/inclusive/resolution/xbj_relres_binned_da.png", "figs/resolutions/binned/DDIS_RelRes_binned_x_DA.png"},
        {"figs/inclusive/resolution/xbj_relres_binned_sigma.png", "figs/resolutions/binned/DDIS_RelRes_binned_x_Sigma.png"},
        {"figs/inclusive/resolution/y_relres_binned_em.png", "figs/resolutions/binned/DDIS_RelRes_binned_y_EM.png"},
        {"figs/inclusive/resolution/y_relres_binned_da.png", "figs/resolutions/binned/DDIS_RelRes_binned_y_DA.png"},
        {"figs/inclusive/resolution/y_relres_binned_sigma.png", "figs/resolutions/binned/DDIS_RelRes_binned_y_Sigma.png"},
        {"figs/diffractive/resolution/t_relres_binned_b0.png", "figs/resolutions/binned/t_resolution_binned_B0.png"},
        {"figs/diffractive/resolution/t_relres_binned_rp.png", "figs/resolutions/binned/t_resolution_binned_RP.png"},

        {"figs/distributions/EPz_distribution.png", "figs/EPz_distribution.png"},
        {"figs/distributions/EPz_distribution_logY.png", "figs/EPz_distribution_logY.png"},
        {"figs/distributions/MX2_distribution.png", "figs/MX2_distribution.png"},
        {"figs/distributions/MX2_distribution_logY.png", "figs/MX2_distribution_logY.png"},
        {"figs/diffractive/histos/MX2_comparison.png", "figs/MX2_comparison.png"},
        {"figs/resolutions/2d_maps/Q2_RelRes_Q2x_EM.png", "figs/Q2_RelRes_Q2x_EM.png"},
        {"figs/resolutions/2d_maps/Q2_RelRes_Q2x_DA.png", "figs/Q2_RelRes_Q2x_DA.png"},
        {"figs/resolutions/2d_maps/Q2_RelRes_Q2x_Sigma.png", "figs/Q2_RelRes_Q2x_Sigma.png"},
        {"figs/resolutions/2d_maps/Q2_RelRes_Q2x_BestMethod.png", "figs/Q2_RelRes_Q2x_BestMethod.png"},
        {"figs/distributions/Q2_hist.png", "figs/Q2_hist.png"},
        {"figs/distributions/beta_distributions_logy.png", "figs/beta_distributions_logy.png"},
        {"figs/distributions/beta_vs_Q2.png", "figs/beta_vs_Q2.png"},
        {"figs/distributions/beta_vs_t.png", "figs/beta_vs_t.png"},
        {"figs/distributions/beta_vs_xpom.png", "figs/beta_vs_xpom.png"},
        {"figs/resolutions/simple/beta_resolution_B0.png", "figs/beta_resolution_B0.png"},
        {"figs/resolutions/simple/beta_resolution_RP.png", "figs/beta_resolution_RP.png"},
        {"figs/cross_sections/d3sigma_vs_beta.png", "figs/d3sigma_vs_beta.png"},
        {"figs/cross_sections/d3sigma_vs_xpom.png", "figs/d3sigma_vs_xpom.png"},
        {"figs/diffractive/histos/dsigma_dt.png", "figs/dsigma_dt.png"},
        {"figs/cross_sections/dsigma_dt_logy_with_fit.png", "figs/dsigma_dt_logy.png"},
        {"figs/distributions/eta_max_distribution.png", "figs/eta_max_distribution.png"},
        {"figs/distributions/eta_max_distribution_logY.png", "figs/eta_max_distribution_logY.png"},
        {"figs/response_matrices/response_matrix_Q2_EM.png", "figs/response_matrix_Q2_EM.png"},
        {"figs/response_matrices/response_matrix_Q2_DA.png", "figs/response_matrix_Q2_DA.png"},
        {"figs/response_matrices/response_matrix_Q2_Esigma.png", "figs/response_matrix_Q2_Esigma.png"},
        {"figs/response_matrices/response_matrix_x_EM.png", "figs/response_matrix_x_EM.png"},
        {"figs/response_matrices/response_matrix_x_DA.png", "figs/response_matrix_x_DA.png"},
        {"figs/response_matrices/response_matrix_x_Sigma.png", "figs/response_matrix_x_Sigma.png"},
        {"figs/response_matrices/response_matrix_y_EM.png", "figs/response_matrix_y_EM.png"},
        {"figs/response_matrices/response_matrix_y_DA.png", "figs/response_matrix_y_DA.png"},
        {"figs/response_matrices/response_matrix_y_Sigma.png", "figs/response_matrix_y_Sigma.png"},
        {"figs/response_matrices/response_matrix_t_B0.png", "figs/response_matrix_t_B0.png"},
        {"figs/response_matrices/response_matrix_t_RP.png", "figs/response_matrix_t_RP.png"},
        {"figs/response_matrices/response_matrix_xL_B0.png", "figs/response_matrix_xL_B0.png"},
        {"figs/response_matrices/response_matrix_xL_RP.png", "figs/response_matrix_xL_RP.png"},
        {"figs/response_matrices/response_matrix_xpom_B0.png", "figs/response_matrix_xpom_B0.png"},
        {"figs/response_matrices/response_matrix_xpom_RP.png", "figs/response_matrix_xpom_RP.png"},
        {"figs/response_matrices/response_matrix_beta_B0.png", "figs/response_matrix_beta_B0.png"},
        {"figs/response_matrices/response_matrix_beta_RP.png", "figs/response_matrix_beta_RP.png"},
        {"figs/response_matrices/response_matrix_MX2.png", "figs/response_matrix_MX2.png"},
        {"figs/response_matrices/response_matrix_t_B0.png", "figs/response_matrix_B0_cutFirstBin.png"},
        {"figs/response_matrices/response_matrix_xL_B0.png", "figs/response_matrix_xL_B0_cutFirstBin.png"},
        {"figs/response_matrices/response_matrix_xpom_B0.png", "figs/response_matrix_xpom_B0_cutFirstBin.png"},
        {"figs/diffractive/histos/t_distributions.png", "figs/t_distributions.png"},
        {"figs/diffractive/histos/t_distributions_logy.png", "figs/t_distributions_logy.png"},
        {"figs/distributions/t_pdf_comparison.png", "figs/t_pdf_comparison.png"},
        {"figs/diffractive/resolution/t_res_b0.png", "figs/t_resolution_B0.png"},
        {"figs/diffractive/resolution/t_res_rp.png", "figs/t_resolution_RP.png"},
        {"figs/distributions/theta_comparison_B0_acceptance.png", "figs/theta_comparison_B0_acceptance.png"},
        {"figs/distributions/theta_comparison_B0_acceptance_logxy.png", "figs/theta_comparison_B0_acceptance_logxy.png"},
        {"figs/diffractive/histos/theta_distributions.png", "figs/theta_distributions.png"},
        {"figs/resolutions/simple/xL_resolution_B0.png", "figs/xL_resolution_B0.png"},
        {"figs/resolutions/simple/xL_resolution_RP.png", "figs/xL_resolution_RP.png"},
        {"figs/resolutions/simple/MX2_resolution.png", "figs/MX2_resolution.png"},
        {"figs/diffractive/histos/xL_distributions.png", "figs/x_L_comparison.png"},
        {"figs/diffractive/histos/xL_distributions.png", "figs/x_L_comparison_logy.png"},
        {"figs/distributions/xpom_2D_comparison_MC.png", "figs/xpom_2D_comparison_MC.png"},
        {"figs/distributions/xpom_2D_comparison_B0.png", "figs/xpom_2D_comparison_B0.png"},
        {"figs/distributions/xpom_2D_comparison_RP.png", "figs/xpom_2D_comparison_RP.png"},
        {"figs/distributions/xpom_comparison_B0_logxy.png", "figs/xpom_B0_comparison_firstBinCut.png"},
        {"figs/distributions/xpom_comparison_B0_logxy.png", "figs/xpom_comparison_B0_logxy.png"},
        {"figs/distributions/xpom_comparison_MC_logxy.png", "figs/xpom_comparison_MC_logxy.png"},
        {"figs/distributions/xpom_comparison_RP_logxy.png", "figs/xpom_comparison_RP_logxy.png"},
        {"figs/distributions/xpom_comparison_all_logxy.png", "figs/xpom_comparison_all_logxy.png"},
        {"figs/distributions/xpom_comparison_all_logxy.png", "figs/xpom_comparison_logx.png"},
        {"figs/distributions/xpom_comparison_all_logxy.png", "figs/xpom_comparison_logxy.png"},
        {"figs/resolutions/simple/xpom_resolution_B0.png", "figs/xpom_resolution_B0.png"},
        {"figs/resolutions/simple/xpom_resolution_RP.png", "figs/xpom_resolution_RP.png"},
        {"figs/diffractive/histos/t_distributions.png", "figs/t_B0_comparison_firstBinCut.png"},
        {"figs/diffractive/histos/t_distributions.png", "figs/t_correlation_combined.png"},
        {"figs/diffractive/histos/xL_distributions.png", "figs/xL_B0_comparison_firstBinCut.png"},
        {"figs/distributions/Q2_hist.png", "figs/Q2_comparison.png"}
    };
    // Keep legacy duplicate plot paths opt-in to reduce clutter in figs/.
    const bool kEnableLegacyPlotAliases = false;
    if (kEnableLegacyPlotAliases) {
        for (const auto& alias : legacyPlotAliases) {
            copyPlotAlias(alias.first, alias.second);
        }
    }

    inputFile->Close();
    delete inputFile;

    return 0;
}
