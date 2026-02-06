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
#include <TMath.h>
#include <TString.h>
#include <TLine.h>
#include <TColor.h>

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

static void PlotXQ2Density(TFile* inputFile) {
    if (!inputFile) return;
    TTree* tree = (TTree*)inputFile->Get("Q2_tree");
    if (!tree) {
        std::cerr << "Warning: Q2_tree not found; skipping x-Q2 density plot." << std::endl;
        return;
    }

    const bool hasRecoBranches = tree->GetBranch("x_EM") && tree->GetBranch("Q2_EM");
    const bool hasTruthBranches = tree->GetBranch("x_truth") && tree->GetBranch("Q2_truth");
    if (!hasRecoBranches && !hasTruthBranches) {
        TH2D* h_existing = (TH2D*)inputFile->Get("xQ2_reco");
        if (!h_existing) {
            h_existing = (TH2D*)inputFile->Get("xQ2_truth");
        }
        if (!h_existing) {
            std::cerr << "Warning: Missing branches x_EM/Q2_EM and x_truth/Q2_truth and no xQ2_reco/xQ2_truth histogram found; "
                         "skipping x-Q2 density plot." << std::endl;
            return;
        }

        TCanvas* c = new TCanvas("c_xq2_density", "x-Q2 Density", 1200, 1000);
        gStyle->SetOptStat(0);
        c->SetRightMargin(0.15);
        c->SetLeftMargin(0.15);
        c->SetTopMargin(0.1);
        c->SetBottomMargin(0.1);
        c->SetLogx();
        c->SetLogy();

        h_existing->SetTitle("");
        h_existing->GetXaxis()->SetTitle("x_{Bj}");
        h_existing->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_existing->GetXaxis()->SetTitleOffset(1.3);
        h_existing->GetYaxis()->SetTitleOffset(1.3);
        h_existing->Draw("COLZ");

        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetNDC();
        latex.SetTextColor(kBlack);
        const std::string simLabel = BuildSimLabel(inputFile);
        latex.DrawLatex(0.2, 0.92, simLabel.c_str());
        latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

        c->Update();
        SaveCanvas(c, "figs/inclusive/histos/xbj_q2_density.png");
        delete c;

        std::cerr << "Note: Q2_tree lacks reco/truth branches; using stored xQ2 histogram from skim output." << std::endl;
        return;
    }

    float x_val = -999.0f;
    float Q2_val = -999.0f;
    if (hasRecoBranches) {
        tree->SetBranchAddress("x_EM", &x_val);
        tree->SetBranchAddress("Q2_EM", &Q2_val);
    } else {
        tree->SetBranchAddress("x_truth", &x_val);
        tree->SetBranchAddress("Q2_truth", &Q2_val);
    }

    const int n_x_bins = 120;
    const int n_q2_bins = 120;
    const double q2_min = 1.0;
    const double q2_max = 1000.0;
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, n_x_bins);
    std::vector<Double_t> q2_bins = GetLogBins(q2_min, q2_max, n_q2_bins);

    TH2D* h = new TH2D("xQ2_density_tmp",
                       "Event Density;x_{Bj};Q^{2} [GeV^{2}]",
                       x_bins.size() - 1, x_bins.data(),
                       q2_bins.size() - 1, q2_bins.data());

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (std::isfinite(x_val) && x_val > 0.0f && x_val < 1.0f &&
            std::isfinite(Q2_val) && Q2_val > 0.0f) {
            h->Fill(x_val, Q2_val);
        }
    }

    TCanvas* c = new TCanvas("c_xq2_density", "x-Q2 Density", 1200, 1000);
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

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.2, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

    c->Update();
    SaveCanvas(c, "figs/inclusive/histos/xbj_q2_density.png");

    delete h;
    delete c;
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

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.2, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
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
static void DrawBinningGridWithCounts(TH2D* hist,
                                      const std::vector<double>& xbins,
                                      const std::vector<double>& ybins,
                                      bool logX,
                                      bool logY) {
    if (!hist || xbins.size() < 2 || ybins.size() < 2) return;

    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double ymin = hist->GetYaxis()->GetXmin();
    const double ymax = hist->GetYaxis()->GetXmax();

    const int lineColor = kRed;
    const int lineWidth = 2;
    const int lineStyle = 1;

    for (double xval : xbins) {
        if (xval <= 0.0) continue;
        if (xval < xmin || xval > xmax) continue;
        TLine* line = new TLine(xval, ymin, xval, ymax);
        line->SetLineColor(lineColor);
        line->SetLineWidth(lineWidth);
        line->SetLineStyle(lineStyle);
        line->Draw("same");
    }

    for (double yval : ybins) {
        if (yval <= 0.0) continue;
        if (yval < ymin || yval > ymax) continue;
        TLine* line = new TLine(xmin, yval, xmax, yval);
        line->SetLineColor(lineColor);
        line->SetLineWidth(lineWidth);
        line->SetLineStyle(lineStyle);
        line->Draw("same");
    }

    TLatex latex;
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.022);
    latex.SetTextAlign(22);
    latex.SetTextFont(42);

    const double tiny = 1e-12;
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

            const double count = hist->Integral(binxlow, binxhigh, binylow, binyhigh);
            if (count > 0.0) {
                latex.DrawLatex(xcenter, ycenter, Form("%d", (int)std::lround(count)));
            }
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
                          bool overlayLogY = false) {
    if (!h3 || edges.size() < 2) return;
    const int nSlices = static_cast<int>(edges.size()) - 1;
    const int nCols = 2;
    const int nRows = (nSlices + nCols - 1) / nCols;

    TCanvas* c = new TCanvas(Form("c_phase_slices_%d", sliceAxis), saveName, 1200, 1000);
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

        if (overlayXBins && overlayYBins) {
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

static void PlotPhaseSpaceSlices(TFile* inputFile) {
    if (!inputFile) return;
    TH3D* h3 = (TH3D*)inputFile->Get("phase3D_reco");
    if (!h3) {
        std::cerr << "Warning: phase3D_reco not found; skipping phase-space slices." << std::endl;
        return;
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

    DrawSliceGrid(h3, 1, q2_edges,
                  "yz",
                  "x_{pom}", "#beta", "Q^{2} [GeV^{2}]",
                  "figs/diffractive/histos/phase_slices_q2.png",
                  true, false);
    DrawSliceGrid(h3, 2, xpom_edges,
                  "xz",
                  "#beta", "Q^{2} [GeV^{2}]", "x_{pom}",
                  "figs/diffractive/histos/phase_slices_xpom.png",
                  false, true);
    DrawSliceGrid(h3, 3, beta_edges,
                  "xy",
                  "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                  "figs/diffractive/histos/phase_slices_beta.png",
                  true, true,
                  &xpom_overlay_bins, &q2_overlay_bins,
                  true, true);
}

static void PlotEfficiencyCorrectedComparison(TFile* inputFile,
                                              const char* effTruthHistName,
                                              const char* recoHistName,
                                              const char* pseudoHistName,
                                              const char* plotTruthHistName,
                                              const char* xLabel,
                                              const char* title,
                                              const char* saveName,
                                              const bool logX,
                                              const char* recoMethodLabel) {
    if (!inputFile) return;
    TH1* h_truth = (TH1*)inputFile->Get(effTruthHistName);
    TH1* h_reco  = (TH1*)inputFile->Get(recoHistName);
    TH1* h_pdata = (TH1*)inputFile->Get(pseudoHistName);
    if (!h_truth || !h_reco || !h_pdata) {
        std::cerr << "Warning: Missing hist(s) for efficiency correction: "
                  << effTruthHistName << ", " << recoHistName << ", " << pseudoHistName << std::endl;
        return;
    }

    TH1D* h_eff = (TH1D*)h_reco->Clone(Form("eff_%s", recoHistName));
    h_eff->SetDirectory(nullptr);
    h_eff->Divide(h_reco, h_truth, 1.0, 1.0, "B");

    TH1D* h_truth_plot = nullptr;
    if (plotTruthHistName && plotTruthHistName[0] != '\0') {
        TH1* h_truth_overlay = (TH1*)inputFile->Get(plotTruthHistName);
        if (!h_truth_overlay) {
            std::cerr << "Warning: Truth overlay histogram " << plotTruthHistName
                      << " not found; using efficiency truth." << std::endl;
        } else {
            h_truth_plot = (TH1D*)h_truth_overlay->Clone(Form("truth_overlay_%s", plotTruthHistName));
            h_truth_plot->SetDirectory(nullptr);
        }
    }
    if (!h_truth_plot) {
        h_truth_plot = (TH1D*)h_truth->Clone(Form("truth_overlay_%s", effTruthHistName));
        h_truth_plot->SetDirectory(nullptr);
    }

    TH1D* h_corr = (TH1D*)h_pdata->Clone(Form("corr_%s", pseudoHistName));
    h_corr->SetDirectory(nullptr);
    const int nbins = h_corr->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        const double eff = h_eff->GetBinContent(i);
        const double pdata = h_pdata->GetBinContent(i);
        const double pdata_err = h_pdata->GetBinError(i);
        if (eff > 0.0) {
            h_corr->SetBinContent(i, pdata / eff);
            h_corr->SetBinError(i, pdata_err / eff);
        } else {
            h_corr->SetBinContent(i, 0.0);
            h_corr->SetBinError(i, 0.0);
        }
    }

    TCanvas* c = new TCanvas(Form("c_eff_%s", pseudoHistName), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();

    h_pdata->SetStats(false);
    h_pdata->SetTitle(title);
    h_pdata->GetXaxis()->SetTitle(xLabel);
    h_pdata->GetYaxis()->SetTitle("Counts");
    h_pdata->GetXaxis()->SetTitleOffset(1.1);
    h_pdata->GetYaxis()->SetTitleOffset(1.2);

    h_pdata->SetMarkerStyle(24);
    h_pdata->SetMarkerSize(1.1);
    h_pdata->SetMarkerColor(kBlack);
    h_pdata->SetLineColor(kBlack);

    h_corr->SetMarkerStyle(20);
    h_corr->SetMarkerSize(1.1);
    h_corr->SetMarkerColor(kBlack);
    h_corr->SetLineColor(kBlack);

    double max_val = 0.0;
    for (int i = 1; i <= h_pdata->GetNbinsX(); ++i) {
        max_val = std::max(max_val, h_pdata->GetBinContent(i) + h_pdata->GetBinError(i));
        max_val = std::max(max_val, h_corr->GetBinContent(i) + h_corr->GetBinError(i));
        max_val = std::max(max_val, h_truth->GetBinContent(i));
    }
    if (max_val > 0.0) {
        h_pdata->SetMaximum(max_val * 1.25);
        h_pdata->SetMinimum(0.0);
    }

    h_truth_plot->SetLineColor(kRed+1);
    h_truth_plot->SetLineWidth(2);
    h_truth_plot->SetMarkerStyle(0);

    h_pdata->Draw("PE");
    h_truth_plot->Draw("HIST SAME");
    h_corr->Draw("PE SAME");

    TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h_truth_plot, "Truth (Pseudo-data)", "l");
    legend->AddEntry(h_pdata, Form("Reco (%s) uncorrected", recoMethodLabel), "p");
    legend->AddEntry(h_corr, Form("Reco (%s) corrected", recoMethodLabel), "p");
    legend->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.2, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
    delete h_eff;
    delete h_corr;
    delete h_truth_plot;
}

static void PlotEfficiencyCorrectedComparisonBR(TFile* inputFile,
                                                const char* effTruthHistNameB0,
                                                const char* recoHistNameB0,
                                                const char* pseudoHistNameB0,
                                                const char* effTruthHistNameRP,
                                                const char* recoHistNameRP,
                                                const char* pseudoHistNameRP,
                                                const char* plotTruthHistNameB0,
                                                const char* plotTruthHistNameRP,
                                                const char* truthOverlayHistName,
                                                const char* xLabel,
                                                const char* title,
                                                const char* saveName,
                                                const bool logX,
                                                const bool logY,
                                                const char* recoMethodLabel,
                                                const bool includeSum) {
    if (!inputFile) return;
    TH1* h_truth_b0 = (TH1*)inputFile->Get(effTruthHistNameB0);
    TH1* h_reco_b0 = (TH1*)inputFile->Get(recoHistNameB0);
    TH1* h_pdata_b0 = (TH1*)inputFile->Get(pseudoHistNameB0);
    TH1* h_truth_rp = (TH1*)inputFile->Get(effTruthHistNameRP);
    TH1* h_reco_rp = (TH1*)inputFile->Get(recoHistNameRP);
    TH1* h_pdata_rp = (TH1*)inputFile->Get(pseudoHistNameRP);

    if (!h_truth_b0 || !h_reco_b0 || !h_pdata_b0 ||
        !h_truth_rp || !h_reco_rp || !h_pdata_rp) {
        std::cerr << "Warning: Missing hist(s) for B0/RP efficiency correction: "
                  << effTruthHistNameB0 << ", "
                  << recoHistNameB0 << ", " << pseudoHistNameB0 << ", "
                  << effTruthHistNameRP << ", " << recoHistNameRP << ", "
                  << pseudoHistNameRP << std::endl;
        return;
    }

    TH1D* h_truth_plot = nullptr;
    if (truthOverlayHistName && truthOverlayHistName[0] != '\0') {
        TH1* h_truth_overlay = (TH1*)inputFile->Get(truthOverlayHistName);
        if (!h_truth_overlay) {
            std::cerr << "Warning: Truth overlay histogram " << truthOverlayHistName
                      << " not found; falling back to pseudo-data truth sum." << std::endl;
        } else {
            h_truth_plot = (TH1D*)h_truth_overlay->Clone(Form("truth_overlay_%s", truthOverlayHistName));
            h_truth_plot->SetDirectory(nullptr);
        }
    }
    if (!h_truth_plot) {
        TH1* h_plot_truth_b0 = (TH1*)inputFile->Get(plotTruthHistNameB0);
        TH1* h_plot_truth_rp = (TH1*)inputFile->Get(plotTruthHistNameRP);
        if (!h_plot_truth_b0 || !h_plot_truth_rp) {
            std::cerr << "Warning: Missing plot-truth hist(s) for B0/RP overlay: "
                      << plotTruthHistNameB0 << ", " << plotTruthHistNameRP
                      << "; using efficiency truth sum." << std::endl;
            h_truth_plot = (TH1D*)h_truth_b0->Clone(Form("truth_sum_%s", effTruthHistNameB0));
            h_truth_plot->SetDirectory(nullptr);
            h_truth_plot->Add(h_truth_rp);
        } else {
            h_truth_plot = (TH1D*)h_plot_truth_b0->Clone(Form("truth_sum_%s", plotTruthHistNameB0));
            h_truth_plot->SetDirectory(nullptr);
            h_truth_plot->Add(h_plot_truth_rp);
        }
    }

    TH1D* h_eff_b0 = (TH1D*)h_reco_b0->Clone(Form("eff_%s", recoHistNameB0));
    h_eff_b0->SetDirectory(nullptr);
    h_eff_b0->Divide(h_reco_b0, h_truth_b0, 1.0, 1.0, "B");

    TH1D* h_eff_rp = (TH1D*)h_reco_rp->Clone(Form("eff_%s", recoHistNameRP));
    h_eff_rp->SetDirectory(nullptr);
    h_eff_rp->Divide(h_reco_rp, h_truth_rp, 1.0, 1.0, "B");

    TH1D* h_corr_b0 = (TH1D*)h_pdata_b0->Clone(Form("corr_%s", pseudoHistNameB0));
    h_corr_b0->SetDirectory(nullptr);
    TH1D* h_corr_rp = (TH1D*)h_pdata_rp->Clone(Form("corr_%s", pseudoHistNameRP));
    h_corr_rp->SetDirectory(nullptr);

    const int nbins = h_corr_b0->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        const double eff_b0 = h_eff_b0->GetBinContent(i);
        const double eff_rp = h_eff_rp->GetBinContent(i);
        const double pdata_b0 = h_pdata_b0->GetBinContent(i);
        const double pdata_rp = h_pdata_rp->GetBinContent(i);
        const double pdata_b0_err = h_pdata_b0->GetBinError(i);
        const double pdata_rp_err = h_pdata_rp->GetBinError(i);

        if (eff_b0 > 0.0) {
            h_corr_b0->SetBinContent(i, pdata_b0 / eff_b0);
            h_corr_b0->SetBinError(i, pdata_b0_err / eff_b0);
        } else {
            h_corr_b0->SetBinContent(i, 0.0);
            h_corr_b0->SetBinError(i, 0.0);
        }

        if (eff_rp > 0.0) {
            h_corr_rp->SetBinContent(i, pdata_rp / eff_rp);
            h_corr_rp->SetBinError(i, pdata_rp_err / eff_rp);
        } else {
            h_corr_rp->SetBinContent(i, 0.0);
            h_corr_rp->SetBinError(i, 0.0);
        }
    }

    TH1D* h_sum_uncorr = nullptr;
    TH1D* h_sum_corr = nullptr;
    if (includeSum) {
        h_sum_uncorr = (TH1D*)h_pdata_b0->Clone(Form("sum_uncorr_%s", pseudoHistNameB0));
        h_sum_uncorr->SetDirectory(nullptr);
        h_sum_uncorr->Add(h_pdata_rp);
        h_sum_corr = (TH1D*)h_corr_b0->Clone(Form("sum_corr_%s", pseudoHistNameB0));
        h_sum_corr->SetDirectory(nullptr);
        h_sum_corr->Add(h_corr_rp);
    }

    TCanvas* c = new TCanvas(Form("c_eff_%s", recoHistNameB0), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h_pdata_b0->SetStats(false);
    h_pdata_b0->SetTitle(title);
    h_pdata_b0->GetXaxis()->SetTitle(xLabel);
    h_pdata_b0->GetYaxis()->SetTitle("Counts");
    h_pdata_b0->GetXaxis()->SetTitleOffset(1.1);
    h_pdata_b0->GetYaxis()->SetTitleOffset(1.2);

    const int color_b0 = kRed + 1;
    const int color_rp = kBlue + 1;

    h_pdata_b0->SetMarkerStyle(24);
    h_pdata_b0->SetMarkerSize(1.0);
    h_pdata_b0->SetMarkerColor(color_b0);
    h_pdata_b0->SetLineColor(color_b0);

    h_corr_b0->SetMarkerStyle(20);
    h_corr_b0->SetMarkerSize(1.0);
    h_corr_b0->SetMarkerColor(color_b0);
    h_corr_b0->SetLineColor(color_b0);

    h_pdata_rp->SetMarkerStyle(24);
    h_pdata_rp->SetMarkerSize(1.0);
    h_pdata_rp->SetMarkerColor(color_rp);
    h_pdata_rp->SetLineColor(color_rp);

    h_corr_rp->SetMarkerStyle(20);
    h_corr_rp->SetMarkerSize(1.0);
    h_corr_rp->SetMarkerColor(color_rp);
    h_corr_rp->SetLineColor(color_rp);

    if (includeSum) {
        h_sum_uncorr->SetMarkerStyle(25);
        h_sum_uncorr->SetMarkerSize(1.0);
        h_sum_uncorr->SetMarkerColor(kGreen + 2);
        h_sum_uncorr->SetLineColor(kGreen + 2);

        h_sum_corr->SetMarkerStyle(21);
        h_sum_corr->SetMarkerSize(1.0);
        h_sum_corr->SetMarkerColor(kGreen + 2);
        h_sum_corr->SetLineColor(kGreen + 2);
    }

    h_truth_plot->SetLineColor(kBlack);
    h_truth_plot->SetLineWidth(2);
    h_truth_plot->SetMarkerStyle(0);

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= h_pdata_b0->GetNbinsX(); ++i) {
        const double vals[] = {
            h_pdata_b0->GetBinContent(i), h_corr_b0->GetBinContent(i),
            h_pdata_rp->GetBinContent(i), h_corr_rp->GetBinContent(i),
            h_truth_plot->GetBinContent(i)
        };
        max_val = std::max(max_val, h_pdata_b0->GetBinContent(i) + h_pdata_b0->GetBinError(i));
        max_val = std::max(max_val, h_corr_b0->GetBinContent(i) + h_corr_b0->GetBinError(i));
        max_val = std::max(max_val, h_pdata_rp->GetBinContent(i) + h_pdata_rp->GetBinError(i));
        max_val = std::max(max_val, h_corr_rp->GetBinContent(i) + h_corr_rp->GetBinError(i));
        if (includeSum) {
            max_val = std::max(max_val, h_sum_uncorr->GetBinContent(i) + h_sum_uncorr->GetBinError(i));
            max_val = std::max(max_val, h_sum_corr->GetBinContent(i) + h_sum_corr->GetBinError(i));
        }
        max_val = std::max(max_val, h_truth_plot->GetBinContent(i));
        for (double v : vals) {
            if (v > 0.0 && v < min_pos) min_pos = v;
        }
        if (includeSum) {
            const double sum_vals[] = {h_sum_uncorr->GetBinContent(i), h_sum_corr->GetBinContent(i)};
            for (double v : sum_vals) {
                if (v > 0.0 && v < min_pos) min_pos = v;
            }
        }
    }
    if (max_val > 0.0) {
        h_pdata_b0->SetMaximum(max_val * 1.25);
        if (logY) {
            const double min_y = (min_pos < std::numeric_limits<double>::max()) ? (min_pos * 0.5) : 1e-6;
            h_pdata_b0->SetMinimum(std::max(min_y, 1e-6));
        } else {
            h_pdata_b0->SetMinimum(0.0);
        }
    }

    h_pdata_b0->Draw("PE");
    h_truth_plot->Draw("HIST SAME");
    h_corr_b0->Draw("PE SAME");
    h_pdata_rp->Draw("PE SAME");
    h_corr_rp->Draw("PE SAME");
    if (includeSum) {
        h_sum_uncorr->Draw("PE SAME");
        h_sum_corr->Draw("PE SAME");
    }

    const bool legend_top_left = (strstr(saveName, "xpom_effcorr") != nullptr) ||
                                 (strstr(saveName, "t_effcorr") != nullptr);
    TLegend* legend = legend_top_left
        ? new TLegend(0.15, 0.7, 0.45, 0.9)
        : new TLegend(0.55, 0.65, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h_truth_plot, "Truth (Pseudo-data)", "l");
    legend->AddEntry(h_pdata_b0, Form("B0 (%s) uncorrected", recoMethodLabel), "p");
    legend->AddEntry(h_corr_b0, Form("B0 (%s) corrected", recoMethodLabel), "p");
    legend->AddEntry(h_pdata_rp, Form("RP (%s) uncorrected", recoMethodLabel), "p");
    legend->AddEntry(h_corr_rp, Form("RP (%s) corrected", recoMethodLabel), "p");
    if (includeSum) {
        legend->AddEntry(h_sum_uncorr, Form("B0+RP (%s) uncorrected", recoMethodLabel), "p");
        legend->AddEntry(h_sum_corr, Form("B0+RP (%s) corrected", recoMethodLabel), "p");
    }
    legend->Draw();

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.2, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
    delete h_eff_b0;
    delete h_eff_rp;
    delete h_corr_b0;
    delete h_corr_rp;
    delete h_truth_plot;
    delete h_sum_uncorr;
    delete h_sum_corr;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <combined.root>" << std::endl;
        return 1;
    }

    gErrorIgnoreLevel = kWarning;

    TString inputFileName = argv[1];

    gStyle->SetOptStat(0);
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

    std::vector<PlotOptions*> plots;
    PlotOptions1D* plot_ptr = nullptr;
    PlotOptionsBinnedRelRes* binned_plot_ptr = nullptr;

    // Create output directories
    gSystem->mkdir("figs", kTRUE);
    gSystem->mkdir("figs/distributions", kTRUE);
    gSystem->mkdir("figs/response_matrices", kTRUE);
    gSystem->mkdir("figs/resolutions/simple", kTRUE);
    gSystem->mkdir("figs/resolutions/binned", kTRUE);
    gSystem->mkdir("figs/resolutions/binned/bins", kTRUE);
    gSystem->mkdir("figs/resolutions/2d_maps", kTRUE);
    gSystem->mkdir("figs/cross_sections", kTRUE);
    gSystem->mkdir("figs/performance", kTRUE);
    gSystem->mkdir("figs/inclusive", kTRUE);
    gSystem->mkdir("figs/inclusive/histos", kTRUE);
    gSystem->mkdir("figs/inclusive/resolution", kTRUE);
    gSystem->mkdir("figs/inclusive/resolution/profile", kTRUE);
    gSystem->mkdir("figs/inclusive/response", kTRUE);
    gSystem->mkdir("figs/inclusive/efficiency", kTRUE);
    gSystem->mkdir("figs/diffractive", kTRUE);
    gSystem->mkdir("figs/diffractive/histos", kTRUE);
    gSystem->mkdir("figs/diffractive/resolution", kTRUE);
    gSystem->mkdir("figs/diffractive/resolution/profile", kTRUE);
    gSystem->mkdir("figs/diffractive/response", kTRUE);
    gSystem->mkdir("figs/diffractive/efficiency", kTRUE);

    // Acceptance/purity plots (if tracking histograms exist)
    TH1D* h_gen_Q2 = (TH1D*)inputFile->Get("h_gen_Q2");
    TH1D* h_gen_and_reco_after_cuts_Q2_EM = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_EM");
    TH1D* h_gen_and_reco_after_cuts_Q2_DA = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_DA");
    TH1D* h_gen_and_reco_after_cuts_Q2_Sigma = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_Sigma");
    TH1D* h_reco_Q2_EM = (TH1D*)inputFile->Get("h_reco_Q2_EM");
    TH1D* h_reco_Q2_DA = (TH1D*)inputFile->Get("h_reco_Q2_DA");
    TH1D* h_reco_Q2_Sigma = (TH1D*)inputFile->Get("h_reco_Q2_Sigma");

    PlotXQ2Density(inputFile);
    PlotEfficiencyCorrectedComparisonBR(inputFile,
                                        "t_truth_mc_B0",
                                        "t_reco_mc_B0",
                                        "t_reco_pdata_B0",
                                        "t_truth_mc_RP",
                                        "t_reco_mc_RP",
                                        "t_reco_pdata_RP",
                                        "t_truth_pdata_B0",
                                        "t_truth_pdata_RP",
                                        "t_truth_pdata_all",
                                        "|t| [GeV^{2}]",
                                        "|t| Efficiency Correction",
                                        "figs/diffractive/efficiency/t_effcorr.png",
                                        true,
                                        false,
                                        "Proton",
                                        false);
    PlotEfficiencyCorrectedComparisonBR(inputFile,
                                        "beta_truth_mc_B0",
                                        "beta_reco_mc_B0",
                                        "beta_reco_pdata_B0",
                                        "beta_truth_mc_RP",
                                        "beta_reco_mc_RP",
                                        "beta_reco_pdata_RP",
                                        "beta_truth_pdata_B0",
                                        "beta_truth_pdata_RP",
                                        "beta_truth_pdata_all",
                                        "#beta",
                                        "#beta Efficiency Correction (EM)",
                                        "figs/diffractive/efficiency/beta_effcorr.png",
                                        false,
                                        false,
                                        "EM",
                                        true);
    PlotEfficiencyCorrectedComparisonBR(inputFile,
                                        "xpom_truth_mc_B0",
                                        "xpom_reco_mc_B0",
                                        "xpom_reco_pdata_B0",
                                        "xpom_truth_mc_RP",
                                        "xpom_reco_mc_RP",
                                        "xpom_reco_pdata_RP",
                                        "xpom_truth_pdata_B0",
                                        "xpom_truth_pdata_RP",
                                        "xpom_truth_pdata_all",
                                        "x_{pom}",
                                        "x_{pom} Efficiency Correction (EM)",
                                        "figs/diffractive/efficiency/xpom_effcorr.png",
                                        true,
                                        false,
                                        "EM",
                                        true);
    PlotEfficiencyCorrectedComparisonBR(inputFile,
                                        "theta_truth_mc_B0",
                                        "theta_reco_mc_B0",
                                        "theta_reco_pdata_B0",
                                        "theta_truth_mc_RP",
                                        "theta_reco_mc_RP",
                                        "theta_reco_pdata_RP",
                                        "theta_truth_pdata_B0",
                                        "theta_truth_pdata_RP",
                                        "theta_truth_pdata_all",
                                        "#theta [mrad]",
                                        "#theta Efficiency Correction",
                                        "figs/diffractive/efficiency/theta_effcorr.png",
                                        false,
                                        true,
                                        "Proton",
                                        false);

    PlotEfficiencyCorrectedComparison(inputFile,
                                      "Q2_truth_mc",
                                      "Q2_reco_mc",
                                      "Q2_reco_pdata",
                                      "Q2_truth_pdata",
                                      "Q^{2} [GeV^{2}]",
                                      "Q^{2} Efficiency Correction (EM)",
                                      "figs/inclusive/efficiency/q2_effcorr.png",
                                      true,
                                      "EM");
    PlotEfficiencyCorrectedComparison(inputFile,
                                      "x_truth_mc",
                                      "x_reco_mc",
                                      "x_reco_pdata",
                                      "x_truth_pdata",
                                      "x_{Bj}",
                                      "x_{Bj} Efficiency Correction (EM)",
                                      "figs/inclusive/efficiency/xbj_effcorr.png",
                                      true,
                                      "EM");
    PlotEfficiencyCorrectedComparison(inputFile,
                                      "y_truth_mc",
                                      "y_reco_mc",
                                      "y_reco_pdata",
                                      "y_truth_pdata",
                                      "y",
                                      "y Efficiency Correction (EM)",
                                      "figs/inclusive/efficiency/y_effcorr.png",
                                      false,
                                      "EM");
    PlotEfficiencyCorrectedComparison(inputFile,
                                      "W2_truth_mc",
                                      "W2_reco_mc",
                                      "W2_reco_pdata",
                                      "W2_truth_pdata",
                                      "W^{2} [GeV^{2}]",
                                      "W^{2} Efficiency Correction (EM)",
                                      "figs/inclusive/efficiency/w2_effcorr.png",
                                      true,
                                      "EM");
    PlotEfficiencyCorrectedComparison(inputFile,
                                      "MX2_truth_mc",
                                      "MX2_reco_mc",
                                      "MX2_reco_pdata",
                                      "MX2_truth_pdata",
                                      "M_{X}^{2} [GeV^{2}]",
                                      "M_{X}^{2} Efficiency Correction",
                                      "figs/diffractive/efficiency/mx2_effcorr.png",
                                      true,
                                      "Hadrons");
    PlotDensityFromHist(inputFile, "beta_Q2_reco", "#beta", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/beta_q2_density.png", false, true);
    PlotDensityFromHist(inputFile, "t_Q2_reco", "|t| [GeV^{2}]", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/t_q2_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_Q2_reco", "x_{pom}", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/xpom_q2_density.png", true, true);
    PlotDensityFromHist(inputFile, "beta_t_reco", "#beta", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/beta_t_density.png", false, true);
    PlotDensityFromHist(inputFile, "xbj_t_reco", "x_{Bj}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xbj_t_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_t_reco", "x_{pom}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xpom_t_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_beta_reco", "x_{pom}", "#beta",
                        "figs/diffractive/histos/xpom_beta_density.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_beta_reco", "x_{Bj}", "#beta",
                        "figs/diffractive/histos/xbj_beta_density.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_xpom_reco", "x_{Bj}", "x_{pom}",
                        "figs/diffractive/histos/xbj_xpom_density.png", true, true);
    PlotPhaseSpaceSlices(inputFile);

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

    // =================================================================
    // Binned resolution plots (Q2/x/y)
    // =================================================================
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_EM",
        "Relative bin by bin resolution (EM);Q^{2}_{MC};#frac{Q^{2}_{EM} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{EM}",
        "",
        {
         {-0.0, 0.0}, {-0.022, 0.02}, {-0.02, 0.02}, {-0.02, 0.02}, {-0.02, 0.02},
         {-0.015, 0.015}, {-0.015, 0.015}, {-0.014, 0.015}, {-0.025, 0.025}, {-0.01, 0.012},
         {-0.027, 0.028}, {-0.018, 0.02}, {-0.022, 0.02}, {-0.02, 0.015}, {-0.018, 0.02},
         {-0.02, 0.015}, {-0.02, 0.017}, {-0.017, 0.02}, {-0.02, 0.02}, {-0.04, 0.04},
         {-0.025, 0.03}, {-0.015, 0.025}, {-0.05, 0.06}
        },
        "figs/inclusive/resolution/q2_relres_binned_em.png",
        "inclusive/resolution/profile/q2_relres_binned_em",
        std::make_pair(5.0, 200),
        true
    ));

    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_DA",
        "Relative bin by bin resolution (DA);Q^{2}_{MC};#frac{Q^{2}_{DA} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{DA}",
        "",
        {
          {-0.0, 0.0}, {-0.005, 0.025}, {-0.01, 0.025}, {-0.006, 0.02}, {-0.01, 0.02},
          {-0.009, 0.02}, {-0, 0}, {-0., 0.}, {-0., 0.}, {-0., 0.},
          {-0.015, 0.03}, {-0.009, 0.02}, {-0.01, 0.02}, {-0.01, 0.02}, {-0.01, 0.02},
          {-0.01, 0.02}, {-0.004, 0.02}, {-0.017, 0.027}, {-0.025, 0.03}, {-0.08, 0.08},
          {-0.05, 0.06}, {-0.05, 0.065}, {-0.05, 0.06}
        },
        "figs/inclusive/resolution/q2_relres_binned_da.png",
        "inclusive/resolution/profile/q2_relres_binned_da",
        std::make_pair(5.0, 200),
        true
    ));

    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_Sigma",
        ";Q^{2}_{MC};#frac{Q^{2}_{#Sigma} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{#Sigma}",
        "",
        {
         {-0.0, 0.0}, {-0.022, 0.02}, {-0.02, 0.02}, {-0.02, 0.02}, {-0.02, 0.02},
         {-0.015, 0.015}, {-0.015, 0.015}, {-0.014, 0.015}, {-0.025, 0.025}, {-0.01, 0.012},
         {-0.027, 0.028}, {-0.018, 0.02}, {-0.022, 0.02}, {-0.02, 0.015}, {-0.018, 0.02},
         {-0.02, 0.015}, {-0.02, 0.017}, {-0.017, 0.02}, {-0.02, 0.02}, {-0.04, 0.04},
         {-0.025, 0.03}, {-0.015, 0.025}, {-0.05, 0.06}
        },
        "figs/inclusive/resolution/q2_relres_binned_sigma.png",
        "inclusive/resolution/profile/q2_relres_binned_sigma",
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
         {-0.015, 0.015}, {-0.015, 0.015}, {0, 0}, {-0.025, 0.025}, {-0.025, 0.02},
         {-0.027, 0.028}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
         {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
         {0, 0}, {0, 0}, {0, 0}
        },
        "figs/inclusive/resolution/xbj_relres_binned_em.png",
        "inclusive/resolution/profile/xbj_relres_binned_em",
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
        "inclusive/resolution/profile/xbj_relres_binned_da",
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
        "inclusive/resolution/profile/xbj_relres_binned_sigma",
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
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0}
        },
        "figs/inclusive/resolution/y_relres_binned_em.png",
        "inclusive/resolution/profile/y_relres_binned_em",
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
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0}
        },
        "figs/inclusive/resolution/y_relres_binned_da.png",
        "inclusive/resolution/profile/y_relres_binned_da",
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
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0}
        },
        "figs/inclusive/resolution/y_relres_binned_sigma.png",
        "inclusive/resolution/profile/y_relres_binned_sigma",
        std::make_pair(0.0, 1.0),
        false
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
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
        {"dsigma_dt_MC", "dsigma_dt_B0", "dsigma_dt_RP"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "pe", "pe"},
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
    // =================================================================
    // Execute all standard plots
    // =================================================================
    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    for (const auto& plot : plots) {
        delete plot;
    }

    inputFile->Close();
    delete inputFile;

    return 0;
}
