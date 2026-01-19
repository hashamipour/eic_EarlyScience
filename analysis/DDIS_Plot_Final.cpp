// Unified DDIS plotter: union of all plots (including previously commented) + optional binning plots
// g++ DDIS_Plot_Final.cpp -o DDIS_Plot_Final $(root-config --cflags --glibs)
// ./DDIS_Plot_Final <combined.root> [binning_inputs.root] [binning_scheme.txt] [binning_output.root]

#include "Plotting.hpp"
#include "BinningUtility.hpp"

#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TString.h>
#include <TF1.h>
#include <TLegend.h>
#include <TBox.h>
#include <TLine.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

// -------------------------
// Binning plot helpers
// -------------------------
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

static CombinedHists BookCombinedHists(const BinningScheme& binning) {
    CombinedHists h;
    const int nBins = binning.nTotalBins();

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
    for (int k = 0; k < nBins; k++) {
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
    for (int i = 1; i <= nBins; i++) {
        double truth = h.truth_count[i - 1];
        for (int j = 1; j <= nBins; j++) {
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
    for (int k = 0; k < nBins; k++) {
        int iQ2 = -1, iBeta = -1, iXpom = -1;
        binning.Get3DIndices(k, iQ2, iBeta, iXpom);
        if (iQ2 < 0 || iQ2 >= binning.nQ2()) continue;
        truth_Q2[iQ2] += h.truth_count[k];
        truth_reco_Q2[iQ2] += h.truth_with_reco[k];
        reco_Q2[iQ2] += h.reco_total[k];
        diag_Q2[iQ2] += h.diag_count[k];
    }

    for (int iQ2 = 0; iQ2 < binning.nQ2(); iQ2++) {
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
    for (int k = 0; k < nBins; k++) {
        int iQ2 = -1, iBeta = -1, iXpom = -1;
        binning.Get3DIndices(k, iQ2, iBeta, iXpom);
        if (iBeta < 0 || iBeta >= nBeta || iXpom < 0 || iXpom >= nXpom) continue;
        int idx = iBeta * nXpom + iXpom;
        truth[idx] += h.truth_count[k];
        truth_reco[idx] += h.truth_with_reco[k];
        reco[idx] += h.reco_total[k];
        diag[idx] += h.diag_count[k];
    }

    for (int iBeta = 0; iBeta < nBeta; iBeta++) {
        for (int iXpom = 0; iXpom < nXpom; iXpom++) {
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
    if (bin < 1) return 1;
    if (bin > nBins) return nBins;
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

    for (size_t ix = 0; ix + 1 < x_edges.size(); ix++) {
        double x1 = x_edges[ix];
        double x2 = x_edges[ix + 1];
        int bx1 = FindBinLow(xaxis, x1);
        int bx2 = FindBinHigh(xaxis, x2);

        for (size_t iy = 0; iy + 1 < y_edges.size(); iy++) {
            double y1 = y_edges[iy];
            double y2 = y_edges[iy + 1];
            int by1 = FindBinLow(yaxis, y1);
            int by2 = FindBinHigh(yaxis, y2);

            double count = h->Integral(bx1, bx2, by1, by2);
            if (count < minCount) continue;

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

    for (size_t ix = 0; ix + 1 < x_edges.size(); ix++) {
        double x1 = x_edges[ix];
        double x2 = x_edges[ix + 1];
        int bx1 = FindBinLow(xaxis, x1);
        int bx2 = FindBinHigh(xaxis, x2);

        for (size_t iy = 0; iy + 1 < y_edges.size(); iy++) {
            double y1 = y_edges[iy];
            double y2 = y_edges[iy + 1];
            int by1 = FindBinLow(yaxis, y1);
            int by2 = FindBinHigh(yaxis, y2);

            double count = h->Integral(bx1, bx2, by1, by2);
            if (count < minCount) continue;

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

    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (Q2_truth > 0 && xpom_truth > 0) {
            h_Q2_xpom->Fill(TMath::Log10(xpom_truth), Q2_truth);
        }
        if (Q2_truth > 0 && beta_truth > 0) {
            h_Q2_beta->Fill(beta_truth, Q2_truth);
        }
        if (beta_truth > 0 && xpom_truth > 0) {
            h_beta_xpom->Fill(TMath::Log10(xpom_truth), beta_truth);
        }
        if (Q2_truth > 0 && y_truth > 0) {
            h_Q2_y->Fill(y_truth, Q2_truth);
        }
    }

    const double minCount = 100.0;

    std::vector<double> xpom_log_edges;
    xpom_log_edges.reserve(binning.xpom_edges.size());
    for (double edge : binning.xpom_edges) {
        xpom_log_edges.push_back(TMath::Log10(edge));
    }

    {
        std::ofstream out("data/binning_scheme_cells.txt", std::ios::trunc);
        if (out.is_open()) {
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

static void RunBinningPlots(const std::string& inputFile,
                            const std::string& binningFile,
                            const std::string& outputFile) {
    BinningScheme binning = ReadBinningScheme(binningFile);
    if (binning.Q2_edges.empty() || binning.beta_edges.empty() || binning.xpom_edges.empty()) {
        std::cerr << "ERROR: Invalid binning scheme: " << binningFile << std::endl;
        return;
    }
    PrintBinningScheme(binning);

    TFile* infile = TFile::Open(inputFile.c_str(), "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "ERROR: Cannot open binning input file " << inputFile << std::endl;
        return;
    }

    TTree* tree = (TTree*)infile->Get("binning");
    if (!tree) {
        std::cerr << "ERROR: Cannot find 'binning' tree in file" << std::endl;
        infile->Close();
        return;
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
    CombinedHists h_all = BookCombinedHists(binning);

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        MethodHists* h = (method == 0) ? &h_b0 : &h_rp;

        int iQ2_true = -1;
        int iBeta_true = -1;
        int iXpom_true = -1;
        int iQ2_reco = -1;
        int iBeta_reco = -1;
        int iXpom_reco = -1;

        if (Q2_truth > 0 && beta_truth > 0 && xpom_truth > 0) {
            iQ2_true = binning.GetQ2Bin(Q2_truth);
            iBeta_true = binning.GetBetaBin(beta_truth);
            iXpom_true = binning.GetXpomBin(xpom_truth);
        }
        if (Q2_reco > 0 && beta_reco > 0 && xpom_reco > 0) {
            iQ2_reco = binning.GetQ2Bin(Q2_reco);
            iBeta_reco = binning.GetBetaBin(beta_reco);
            iXpom_reco = binning.GetXpomBin(xpom_reco);
        }

        int k_true = (iQ2_true >= 0 && iBeta_true >= 0 && iXpom_true >= 0)
            ? binning.GetLinearBin(iQ2_true, iBeta_true, iXpom_true)
            : -1;
        int k_reco = (iQ2_reco >= 0 && iBeta_reco >= 0 && iXpom_reco >= 0)
            ? binning.GetLinearBin(iQ2_reco, iBeta_reco, iXpom_reco)
            : -1;

        if (k_true >= 0) h_all.truth_count[k_true] += 1.0;
        if (k_reco >= 0) h_all.reco_total[k_reco] += 1.0;
        if (k_true >= 0 && k_reco >= 0) {
            h_all.truth_with_reco[k_true] += 1.0;
            h_all.h_migration->Fill(k_true + 1, k_reco + 1);
            if (k_true == k_reco) {
                h_all.diag_count[k_true] += 1.0;
            }
        }

        if (k_true >= 0 && k_reco >= 0) {
            if (Q2_truth != 0.0f) {
                double res_Q2 = (Q2_reco - Q2_truth) / Q2_truth;
                h_all.h_res_Q2->Fill(Q2_truth, res_Q2);
            }
            if (beta_truth != 0.0f) {
                double res_beta = (beta_reco - beta_truth) / beta_truth;
                h->h_res_beta->Fill(beta_truth, res_beta);

                if (xpom_truth != 0.0f) {
                    double res_xpom = (xpom_reco - xpom_truth) / xpom_truth;
                    h->h_res_xpom->Fill(xpom_truth, res_xpom);
                    h->h_res_beta_vs_xpom->Fill(res_xpom, res_beta);
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

    outfile->Close();

    const std::string outdir = "figs/binning";
    const std::string jpgdir = "figs/binning/jpg";
    gSystem->mkdir(outdir.c_str(), true);
    gSystem->mkdir(jpgdir.c_str(), true);
    SavePlots(h_b0, outdir, jpgdir);
    SavePlots(h_rp, outdir, jpgdir);
    SaveCombinedPlots(h_all, outdir, jpgdir);
    SaveBinningSchemePlots(tree, binning, outdir, jpgdir);

    infile->Close();
    std::cout << "Wrote binning plots to " << outdir << " and histograms to " << outputFile << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <combined.root> [binning_inputs.root] [binning_scheme.txt] [binning_output.root]" << std::endl;
        return 1;
    }

    gErrorIgnoreLevel = kWarning;

    TString inputFileName = argv[1];
    std::string binningInput = (argc > 2) ? argv[2] : "";
    std::string binningScheme = (argc > 3) ? argv[3] : "data/binning_scheme.txt";
    std::string binningOutput = (argc > 4) ? argv[4] : "DDIS_BinningPlots.root";

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

    // Acceptance/purity plots (if tracking histograms exist)
    TH1D* h_gen_Q2 = (TH1D*)inputFile->Get("h_gen_Q2");
    TH1D* h_gen_and_reco_after_cuts_Q2_EM = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_EM");
    TH1D* h_gen_and_reco_after_cuts_Q2_DA = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_DA");
    TH1D* h_gen_and_reco_after_cuts_Q2_Sigma = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_Sigma");
    TH1D* h_reco_Q2_EM = (TH1D*)inputFile->Get("h_reco_Q2_EM");
    TH1D* h_reco_Q2_DA = (TH1D*)inputFile->Get("h_reco_Q2_DA");
    TH1D* h_reco_Q2_Sigma = (TH1D*)inputFile->Get("h_reco_Q2_Sigma");

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
        "figs/distributions/Q2_hist.png",
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
        "figs/distributions/Q2_pdf.png",
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
        "figs/distributions/x_hist.png",
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
        "figs/distributions/x_pdf.png",
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
        "figs/distributions/y_hist.png",
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
        "figs/distributions/y_pdf.png",
        false,
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // =================================================================
    // Overall resolution plots (Q2/x/y)
    // =================================================================
    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_EM",
        "#frac{Q^{2}_{EM} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Counts",
        -0.01, 0.01,
        "figs/resolutions/simple/DDIS_Q2RelRes_EM.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_DA",
        "#frac{Q^{2}_{DA} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Counts",
        -0.005, 0.02,
        "figs/resolutions/simple/DDIS_Q2RelRes_DA.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_Sigma",
        "#frac{Q^{2}_{#Sigma} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Counts",
        -0.01, 0.01,
        "figs/resolutions/simple/DDIS_Q2RelRes_Sigma.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_EM",
        "#frac{x_{EM} - x_{MC}}{ x_{MC}}",
        "Counts",
        -0.025, 0.02,
        "figs/resolutions/simple/DDIS_RelRes_xBj_EM.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_DA",
        "#frac{x_{DA} - x_{MC}}{X_{MC}}",
        "Counts",
        0., 0.,
        "figs/resolutions/simple/DDIS_RelRes_xBj_DA.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_Sigma",
        "#frac{x_{#Sigma} - x_{MC}}{X_{MC}}",
        "Counts",
        0., 0.,
        "figs/resolutions/simple/DDIS_RelRes_x_Sigma.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_EM",
        "#frac{y_{EM} - y_{MC}}{ y_{MC}}",
        "Counts",
        -0.009, 0.009,
        "figs/resolutions/simple/DDIS_RelRes_y_EM.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_DA",
        "#frac{y_{DA} - y_{MC}}{y_{MC}}",
        "Counts",
        0., 0.,
        "figs/resolutions/simple/DDIS_RelRes_y_DA.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_Sigma",
        "#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "Counts",
        0., 0.,
        "figs/resolutions/simple/DDIS_RelRes_y_Sigma.png"
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
        "figs/resolutions/binned/DDIS_Q2RelRes_binned_EM.png",
        "DDIS_Q2RelRes_binned_EM",
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
        "figs/resolutions/binned/DDIS_Q2RelRes_binned_DA.png",
        "DDIS_Q2RelRes_binned_DA",
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
        "figs/resolutions/binned/DDIS_Q2RelRes_binned_Sigma.png",
        "DDIS_Q2RelRes_binned_Sigma",
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
        "figs/resolutions/binned/DDIS_RelRes_binned_x_EM.png",
        "DDIS_RelRes_binned_x_EM",
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
        "figs/resolutions/binned/DDIS_RelRes_binned_x_DA.png",
        "DDIS_RelRes_binned_x_DA",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.75, 0.1, 0.9, 0.25);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_Sigma",
        ";x_{MC};#frac{x_{#Sigma} - x_{MC}}{x_{MC}}",
        "x_Sigma",
        "",
        {
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0}
        },
        "figs/resolutions/binned/DDIS_RelRes_binned_x_Sigma.png",
        "DDIS_RelRes_binned_x_Sigma",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.75, 0.3, 0.9);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_EM",
        ";y_{MC};#frac{y_{EM} - y_{MC}}{y_{MC}}",
        "y_{EM}",
        "",
        {
            {-0.,0.},{-0.,0.},{-0.,0.},{-0.08,0.08},
            {-0.06,0.07},{-0.04,0.045},{-0.025,0.025},{-0.03,0.035},
            {-0.03,0.036},{-0.02,0.02},{-0.018,0.018},{-0.015,0.015},
            {-0.015,0.015},{-0.01,0.015},{-0.01,0.01},{-0.011,0.011}
        },
        "figs/resolutions/binned/DDIS_RelRes_binned_y_EM.png",
        "DDIS_RelRes_binned_y_EM",
        std::make_pair(0.0, 1.0)
    );
    binned_plot_ptr->SetLegendPosition(0.75, 0.75, 0.9, 0.9);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_DA",
        ";y_{MC};#frac{y_{DA} - y_{MC}}{y_{MC}}",
        "y_{DA}",
        "",
        {
            {0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0},
            {0,0},{0,0},{0,0},{0,0}
        },
        "figs/resolutions/binned/DDIS_RelRes_binned_y_DA.png",
        "DDIS_RelRes_binned_y_DA"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_Sigma",
        ";y_{MC};#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "y_Sigma",
        "",
        {
            {0,0},{0,0},{-0,0.},{-0,0},
            {0,0},{0,0},{-0,0.},{-0,0},
            {0,0},{0,0},{-0,0.},{-0,0},
            {0,0},{0,0},{-0,0.},{-0,0}
        },
        "figs/resolutions/binned/DDIS_RelRes_binned_y_Sigma.png",
        "DDIS_RelRes_binned_y_Sigma"
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.75, 0.3, 0.9);
    plots.push_back(binned_plot_ptr);

    // =================================================================
    // Response matrices (Q2/x/y)
    // =================================================================
    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_EM",
        "Q^{2} (true) [GeV]",
        "Q^{2} (EM) [GeV]",
        "figs/response_matrices/response_matrix_Q2_EM.png",
        true,
        true,
        {1.0, 300},
        {1.0, 300}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_DA",
        "Q^{2} (true) [GeV]",
        "Q^{2} (DA) [GeV]",
        "figs/response_matrices/response_matrix_Q2_DA.png",
        true,
        true,
        {1.0, 300},
        {1.0, 300}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_Sigma",
        "Q^{2} (true) [GeV]",
        "Q^{2} (Sigma) [GeV]",
        "figs/response_matrices/response_matrix_Q2_Sigma.png",
        true,
        true,
        {1.0, 300},
        {1.0, 300}
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
        {1e-3, 0.5}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_EM",
        "y (true)",
        "y (EM)",
        "figs/response_matrices/response_matrix_y_EM.png",
        false,
        false,
        {0., 1.},
        {0., 1.}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_DA",
        "y (true)",
        "y (DA)",
        "figs/response_matrices/response_matrix_y_DA.png",
        false,
        false,
        {0., 1.},
        {0., 1.}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_Sigma",
        "y (true)",
        "y (Sigma)",
        "figs/response_matrices/response_matrix_y_Sigma.png",
        false,
        false,
        {0., 1.},
        {0., 1.}
    ));

    // =================================================================
    // Hadronic final state distributions
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Hadronic Final State E-p_{z}",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/distributions/EPz_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Hadronic Final State E-p_{z}",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/distributions/EPz_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Counts",
        "figs/distributions/eta_max_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Counts",
        "figs/distributions/eta_max_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.3, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_MX2_truth", "h_MX2"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Hadronic Invariant Mass Squared",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/distributions/MX2_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"h_MX2_truth", "h_MX2"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Hadronic Invariant Mass Squared",
        "M_{X}^{2} [GeV^{2}]",
        "Counts",
        "figs/distributions/MX2_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    // =================================================================
    // Mandelstam t plots (all previously commented now enabled)
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| Distributions",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/distributions/t_distributions.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| Distributions",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/distributions/t_distributions_logy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
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

    plot_ptr = new PlotOptions1D(
        {"theta_MC", "theta_B0", "theta_RP"},
        {"MC Truth", "Reco B0", "Reco RP"},
        {"hist", "pe", "pe"},
        "Proton Scattering Angles",
        "#theta [mrad]",
        "Counts",
        "figs/distributions/theta_distributions.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsRelRes(
        "t_res_B0",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/t_resolution_B0.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "t_res_RP",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Counts",
        -0.5, 0.5,
        "figs/resolutions/simple/t_resolution_RP.png"
    ));

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

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"t_corr_B0_graph", "t_corr_RP_graph"},
        {"B0 Protons (5.5-20 mrad)", "Roman Pot (<5 mrad)"},
        {kBlue, kCyan + 1},
        {20, 21},
        "Combined Truth vs Reco |t| Correlation",
        "Truth |t| [GeV^{2}]",
        "Reco |t| [GeV^{2}]",
        "figs/response_matrices/t_correlation_combined.png",
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}]",
        "figs/response_matrices/response_matrix_t_B0.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_RP",
        "Truth |t| [GeV^{2}]",
        "RP Reco |t| [GeV^{2}]",
        "figs/response_matrices/response_matrix_t_RP.png",
        true,
        true,
        {0.0, 0.5},
        {0.0, 0.5}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_B0",
        "Truth x_L",
        "B0 Reco x_L",
        "figs/response_matrices/response_matrix_xL_B0.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_RP",
        "Truth x_L",
        "RP Reco x_L",
        "figs/response_matrices/response_matrix_xL_RP.png",
        false,
        false,
        {0.0, 0.5},
        {0.0, 0.5}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0",
        "Truth x_{pom}",
        "B0 Reco x_{pom}",
        "figs/response_matrices/response_matrix_xpom_B0.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_RP",
        "Truth x_{pom}",
        "RP Reco x_{pom}",
        "figs/response_matrices/response_matrix_xpom_RP.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    // x_L comparisons
    plot_ptr = new PlotOptions1D(
        {"xL_MC", "xL_B0", "xL_RP"},
        {"Truth", "B0", "Roman Pot"},
        {"hist", "hist", "hist"},
        "x_{L}",
        "x_{L}",
        "Counts",
        "figs/distributions/x_L_comparison.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xL_MC", "xL_B0", "xL_RP"},
        {"Truth", "B0", "Roman Pot"},
        {"hist", "hist", "hist"},
        "x_{L}",
        "x_{L}",
        "Counts",
        "figs/distributions/x_L_comparison_logy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    plots.push_back(plot_ptr);

    // x_pom comparisons
    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_B0", "xpom_RP"},
        {"Truth", "B0", "Roman Pot"},
        {"hist", "hist", "hist"},
        "x_{pom}",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_logxy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_B0", "xpom_RP"},
        {"Truth", "B0", "Roman Pot"},
        {"hist", "hist", "hist"},
        "x_{pom}",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_logx.png",
        true,
        false
    );
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
    plots.push_back(plot_ptr);

    // Beta distributions (with sum) drawn manually later

    plots.push_back(new PlotOptionsRelRes(
        "beta_res_B0",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Counts",
        -999., -999.,
        "figs/resolutions/simple/beta_resolution_B0.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "beta_res_RP",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Counts",
        -999., -999.,
        "figs/resolutions/simple/beta_resolution_RP.png"
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
        "beta_vs_Q2",
        "Q^{2} [GeV^{2}]",
        "#beta",
        "figs/distributions/beta_vs_Q2.png",
        true,
        false,
        {0.1, 100.0},
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
        {0.0, 2.0},
        {0.0, 1.0}
    ));

    // B0 cut-first-bin comparisons
    plot_ptr = new PlotOptions1D(
        {"t_B0", "t_B0_cutFirstBin"},
        {"B0 (all bins)", "B0 (first bin cut)"},
        {"hist", "hist"},
        "B0 |t| Distribution: Effect of First Bin Cut",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/distributions/t_B0_comparison_firstBinCut.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xL_B0", "xL_B0_cutFirstBin"},
        {"B0 (all bins)", "B0 (first bin cut)"},
        {"hist", "hist"},
        "B0 x_{L} Distribution: Effect of First Bin Cut",
        "x_{L}",
        "Counts",
        "figs/distributions/xL_B0_comparison_firstBinCut.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_B0", "xpom_B0_cutFirstBin"},
        {"B0 (all bins)", "B0 (first bin cut)"},
        {"hist", "hist"},
        "B0 x_{pom} Distribution: Effect of First Bin Cut",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_B0_comparison_firstBinCut.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.15, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0_cutFirstBin",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}] (first bin cut)",
        "figs/response_matrices/response_matrix_B0_cutFirstBin.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_B0_cutFirstBin",
        "Truth x_L",
        "B0 Reco x_L (first bin cut)",
        "figs/response_matrices/response_matrix_xL_B0_cutFirstBin.png",
        false,
        false,
        {0.0, 2.0},
        {0.0, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0_cutFirstBin",
        "Truth x_{pom}",
        "B0 Reco x_{pom} (first bin cut)",
        "figs/response_matrices/response_matrix_xpom_B0_cutFirstBin.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    // x_pom comparison: definition vs x_L
    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_def_MC"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
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
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
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
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
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

    // x_pom sum overlay drawn manually later

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_MC",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_MC.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_B0",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_B0.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_RP",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_RP.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    // Theta acceptance plots
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

    // Binned resolution for diffractive variables
    PlotOptionsBinnedRelRes* plot_t_B0 = new PlotOptionsBinnedRelRes(
        "t_RelRes_binned_B0",
        "B0 |t| Resolution vs Truth |t|",
        "|t|_{truth} [GeV^{2}]",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/t_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/t_B0",
        {0.0, 2.0},
        false
    );
    plot_t_B0->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_t_B0);

    PlotOptionsBinnedRelRes* plot_t_RP = new PlotOptionsBinnedRelRes(
        "t_RelRes_binned_RP",
        "RP |t| Resolution vs Truth |t|",
        "|t|_{truth} [GeV^{2}]",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/t_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/t_RP",
        {0.0, 0.5},
        false
    );
    plot_t_RP->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_t_RP);

    PlotOptionsBinnedRelRes* plot_xL_B0 = new PlotOptionsBinnedRelRes(
        "xL_RelRes_binned_B0",
        "B0 x_{L} Resolution vs Truth x_{L}",
        "x_{L,truth}",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/xL_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/xL_B0",
        {0.75, 1.05},
        false
    );
    plot_xL_B0->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_xL_B0);

    PlotOptionsBinnedRelRes* plot_xL_RP = new PlotOptionsBinnedRelRes(
        "xL_RelRes_binned_RP",
        "RP x_{L} Resolution vs Truth x_{L}",
        "x_{L,truth}",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/xL_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/xL_RP",
        {0.75, 1.05},
        false
    );
    plot_xL_RP->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_xL_RP);

    PlotOptionsBinnedRelRes* plot_xpom_B0 = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_B0",
        "B0 x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/xpom_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/xpom_B0",
        {1e-4, 0.4},
        true
    );
    plot_xpom_B0->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_xpom_B0);

    PlotOptionsBinnedRelRes* plot_xpom_RP = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_RP",
        "RP x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/xpom_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/xpom_RP",
        {1e-4, 0.4},
        true
    );
    plot_xpom_RP->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_xpom_RP);

    PlotOptionsBinnedRelRes* plot_beta_B0 = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_B0",
        "B0 #beta Resolution vs Truth #beta",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/beta_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/beta_B0",
        {0.0, 1.0},
        false
    );
    plot_beta_B0->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_beta_B0);

    PlotOptionsBinnedRelRes* plot_beta_RP = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_RP",
        "RP #beta Resolution vs Truth #beta",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.3, 0.3}},
        "figs/resolutions/binned/beta_resolution_binned_RP.png",
        "figs/resolutions/binned/bins/beta_RP",
        {0.0, 1.0},
        false
    );
    plot_beta_RP->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_beta_RP);

    PlotOptionsBinnedRelRes* plot_MX2 = new PlotOptionsBinnedRelRes(
        "MX2_RelRes_binned",
        "M_{X}^{2} Resolution vs Truth M_{X}^{2}",
        "M_{X,truth}^{2} [GeV^{2}]",
        "(M_{X,reco}^{2} - M_{X,truth}^{2})/M_{X,truth}^{2}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/MX2_resolution_binned.png",
        "figs/resolutions/binned/bins/MX2",
        {0.0, 200.0},
        false
    );
    plot_MX2->SetLegendPosition(0.6, 0.7, 0.85, 0.85);
    plots.push_back(plot_MX2);

    // Differential cross sections
    TH1D* h_dsigma_dt_B0_temp = (TH1D*)inputFile->Get("dsigma_dt_B0");
    TH1D* h_dsigma_dt_RP_temp = (TH1D*)inputFile->Get("dsigma_dt_RP");
    if (h_dsigma_dt_B0_temp && h_dsigma_dt_RP_temp) {
        TH1D* h_dsigma_dt_Sum = (TH1D*)h_dsigma_dt_B0_temp->Clone("dsigma_dt_Sum");
        h_dsigma_dt_Sum->SetTitle("B0+RP Sum d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]");
        h_dsigma_dt_Sum->Add(h_dsigma_dt_RP_temp);
    }

    plot_ptr = new PlotOptions1D(
        {"dsigma_dt_MC", "dsigma_dt_B0", "dsigma_dt_RP", "dsigma_dt_Sum"},
        {"MC Truth", "B0 Reco", "RP Reco", "B0+RP Sum"},
        {"hist", "pe", "pe", "pe"},
        "Differential Cross Section d#sigma/dt",
        "|t| [GeV^{2}]",
        "d#sigma/dt [nb/GeV^{2}]",
        "figs/cross_sections/dsigma_dt.png",
        true,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.65, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"dsigma_dt_MC", "dsigma_dt_B0", "dsigma_dt_RP", "dsigma_dt_Sum"},
        {"MC Truth", "B0 Reco", "RP Reco", "B0+RP Sum"},
        {"hist", "pe", "pe", "pe"},
        "Differential Cross Section d#sigma/dt",
        "|t| [GeV^{2}]",
        "d#sigma/dt [nb/GeV^{2}]",
        "figs/cross_sections/dsigma_dt_logy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.65, 0.3, 0.9);
    plots.push_back(plot_ptr);

    // Additional debug plot (Q2_EICRecon vs Q2_calc)
    plot_ptr = new PlotOptions1D(
        {"Q2_EICRecon", "Q2_calc"},
        {"EICRecon", "Our Calculation"},
        {"hist", "hist"},
        "Q^{2} Comparison",
        "Q^{2} [GeV^{2}]",
        "Counts",
        "figs/distributions/Q2_comparison.png",
        false,
        true
    );
    plots.push_back(plot_ptr);

    // MX2 comparison (legacy names)
    plot_ptr = new PlotOptions1D(
        {"MX2_MC", "MX2_eX"},
        {"Truth", "Reco"},
        {"hist", "hist"},
        "M_{X}^{2}",
        "M_{X}^{2}",
        "Counts",
        "figs/distributions/MX2_comparison.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.75, 0.7, 0.95, 0.9);
    plots.push_back(plot_ptr);

    // =================================================================
    // Execute all standard plots
    // =================================================================
    for (const auto& plot : plots) {
        plot->Plot(inputFile);
    }
    for (const auto& plot : plots) {
        delete plot;
    }

    // Manual overlays that require derived sum histograms
    {
        TH1D* h_beta_MC = (TH1D*)inputFile->Get("beta_MC");
        TH1D* h_beta_B0_local = (TH1D*)inputFile->Get("beta_B0");
        TH1D* h_beta_RP_local = (TH1D*)inputFile->Get("beta_RP");
        if (h_beta_MC && h_beta_B0_local && h_beta_RP_local) {
            TH1D* h_beta_Sum = (TH1D*)h_beta_B0_local->Clone("beta_Sum_manual");
            h_beta_Sum->Add(h_beta_RP_local);

            TCanvas* c_beta = new TCanvas("c_beta_sum", "beta distributions", 900, 700);
            c_beta->SetLogy();
            h_beta_MC->SetLineColor(kBlack);
            h_beta_MC->Draw("HIST");
            h_beta_B0_local->SetMarkerStyle(20);
            h_beta_B0_local->SetMarkerColor(kRed);
            h_beta_B0_local->Draw("PE SAME");
            h_beta_RP_local->SetMarkerStyle(20);
            h_beta_RP_local->SetMarkerColor(kBlue);
            h_beta_RP_local->Draw("PE SAME");
            h_beta_Sum->SetMarkerStyle(20);
            h_beta_Sum->SetMarkerColor(kOrange + 7);
            h_beta_Sum->Draw("PE SAME");

            TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.9);
            leg->AddEntry(h_beta_MC, "MC Truth", "l");
            leg->AddEntry(h_beta_B0_local, "B0 Reco", "pe");
            leg->AddEntry(h_beta_RP_local, "RP Reco", "pe");
            leg->AddEntry(h_beta_Sum, "B0+RP Sum", "pe");
            leg->Draw();

            c_beta->SaveAs("figs/distributions/beta_distributions_logy.png");
            delete c_beta;
        }
    }

    {
        TH1D* h_xpom_def_MC = (TH1D*)inputFile->Get("xpom_def_MC");
        TH1D* h_xpom_def_B0_local = (TH1D*)inputFile->Get("xpom_def_B0");
        TH1D* h_xpom_def_RP_local = (TH1D*)inputFile->Get("xpom_def_RP");
        if (h_xpom_def_MC && h_xpom_def_B0_local && h_xpom_def_RP_local) {
            TH1D* h_xpom_def_Sum = (TH1D*)h_xpom_def_B0_local->Clone("xpom_def_Sum_manual");
            h_xpom_def_Sum->Add(h_xpom_def_RP_local);

            TCanvas* c_xpom = new TCanvas("c_xpom_sum", "x_pom comparison", 900, 700);
            c_xpom->SetLogx();
            c_xpom->SetLogy();
            h_xpom_def_MC->SetLineColor(kBlack);
            h_xpom_def_MC->Draw("HIST");
            h_xpom_def_B0_local->SetMarkerStyle(20);
            h_xpom_def_B0_local->SetMarkerColor(kRed);
            h_xpom_def_B0_local->Draw("PE SAME");
            h_xpom_def_RP_local->SetMarkerStyle(20);
            h_xpom_def_RP_local->SetMarkerColor(kBlue);
            h_xpom_def_RP_local->Draw("PE SAME");
            h_xpom_def_Sum->SetMarkerStyle(20);
            h_xpom_def_Sum->SetMarkerColor(kOrange + 7);
            h_xpom_def_Sum->Draw("PE SAME");

            TLegend* leg = new TLegend(0.1, 0.7, 0.5, 0.9);
            leg->AddEntry(h_xpom_def_MC, "x_{pom} MC", "l");
            leg->AddEntry(h_xpom_def_B0_local, "x_{pom} B0 Reco", "pe");
            leg->AddEntry(h_xpom_def_RP_local, "x_{pom} RP Reco", "pe");
            leg->AddEntry(h_xpom_def_Sum, "B0+RP Sum", "pe");
            leg->Draw();

            c_xpom->SaveAs("figs/distributions/xpom_comparison_all_logxy.png");
            delete c_xpom;
        }
    }

    // =================================================================
    // 2D resolution circle maps
    // =================================================================
    auto createCirclePlot = [&](const char* histName, const char* saveName,
                                const char* xTitle, const char* yTitle,
                                bool logX = true, bool logY = true) {
        TProfile2D* prof = (TProfile2D*)inputFile->Get(histName);
        if (!prof) {
            std::cout << "Warning: Histogram " << histName << " not found, skipping..." << std::endl;
            return;
        }

        TCanvas* c = new TCanvas("c_circle", "Circle Plot", 1200, 900);
        c->SetLogx(logX);
        c->SetLogy(logY);
        c->SetGrid();
        c->SetRightMargin(0.15);

        TH2D* frame = new TH2D("frame", Form(";%s;%s", xTitle, yTitle),
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

        c->SaveAs(saveName);
        delete c;
        delete frame;
    };

    auto createBestMethodPlot = [&](const char* histName_EM, const char* histName_DA, const char* histName_Sigma,
                                     const char* saveName, const char* xTitle, const char* yTitle,
                                     bool logX = true, bool logY = true) {
        TProfile2D* prof_EM = (TProfile2D*)inputFile->Get(histName_EM);
        TProfile2D* prof_DA = (TProfile2D*)inputFile->Get(histName_DA);
        TProfile2D* prof_Sigma = (TProfile2D*)inputFile->Get(histName_Sigma);

        if (!prof_EM || !prof_DA || !prof_Sigma) {
            std::cout << "Warning: One or more histograms for best method plot not found, skipping..." << std::endl;
            return;
        }

        TCanvas* c = new TCanvas("c_best", "Best Method Plot", 1200, 900);
        c->SetLogx(logX);
        c->SetLogy(logY);
        c->SetGrid();
        c->SetRightMargin(0.15);

        TH2D* frame = new TH2D("frame_best", Form(";%s;%s", xTitle, yTitle),
                               prof_EM->GetNbinsX(), prof_EM->GetXaxis()->GetXmin(), prof_EM->GetXaxis()->GetXmax(),
                               prof_EM->GetNbinsY(), prof_EM->GetYaxis()->GetXmin(), prof_EM->GetYaxis()->GetXmax());
        frame->Draw();

        for (int ix = 1; ix <= prof_EM->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= prof_EM->GetNbinsY(); iy++) {
                double rms_EM = prof_EM->GetBinError(ix, iy);
                double rms_DA = prof_DA->GetBinError(ix, iy);
                double rms_Sigma = prof_Sigma->GetBinError(ix, iy);
                double entries = prof_EM->GetBinEntries(prof_EM->GetBin(ix, iy));

                if (entries < 10) continue;

                double best_rms = rms_EM;
                int best_method = 0;
                if (rms_DA > 0 && rms_DA < best_rms) { best_rms = rms_DA; best_method = 1; }
                if (rms_Sigma > 0 && rms_Sigma < best_rms) { best_rms = rms_Sigma; best_method = 2; }

                double x = prof_EM->GetXaxis()->GetBinCenter(ix);
                double y = prof_EM->GetYaxis()->GetBinCenter(iy);

                TMarker* marker = new TMarker(x, y, 20);
                marker->SetMarkerSize(1.5);
                if (best_method == 0) marker->SetMarkerColor(kRed);
                else if (best_method == 1) marker->SetMarkerColor(kBlue);
                else marker->SetMarkerColor(kGreen + 2);
                marker->Draw();
            }
        }

        TLegend* leg = new TLegend(0.7, 0.85, 0.9, 0.95);
        TMarker* m1 = new TMarker(0, 0, 20); m1->SetMarkerColor(kRed);
        TMarker* m2 = new TMarker(0, 0, 20); m2->SetMarkerColor(kBlue);
        TMarker* m3 = new TMarker(0, 0, 20); m3->SetMarkerColor(kGreen + 2);
        leg->AddEntry(m1, "EM Best", "p");
        leg->AddEntry(m2, "DA Best", "p");
        leg->AddEntry(m3, "Sigma Best", "p");
        leg->Draw();

        c->SaveAs(saveName);
        delete c;
        delete frame;
    };

    createCirclePlot("Q2_RelRes_vs_xy_EM", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_EM.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_DA", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_DA.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_Sigma", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_Sigma.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    // Additional circle maps for x_Bj and inelasticity
    createCirclePlot("x_RelRes_vs_xQ2_EM", "figs/resolutions/2d_maps/x_RelRes_xQ2_EM.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_DA", "figs/resolutions/2d_maps/x_RelRes_xQ2_DA.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_Sigma", "figs/resolutions/2d_maps/x_RelRes_xQ2_Sigma.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    createCirclePlot("y_RelRes_vs_xQ2_EM", "figs/resolutions/2d_maps/y_RelRes_xQ2_EM.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_DA", "figs/resolutions/2d_maps/y_RelRes_xQ2_DA.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_Sigma", "figs/resolutions/2d_maps/y_RelRes_xQ2_Sigma.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    createCirclePlot("t_RelRes_vs_xpomQ2_B0", "figs/resolutions/2d_maps/t_RelRes_xpomQ2_B0.png",
                     "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("t_RelRes_vs_xpomQ2_RP", "figs/resolutions/2d_maps/t_RelRes_xpomQ2_RP.png",
                     "x_{pom}", "Q^{2} [GeV^{2}]", true, true);

    createCirclePlot("xpom_RelRes_vs_xpomQ2_B0", "figs/resolutions/2d_maps/xpom_RelRes_xpomQ2_B0.png",
                     "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("xpom_RelRes_vs_xpomQ2_RP", "figs/resolutions/2d_maps/xpom_RelRes_xpomQ2_RP.png",
                     "x_{pom}", "Q^{2} [GeV^{2}]", true, true);

    createCirclePlot("beta_RelRes_vs_betaQ2_B0", "figs/resolutions/2d_maps/beta_RelRes_betaQ2_B0.png",
                     "#beta", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("beta_RelRes_vs_betaQ2_RP", "figs/resolutions/2d_maps/beta_RelRes_betaQ2_RP.png",
                     "#beta", "Q^{2} [GeV^{2}]", false, true);

    createCirclePlot("xL_RelRes_vs_xLQ2_B0", "figs/resolutions/2d_maps/xL_RelRes_xLQ2_B0.png",
                     "x_{L}", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("xL_RelRes_vs_xLQ2_RP", "figs/resolutions/2d_maps/xL_RelRes_xLQ2_RP.png",
                     "x_{L}", "Q^{2} [GeV^{2}]", false, true);

    createBestMethodPlot("Q2_RelRes_vs_xy_EM", "Q2_RelRes_vs_xy_DA", "Q2_RelRes_vs_xy_Sigma",
                         "figs/resolutions/2d_maps/Q2_RelRes_Q2x_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("x_RelRes_vs_xQ2_EM", "x_RelRes_vs_xQ2_DA", "x_RelRes_vs_xQ2_Sigma",
                         "figs/resolutions/2d_maps/x_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("y_RelRes_vs_xQ2_EM", "y_RelRes_vs_xQ2_DA", "y_RelRes_vs_xQ2_Sigma",
                         "figs/resolutions/2d_maps/y_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    // MX2 2D heatmap
    TProfile2D* prof_MX2 = (TProfile2D*)inputFile->Get("MX2_RelRes_vs_MX2Q2");
    if (prof_MX2) {
        TCanvas* c_MX2 = new TCanvas("c_MX2", "MX2 Resolution", 1200, 900);
        c_MX2->SetRightMargin(0.15);
        c_MX2->SetGrid();
        prof_MX2->SetTitle("M_{X}^{2} Resolution vs (M_{X}^{2}, Q^{2});M_{X}^{2} [GeV^{2}];Q^{2} [GeV^{2}]");
        prof_MX2->Draw("COLZ TEXT");
        c_MX2->SaveAs("figs/resolutions/2d_maps/MX2_RelRes_MX2Q2.png");
        delete c_MX2;
    }

    // E-pz 2D correlation
    TH2D* h_EPz_2D = (TH2D*)inputFile->Get("h_EPz_2D");
    if (h_EPz_2D) {
        TCanvas* c_EPz = new TCanvas("c_EPz_2D", "E-pz Correlation", 1200, 900);
        c_EPz->SetRightMargin(0.15);
        c_EPz->SetGrid();
        gStyle->SetPalette(kBird);
        h_EPz_2D->SetTitle("E-p_{z} Truth vs Reco;#Sigma(E-p_{z})_{truth} [GeV];#Sigma(E-p_{z})_{reco} [GeV]");
        h_EPz_2D->Draw("COLZ");
        c_EPz->SaveAs("figs/distributions/EPz_2D.png");
        delete c_EPz;
    }

    // Differential cross section with exponential fit
    TH1D* h_dsigma_dt_MC = (TH1D*)inputFile->Get("dsigma_dt_MC");
    TH1D* h_dsigma_dt_B0 = (TH1D*)inputFile->Get("dsigma_dt_B0");
    TH1D* h_dsigma_dt_RP = (TH1D*)inputFile->Get("dsigma_dt_RP");
    TH1D* h_dsigma_dt_Sum = (TH1D*)gDirectory->Get("dsigma_dt_Sum");

    TF1* fit_exp = nullptr;
    double b_value = 0.0;
    double b_error = 0.0;
    if (h_dsigma_dt_MC) {
        fit_exp = new TF1("fit_exp", "[0]*TMath::Exp(-[1]*x)", 0.01, 2.0);
        fit_exp->SetParameters(1000, 5.0);
        fit_exp->SetParNames("A", "b");
        fit_exp->SetLineColor(kMagenta + 2);
        fit_exp->SetLineWidth(3);
        fit_exp->SetLineStyle(2);
        h_dsigma_dt_MC->Fit(fit_exp, "RSQ");
        b_value = fit_exp->GetParameter(1);
        b_error = fit_exp->GetParError(1);
    }

    if (h_dsigma_dt_MC && h_dsigma_dt_B0 && h_dsigma_dt_RP && h_dsigma_dt_Sum && fit_exp) {
        TCanvas* c_dsigma_linear = new TCanvas("c_dsigma_dt", "Differential Cross Section", 1200, 900);
        c_dsigma_linear->SetLogx();
        c_dsigma_linear->SetGrid();

        h_dsigma_dt_MC->SetLineColor(kBlack);
        h_dsigma_dt_MC->SetLineWidth(2);
        h_dsigma_dt_MC->SetTitle("Differential Cross Section d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]");
        h_dsigma_dt_MC->Draw("HIST");

        h_dsigma_dt_B0->SetMarkerStyle(20);
        h_dsigma_dt_B0->SetMarkerColor(kRed);
        h_dsigma_dt_B0->SetLineColor(kRed);
        h_dsigma_dt_B0->Draw("PE SAME");

        h_dsigma_dt_RP->SetMarkerStyle(20);
        h_dsigma_dt_RP->SetMarkerColor(kBlue);
        h_dsigma_dt_RP->SetLineColor(kBlue);
        h_dsigma_dt_RP->Draw("PE SAME");

        h_dsigma_dt_Sum->SetMarkerStyle(20);
        h_dsigma_dt_Sum->SetMarkerColor(kOrange + 7);
        h_dsigma_dt_Sum->SetLineColor(kOrange + 7);
        h_dsigma_dt_Sum->Draw("PE SAME");

        fit_exp->Draw("SAME");

        TLegend* leg1 = new TLegend(0.6, 0.55, 0.85, 0.9);
        leg1->AddEntry(h_dsigma_dt_MC, "MC Truth", "l");
        leg1->AddEntry(h_dsigma_dt_B0, "B0 Reco", "pe");
        leg1->AddEntry(h_dsigma_dt_RP, "RP Reco", "pe");
        leg1->AddEntry(h_dsigma_dt_Sum, "B0+RP Sum", "pe");
        leg1->AddEntry(fit_exp, Form("Fit: e^{-b|t|}, b = %.2f #pm %.2f GeV^{-2}", b_value, b_error), "l");
        leg1->Draw();

        c_dsigma_linear->SaveAs("figs/cross_sections/dsigma_dt_with_fit.png");
        delete c_dsigma_linear;

        TCanvas* c_dsigma_logy = new TCanvas("c_dsigma_dt_logy", "Differential Cross Section", 1200, 900);
        c_dsigma_logy->SetLogx();
        c_dsigma_logy->SetLogy();
        c_dsigma_logy->SetGrid();

        h_dsigma_dt_MC->Draw("HIST");
        h_dsigma_dt_B0->Draw("PE SAME");
        h_dsigma_dt_RP->Draw("PE SAME");
        h_dsigma_dt_Sum->Draw("PE SAME");
        fit_exp->Draw("SAME");

        TLegend* leg2 = new TLegend(0.6, 0.55, 0.85, 0.9);
        leg2->AddEntry(h_dsigma_dt_MC, "MC Truth", "l");
        leg2->AddEntry(h_dsigma_dt_B0, "B0 Reco", "pe");
        leg2->AddEntry(h_dsigma_dt_RP, "RP Reco", "pe");
        leg2->AddEntry(h_dsigma_dt_Sum, "B0+RP Sum", "pe");
        leg2->AddEntry(fit_exp, Form("Fit: e^{-b|t|}, b = %.2f #pm %.2f GeV^{-2}", b_value, b_error), "l");
        leg2->Draw();

        c_dsigma_logy->SaveAs("figs/cross_sections/dsigma_dt_logy_with_fit.png");
        delete c_dsigma_logy;
    }

    // Triple differential cross section slices
    TH3D* h_d3sigma_MC = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_MC");
    TH3D* h_d3sigma_B0 = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_B0");
    TH3D* h_d3sigma_RP = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_RP");

    if (h_d3sigma_MC && h_d3sigma_B0 && h_d3sigma_RP) {
        int n_Q2_bins = h_d3sigma_MC->GetNbinsX();
        int n_beta_bins = h_d3sigma_MC->GetNbinsY();
        int n_xpom_bins = h_d3sigma_MC->GetNbinsZ();

        TCanvas* c_beta = new TCanvas("c_d3sigma_beta", "d3sigma vs beta", 2400, 1600);
        c_beta->Divide(n_Q2_bins, n_xpom_bins);

        for (int i_Q2 = 1; i_Q2 <= n_Q2_bins; i_Q2++) {
            for (int i_xpom = 1; i_xpom <= n_xpom_bins; i_xpom++) {
                int pad_num = (i_xpom - 1) * n_Q2_bins + i_Q2;
                c_beta->cd(pad_num);
                gPad->SetLogy();

                TH1D* proj_MC = h_d3sigma_MC->ProjectionY(Form("proj_beta_MC_%d_%d", i_Q2, i_xpom), i_Q2, i_Q2, i_xpom, i_xpom);
                TH1D* proj_B0 = h_d3sigma_B0->ProjectionY(Form("proj_beta_B0_%d_%d", i_Q2, i_xpom), i_Q2, i_Q2, i_xpom, i_xpom);
                TH1D* proj_RP = h_d3sigma_RP->ProjectionY(Form("proj_beta_RP_%d_%d", i_Q2, i_xpom), i_Q2, i_Q2, i_xpom, i_xpom);

                proj_MC->SetLineColor(kBlack);
                proj_B0->SetMarkerColor(kRed);
                proj_B0->SetMarkerStyle(20);
                proj_RP->SetMarkerColor(kBlue);
                proj_RP->SetMarkerStyle(20);

                proj_MC->Draw("HIST");
                proj_B0->Draw("PE SAME");
                proj_RP->Draw("PE SAME");
            }
        }

        c_beta->SaveAs("figs/cross_sections/d3sigma_vs_beta.png");
        delete c_beta;

        TCanvas* c_xpom = new TCanvas("c_d3sigma_xpom", "d3sigma vs xpom", 2400, 1000);
        c_xpom->Divide(n_Q2_bins, n_beta_bins);

        for (int i_Q2 = 1; i_Q2 <= n_Q2_bins; i_Q2++) {
            for (int i_beta = 1; i_beta <= n_beta_bins; i_beta++) {
                int pad_num = (i_beta - 1) * n_Q2_bins + i_Q2;
                c_xpom->cd(pad_num);
                gPad->SetLogx();
                gPad->SetLogy();

                TH1D* proj_MC = h_d3sigma_MC->ProjectionZ(Form("proj_xpom_MC_%d_%d", i_Q2, i_beta), i_Q2, i_Q2, i_beta, i_beta);
                TH1D* proj_B0 = h_d3sigma_B0->ProjectionZ(Form("proj_xpom_B0_%d_%d", i_Q2, i_beta), i_Q2, i_Q2, i_beta, i_beta);
                TH1D* proj_RP = h_d3sigma_RP->ProjectionZ(Form("proj_xpom_RP_%d_%d", i_Q2, i_beta), i_Q2, i_Q2, i_beta, i_beta);

                proj_MC->SetLineColor(kBlack);
                proj_B0->SetMarkerColor(kRed);
                proj_B0->SetMarkerStyle(20);
                proj_RP->SetMarkerColor(kBlue);
                proj_RP->SetMarkerStyle(20);

                proj_MC->Draw("HIST");
                proj_B0->Draw("PE SAME");
                proj_RP->Draw("PE SAME");
            }
        }

        c_xpom->SaveAs("figs/cross_sections/d3sigma_vs_xpom.png");
        delete c_xpom;
    }

    inputFile->Close();
    delete inputFile;

    if (!binningInput.empty()) {
        RunBinningPlots(binningInput, binningScheme, binningOutput);
    }

    return 0;
}
