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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "BinningUtility.hpp"

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
        std::cerr << "Usage: " << argv[0] << " <DDIS_BinningInputs.root> [binning_scheme.txt] [output.root]" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string binningFile = (argc > 2) ? argv[2] : "binning_scheme.txt";
    std::string outputFile = (argc > 3) ? argv[3] : "DDIS_BinningPlots.root";

    BinningScheme binning = ReadBinningScheme(binningFile);
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
    CombinedHists h_all = BookCombinedHists(binning);

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

        if(Q2_truth > 0 && beta_truth > 0 && xpom_truth > 0) {
            iQ2_true = binning.GetQ2Bin(Q2_truth);
            iBeta_true = binning.GetBetaBin(beta_truth);
            iXpom_true = binning.GetXpomBin(xpom_truth);
        }
        if(Q2_reco > 0 && beta_reco > 0 && xpom_reco > 0) {
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

        if(k_true >= 0) {
            h_all.truth_count[k_true] += 1.0;
        }
        if(k_reco >= 0) {
            h_all.reco_total[k_reco] += 1.0;
        }
        if(k_true >= 0 && k_reco >= 0) {
            h_all.truth_with_reco[k_true] += 1.0;
            h_all.h_migration->Fill(k_true + 1, k_reco + 1);
            if(k_true == k_reco) {
                h_all.diag_count[k_true] += 1.0;
            }
        }

        if(k_true >= 0 && k_reco >= 0) {
            if(Q2_truth != 0.0f) {
                double res_Q2 = (Q2_reco - Q2_truth) / Q2_truth;
                h_all.h_res_Q2->Fill(Q2_truth, res_Q2);
            }
            if(beta_truth != 0.0f) {
                double res_beta = (beta_reco - beta_truth) / beta_truth;
                h->h_res_beta->Fill(beta_truth, res_beta);

                if(xpom_truth != 0.0f) {
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
    std::cout << "Wrote plots to " << outdir << " and histograms to " << outputFile << std::endl;

    return 0;
}
