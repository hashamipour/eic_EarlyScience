// Inclusive DIS plotter
// g++ DDIS_Plot_Final.cpp -o DDIS_Plot_Final $(root-config --cflags --glibs)
// ./DDIS_Plot_Final <combined.root>

#include "Plotting.hpp"

#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TTree.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

static void PlotXQ2Density(TFile* inputFile) {
    if (!inputFile) return;
    TTree* tree = (TTree*)inputFile->Get("Q2_tree");
    if (!tree) {
        std::cerr << "Warning: Q2_tree not found; skipping x-Q2 density plot." << std::endl;
        return;
    }

    if (!tree->GetBranch("x_truth") || !tree->GetBranch("Q2_truth")) {
        TH2D* h_existing = (TH2D*)inputFile->Get("xQ2_truth");
        if (!h_existing) {
            std::cerr << "Warning: Missing branches x_truth/Q2_truth and no xQ2_truth histogram found; "
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

        std::cerr << "Note: Q2_tree lacks Q2_truth; using stored xQ2_truth histogram from skim output. "
                     "Re-skim with Q2_truth in the tree for finer/unbinned x-Q2 plotting." << std::endl;
        return;
    }

    float x_truth = -999.0f;
    float Q2_truth = -999.0f;
    tree->SetBranchAddress("x_truth", &x_truth);
    tree->SetBranchAddress("Q2_truth", &Q2_truth);

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
        if (std::isfinite(x_truth) && x_truth > 0.0f && x_truth < 1.0f &&
            std::isfinite(Q2_truth) && Q2_truth > 0.0f) {
            h->Fill(x_truth, Q2_truth);
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
        std::cerr << "Warning: Histogram " << histName << " not found; skipping density plot." << std::endl;
        return;
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
    gSystem->mkdir("figs/diffractive", kTRUE);
    gSystem->mkdir("figs/diffractive/histos", kTRUE);
    gSystem->mkdir("figs/diffractive/resolution", kTRUE);
    gSystem->mkdir("figs/diffractive/resolution/profile", kTRUE);
    gSystem->mkdir("figs/diffractive/response", kTRUE);

    // Acceptance/purity plots (if tracking histograms exist)
    TH1D* h_gen_Q2 = (TH1D*)inputFile->Get("h_gen_Q2");
    TH1D* h_gen_and_reco_after_cuts_Q2_EM = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_EM");
    TH1D* h_gen_and_reco_after_cuts_Q2_DA = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_DA");
    TH1D* h_gen_and_reco_after_cuts_Q2_Sigma = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_Sigma");
    TH1D* h_reco_Q2_EM = (TH1D*)inputFile->Get("h_reco_Q2_EM");
    TH1D* h_reco_Q2_DA = (TH1D*)inputFile->Get("h_reco_Q2_DA");
    TH1D* h_reco_Q2_Sigma = (TH1D*)inputFile->Get("h_reco_Q2_Sigma");

    PlotXQ2Density(inputFile);
    PlotDensityFromHist(inputFile, "beta_Q2_truth", "#beta", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/beta_q2_density.png", false, true);
    PlotDensityFromHist(inputFile, "t_Q2_truth", "|t| [GeV^{2}]", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/t_q2_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_Q2_truth", "x_{pom}", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/histos/xpom_q2_density.png", true, true);
    PlotDensityFromHist(inputFile, "beta_t_truth", "#beta", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/beta_t_density.png", false, true);
    PlotDensityFromHist(inputFile, "xbj_t_truth", "x_{Bj}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xbj_t_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_t_truth", "x_{pom}", "|t| [GeV^{2}]",
                        "figs/diffractive/histos/xpom_t_density.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_beta_truth", "x_{pom}", "#beta",
                        "figs/diffractive/histos/xpom_beta_density.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_beta_truth", "x_{Bj}", "#beta",
                        "figs/diffractive/histos/xbj_beta_density.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_xpom_truth", "x_{Bj}", "x_{pom}",
                        "figs/diffractive/histos/xbj_xpom_density.png", true, true);

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
