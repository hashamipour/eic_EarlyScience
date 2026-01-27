//g++ -o DDIS_Plot_Combined $(root-config --cflags --glibs) DDIS_Plot_Combined.cpp Plotting.cpp
// FINAL COMPREHENSIVE COMBINED PLOTTING SCRIPT
// Merges ALL plotting logic from:
//   - DDIS_Plot_Combined.cpp (primary comprehensive file)
//   - DDIS_Plots_Q2_with_xy.cpp (Q2/xy analysis with commented plots)
//   - DDIS_Plot_t.cpp (Mandelstam t analysis with commented plots)
//
// This version includes the UNION of all plots (active + commented) from all three input files
// Preserves output directory structure exactly as in individual scripts

#include "Plotting.hpp"
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include "TError.h"
#include <TCanvas.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TProfile2D.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TString.h>
#include <TF1.h>
#include <TLegend.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    gErrorIgnoreLevel = kWarning; // Suppress ROOT infos

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

    std::cout << "\n========================================" << std::endl;
    std::cout << "COMPREHENSIVE COMBINED PLOTTING SCRIPT" << std::endl;
    std::cout << "========================================\n" << std::endl;

    //=================================================================
    // SECTION 1: Q2/XY ANALYSIS PLOTS
    // Source: DDIS_Plot_Combined.cpp + DDIS_Plots_Q2_with_xy.cpp
    //=================================================================

    std::cout << "Adding Q2/xy analysis plots..." << std::endl;

    // Q2 1D histogram plots (counts)
    plot_ptr = new PlotOptions1D(
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe","pe"},
        "Q^{2} Reconstruction Methods",
        "Q^{2}",
        "# of events",
        "figs/distributions/Q2_hist.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // Q2 PDF comparison
    plot_ptr = new PlotOptions1D(
        {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA","h_Q2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe","pe"},
        "Q^{2} PDF Comparison",
        "Q^{2}",
        "PDF",
        "figs/distributions/Q2_pdf.png",
        true,  // logX
        true,  // logY
        true   // normalize to PDF
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // x_Bj distributions (counts)
    plot_ptr = new PlotOptions1D(
        {"x_truth", "x_EM", "x_DA", "x_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{Bj} Reconstruction Methods",
        "x_{Bj}",
        "# of events",
        "figs/distributions/x_hist.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // x_Bj PDF comparison
    plot_ptr = new PlotOptions1D(
        {"x_truth", "x_EM", "x_DA", "x_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{Bj} PDF Comparison",
        "x_{Bj}",
        "PDF",
        "figs/distributions/x_pdf.png",
        true,  // logX
        true,  // logY
        true   // normalize to PDF
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // y (inelasticity) distributions (counts)
    plot_ptr = new PlotOptions1D(
        {"y_truth", "y_EM", "y_DA", "y_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "y (inelasticity) Reconstruction Methods",
        "y",
        "# of events",
        "figs/distributions/y_hist.png",
        false,  // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // y PDF comparison
    plot_ptr = new PlotOptions1D(
        {"y_truth", "y_EM", "y_DA", "y_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe"},
        "y (inelasticity) PDF Comparison",
        "y",
        "PDF",
        "figs/distributions/y_pdf.png",
        false, // logX
        false, // logY
        true   // normalize to PDF
    );
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    //=================================================================
    // SECTION 1A: Q2/X/Y OVERALL RELATIVE RESOLUTION (from Q2_with_xy)
    // These were commented in the original but should be included
    //=================================================================

    std::cout << "Adding Q2/x/y overall resolution plots..." << std::endl;

    // Q2 Relative Resolution plots (PlotOptionsRelRes)
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

    // x_Bj Relative Resolution plots
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
        0., 0.,  // Skip fitting, just save histogram
        "figs/resolutions/simple/DDIS_RelRes_xBj_DA.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_Sigma",
        "#frac{x_{#Sigma} - x_{MC}}{X_{MC}}",
        "Counts",
        0., 0.,  // Skip fitting
        "figs/resolutions/simple/DDIS_RelRes_x_Sigma.png"
    ));

    // y Relative Resolution plots
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
        0., 0.,  // Skip fitting
        "figs/resolutions/simple/DDIS_RelRes_y_DA.png"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_Sigma",
        "#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "Counts",
        0., 0.,  // Skip fitting
        "figs/resolutions/simple/DDIS_RelRes_y_Sigma.png"
    ));

    //=================================================================
    // SECTION 1B: Q2/X/Y BINNED RELATIVE RESOLUTION (from Q2_with_xy)
    // These were commented in the original but included for completeness
    //=================================================================

    std::cout << "Adding Q2/x/y binned resolution plots..." << std::endl;

    // Q2 Binned Resolution (EM method)
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_EM",
        "Relative bin by bin resolution (EM);Q^{2}_{MC};#frac{Q^{2}_{EM} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{EM}",
        "",
        {
         {-0.0, 0.0},/*2*/ {-0.022, 0.02},{-0.02, 0.02},{-0.02, 0.02},{-0.02, 0.02},
         {-0.015, 0.015},/*7*/{-0.015, 0.015},/*8*/{-0.014, 0.015},{-0.025, 0.025},{-0.01, 0.012},
         {-0.027, 0.028},{-0.018, 0.02},{-0.022, 0.02},{-0.02, 0.015},{-0.018, 0.02},
         {-0.02, 0.015},{-0.02, 0.017},{-0.017, 0.02},{-0.02, 0.02},{-0.04, 0.04},
         {-0.025, 0.03},{-0.015, 0.025},{-0.05, 0.06}
        },
        "figs/resolutions/binned/DDIS_Q2RelRes_binned_EM.png",
        "DDIS_Q2RelRes_binned_EM",
        std::make_pair(5.0, 200),
        true
    ));

    // Q2 Binned Resolution (DA method)
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_DA",
        "Relative bin by bin resolution (DA);Q^{2}_{MC};#frac{Q^{2}_{DA} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{DA}",
        "",
        {
          {-0.0, 0.0}, {-0.005, 0.025},{-0.01, 0.025},/*4*/{-0.006, 0.02},{-0.01, 0.02},
          {-0.009, 0.02},/*7*/{-0, 0},{-0., 0.},{-0., 0.},{-0., 0.},
          {-0.015, 0.03},{-0.009, 0.02},{-0.01, 0.02},/*14*/{-0.01, 0.02},{-0.01, 0.02},
          {-0.01, 0.02},/*17*/{-0.004, 0.02},{-0.017, 0.027},{-0.025, 0.03},{-0.08, 0.08},
          {-0.05, 0.06},{-0.05, 0.065},{-0.05, 0.06}
        },
        "figs/resolutions/binned/DDIS_Q2RelRes_binned_DA.png",
        "DDIS_Q2RelRes_binned_DA",
        std::make_pair(5.0, 200),
        true
    ));

    // Q2 Binned Resolution (Sigma method)
    plots.push_back(new PlotOptionsBinnedRelRes(
        "Q2_RelRes_binned_Sigma",
        ";Q^{2}_{MC};#frac{Q^{2}_{#Sigma} - Q^{2}_{MC}}{Q^{2}_{MC}}",
        "Q^{2}_{#Sigma}",
        "",
        {
         {-0.0, 0.0},/*2*/ {-0.022, 0.02},{-0.02, 0.02},{-0.02, 0.02},{-0.02, 0.02},
         {-0.015, 0.015},/*7*/{-0.015, 0.015},/*8*/{-0.014, 0.015},{-0.025, 0.025},{-0.01, 0.012},
         {-0.027, 0.028},{-0.018, 0.02},{-0.022, 0.02},{-0.02, 0.015},{-0.018, 0.02},
         {-0.02, 0.015},{-0.02, 0.017},{-0.017, 0.02},{-0.02, 0.02},{-0.04, 0.04},
         {-0.025, 0.03},{-0.015, 0.025},{-0.05, 0.06}
        },
        "figs/resolutions/binned/DDIS_Q2RelRes_binned_Sigma.png",
        "DDIS_Q2RelRes_binned_Sigma",
        std::make_pair(5.0, 200),
        true
    ));

    // x_Bj Binned Resolution (EM method)
    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "x_RelRes_binned_EM",
        ";x_{MC};#frac{x_{EM} - x_{MC}}{x_{MC}}",
        "x_{EM}",
        "",
        {
         {-0.0, 0.0},/*2*/ {-0.022, 0.02},{-0.02, 0.02},/*4*/{-0.02, 0.02},{-0.02, 0.02},
         {-0.015, 0.015},/*7*/{-0.015, 0.015},/*8*/{0,0},{-0.025, 0.025},{-0.025, 0.02},
         {-0.027, 0.028},/*12*/{0, 0},{0, 0},{0, 0},{0, 0},
         {0, 0},{0, 0},{0, 0},{0, 0},{0, 0},
         {0, 0},{0, 0},{0, 0}
        },
        "figs/resolutions/binned/DDIS_RelRes_binned_x_EM.png",
        "DDIS_RelRes_binned_x_EM",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    plots.push_back(binned_plot_ptr);

    // x_Bj Binned Resolution (DA method)
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

    // x_Bj Binned Resolution (Sigma method)
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

    // y Binned Resolution (EM method)
    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_EM",
        ";y_{MC};#frac{y_{EM} - y_{MC}}{y_{MC}}",
        "y_{EM}",
        "",
        {
            {-0.,0.},{-0.,0.},{-0.,0.},/*4*/{-0.08,0.08},
            {-0.06,0.07},/*6*/{-0.04,0.045},{-0.025,0.025},/*8*/{-0.03,0.035},
            {-0.03,0.036},/*10*/{-0.02,0.02},{-0.018,0.018},/*12*/{-0.015,0.015},
            {-0.015,0.015},/*14*/{-0.01,0.015},{-0.01,0.01},/*16*/{-0.011,0.011}
        },
        "figs/resolutions/binned/DDIS_RelRes_binned_y_EM.png",
        "DDIS_RelRes_binned_y_EM",
        std::make_pair(0.0, 1.0)
    );
    binned_plot_ptr->SetLegendPosition(0.75, 0.75, 0.9, 0.9);
    plots.push_back(binned_plot_ptr);

    // y Binned Resolution (DA method)
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

    // y Binned Resolution (Sigma method)
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

    //=================================================================
    // SECTION 1C: Q2/X/Y RESPONSE MATRICES (Correlation Plots)
    //=================================================================

    std::cout << "Adding Q2/x/y response matrices..." << std::endl;

    // Q2 Response Matrices
    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_EM",
        "Q^{2} (true) [GeV]",
        "Q^{2} (EM) [GeV]",
        "figs/response_matrices/response_matrix_Q2_EM.png",
        true,  // isLogX
        true,  // isLogY
        std::make_pair(1.0, 300),
        std::make_pair(1.0, 300)
    ));

    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_DA",
        "Q^{2} (true) [GeV]",
        "Q^{2} (DA) [GeV]",
        "figs/response_matrices/response_matrix_Q2_DA.png",
        true,  // isLogX
        true,  // isLogY
        std::make_pair(1.0, 300),
        std::make_pair(1.0, 300)
    ));

    plots.push_back(new PlotOptionsResponseMatrix (
        "Corr_Q2_Sigma",
        "Q^{2} (true) [GeV]",
        "Q^{2} (Sigma) [GeV]",
        "figs/response_matrices/response_matrix_Q2_Esigma.png",
        true,  // isLogX
        true,  // isLogY
        std::make_pair(1.0, 300),
        std::make_pair(1.0, 300)
    ));

    // x_Bj Response Matrices
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_EM",
        "x_{Bj} (true)",
        "x_{Bj} (EM)",
        "figs/response_matrices/response_matrix_x_EM.png",
        true,  // isLogX
        true,   // isLogY
        {1e-3,0.3},
        {1e-3,0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_DA",
        "x_{Bj} (true)",
        "x_{Bj} (DA)",
        "figs/response_matrices/response_matrix_x_DA.png",
        true,  // isLogX
        true,   // isLogY
        {1e-3,0.3},
        {1e-3,0.3}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_Sigma",
        "x_{Bj} (true)",
        "x_{Bj} (Sigma)",
        "figs/response_matrices/response_matrix_x_Sigma.png",
        true,  // isLogX
        true,   // isLogY
        {1e-3,0.3},
        {1e-3,0.5}
    ));

    // y Response Matrices
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_EM",
        "y (true)",
        "y (EM)",
        "figs/response_matrices/response_matrix_y_EM.png",
        false, // isLogX
        false,  // isLogY
        {0.,1.},
        {0.,1.}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_DA",
        "y (true)",
        "y (DA)",
        "figs/response_matrices/response_matrix_y_DA.png",
        false, // isLogX
        false,  // isLogY
        {0.,1.},
        {0.,1.}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_Sigma",
        "y (true)",
        "y (Sigma)",
        "figs/response_matrices/response_matrix_y_Sigma.png",
        false, // isLogX
        false,  // isLogY
        {0.,1.},
        {0.,1.}
    ));

    //=================================================================
    // SECTION 1D: HADRONIC FINAL STATE VARIABLES (E-pz, eta_max, MX2)
    //=================================================================

    std::cout << "Adding hadronic final state variable plots..." << std::endl;

    // E-pz distribution plots
    plot_ptr = new PlotOptions1D(
        {"h_EPz_truth", "h_EPz"},
        {"MC Truth"   , "Reconstruction"},
        {"hist"       , "E1"},
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
        {"MC Truth"   , "Reconstruction"},
        {"hist"       , "E1"},
        "Hadronic Final State E-p_{z}",
        "#Sigma(E-p_{z}) [GeV]",
        "Counts",
        "figs/distributions/EPz_distribution_logY.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.2, 0.7, 0.4, 0.9);
    plots.push_back(plot_ptr);

    // eta_max distribution plots
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

    // M_X^2 distribution plots
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

    //=================================================================
    // SECTION 2: MANDELSTAM t ANALYSIS PLOTS
    // Source: DDIS_Plot_Combined.cpp + DDIS_Plot_t.cpp
    //=================================================================

    std::cout << "Adding Mandelstam t analysis plots..." << std::endl;

    // |t| distributions
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

    // |t| distributions with log y
    plot_ptr = new PlotOptions1D(
        {"t_MC", "t_B0", "t_RP_histo"},
        {"MC Truth", "Reco B0 (5.5-20 mrad)", "Reco RP (<5 mrad)"},
        {"hist", "pe", "pe"},
        "Mandelstam |t| Distributions",
        "|t| [GeV^{2}]",
        "Counts",
        "figs/distributions/t_distributions_logy.png",
        false,
        true // logY
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    // |t| PDF comparison
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
        true // normalize to PDF
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    // Angular distributions
    plot_ptr = new PlotOptions1D(
        {"theta_MC", "theta_B0", "theta_RP"},
        {"MC Truth", "Reco B0", "Reco RP"},
        {"hist", "pe", "pe"},
        "Proton Scattering Angles",
        "#theta [mrad]",
        "Counts",
        "figs/distributions/theta_distributions.png",
        false,
        true  // Log y for better visibility
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    //=================================================================
    // SECTION 2A: DIFFRACTIVE VARIABLE OVERALL RESOLUTION PLOTS
    //=================================================================

    std::cout << "Adding diffractive variable resolution plots..." << std::endl;

    // |t| overall resolution plots (B0 and RP)
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

    // x_L overall resolution plots (B0 and RP)
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

    // x_pom overall resolution plots (B0 and RP)
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

    //=================================================================
    // SECTION 2B: DIFFRACTIVE VARIABLE BINNED RESOLUTION PLOTS
    //=================================================================

    std::cout << "Adding diffractive variable binned resolution plots..." << std::endl;

    // |t| binned resolution (B0 and RP)
    PlotOptionsBinnedRelRes* plot_t_B0 = new PlotOptionsBinnedRelRes(
        "t_RelRes_binned_B0",
        "B0 |t| Resolution vs Truth |t|",
        "|t|_{truth} [GeV^{2}]",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        {{-0.3, 0.3}}, // Single fit range for all bins
        "figs/resolutions/binned/t_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/t_B0",
        {0.0, 2.0},
        false // linear x-axis for |t|
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

    // x_L binned resolution (B0 and RP)
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

    // x_pom binned resolution (B0 and RP)
    PlotOptionsBinnedRelRes* plot_xpom_B0 = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_B0",
        "B0 x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{-0.5, 0.5}},
        "figs/resolutions/binned/xpom_resolution_binned_B0.png",
        "figs/resolutions/binned/bins/xpom_B0",
        {1e-4, 0.4},
        true // log x-axis for x_pom
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

    // beta binned resolution (B0 and RP)
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

    // M_X^2 binned resolution
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

    //=================================================================
    // SECTION 2C: DIFFRACTIVE VARIABLE RESPONSE MATRICES
    //=================================================================

    std::cout << "Adding diffractive variable response matrices..." << std::endl;

    // Response matrices for t
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

    // x_L response matrices
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

    // x_pom correlation matrices (from x_L method)
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0",
        "Truth x_{pom} (1-x_L)",
        "B0 Reco x_{pom} (1-x_L)",
        "figs/response_matrices/response_matrix_xpom_B0.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_RP",
        "Truth x_{pom} (1-x_L)",
        "RP Reco x_{pom} (1-x_L)",
        "figs/response_matrices/response_matrix_xpom_RP.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    //=================================================================
    // SECTION 3: X_POM COMPARISON PLOTS (Definition vs from x_L)
    //=================================================================

    std::cout << "Adding x_pom comparison plots..." << std::endl;

    // x_pom comparison: 1D overlays
    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_def_MC"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "MC Truth x_{pom} Comparison",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_MC_logxy.png",
        true,   // isLogX
        true    // isLogY
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
        true,   // isLogX
        true    // isLogY
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
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_def_MC", "xpom_def_B0", "xpom_def_RP", "xpom_def_Sum"},
        {"x_{pom} MC", "x_{pom} B0 Reco", "x_{pom} RP Reco", "B0+RP Sum"},
        {"hist", "E1", "E1", "E1"},
        "x_{pom} Comparison: All Methods",
        "x_{pom}",
        "Counts",
        "figs/distributions/xpom_comparison_all_logxy.png",
        true,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.5, 0.9);
    plots.push_back(plot_ptr);

    // x_pom comparison: 2D correlation plots
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_MC",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_MC.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_B0",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_B0.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_RP",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}-t)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/distributions/xpom_2D_comparison_RP.png",
        true,   // logX
        true,   // logY
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));

    //=================================================================
    // SECTION 4: PROTON THETA ANGLE ACCEPTANCE
    //=================================================================

    std::cout << "Adding proton theta acceptance plots..." << std::endl;

    // Theta angle comparison: all truth-seeded protons vs B0 accepted
    plot_ptr = new PlotOptions1D(
        {"theta_all_TS", "theta_B0"},
        {"All Truth-Seeded Protons", "B0 Accepted Protons"},
        {"hist", "hist"},
        "Proton Scattering Angle Distribution",
        "#theta [mrad]",
        "Counts",
        "figs/distributions/theta_comparison_B0_acceptance.png",
        false,   // isLogX
        false   // isLogY
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
        false,   // isLogX
        true    // isLogY
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    //=================================================================
    // SECTION 5: BETA ANALYSIS PLOTS (beta = x_Bj / x_pom)
    //=================================================================

    std::cout << "Adding beta analysis plots..." << std::endl;

    // Beta distributions with log y-axis for better visibility
    PlotOptions1D* plot_beta_logy = new PlotOptions1D(
        {"beta_MC", "beta_B0", "beta_RP", "beta_Sum"},
        {"MC Truth", "B0 Reco", "RP Reco", "B0+RP Sum"},
        {"hist", "E1", "E1", "E1"},
        "#beta = x_{Bj} / x_{pom} Distributions (from x_{pom} definition)",
        "#beta",
        "Counts",
        "figs/distributions/beta_distributions_logy.png",
        false,  // isLogX
        true    // isLogY
    );
    plot_beta_logy->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plot_beta_logy->SetRangeY(1, 4e4);
    plots.push_back(plot_beta_logy);

    // Beta resolution plots
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

    // Beta response matrices
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_B0",
        "Truth #beta",
        "B0 Reco #beta",
        "figs/response_matrices/response_matrix_beta_B0.png",
        false,  // isLogX
        false,  // isLogY
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_RP",
        "Truth #beta",
        "RP Reco #beta",
        "figs/response_matrices/response_matrix_beta_RP.png",
        false,  // isLogX
        false,  // isLogY
        {0.0, 1.0},
        {0.0, 1.0}
    ));

    // Beta physics correlations
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_Q2",
        "Q^{2} [GeV^{2}]",
        "#beta",
        "figs/distributions/beta_vs_Q2.png",
        true,   // logX (Q² often on log scale)
        false,  // logY
        {0.1, 100.0},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_xpom",
        "x_{pom}",
        "#beta",
        "figs/distributions/beta_vs_xpom.png",
        true,   // logX
        false,  // logY
        {1e-4, 0.4},
        {0.0, 1.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_t",
        "|t| [GeV^{2}]",
        "#beta",
        "figs/distributions/beta_vs_t.png",
        false,  // logX
        false,  // logY
        {0.0, 2.0},
        {0.0, 1.0}
    ));

    //=================================================================
    // SECTION 6: 2D RESOLUTION CIRCLE MAPS (TProfile2D visualization)
    //=================================================================

    std::cout << "Creating 2D resolution circle maps..." << std::endl;

    // Create output directories
    gSystem->mkdir("figs", kTRUE);
    gSystem->mkdir("figs/resolutions/2d_maps", kTRUE);

    // Helper function to create circle plots from TProfile2D
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

        // Create empty histogram for frame
        TH2D* frame = new TH2D("frame", Form(";%s;%s", xTitle, yTitle),
                                prof->GetNbinsX(), prof->GetXaxis()->GetXmin(), prof->GetXaxis()->GetXmax(),
                                prof->GetNbinsY(), prof->GetYaxis()->GetXmin(), prof->GetYaxis()->GetXmax());
        frame->Draw();

        // Draw circles with size proportional to resolution
        for (int ix = 1; ix <= prof->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= prof->GetNbinsY(); iy++) {
                double entries = prof->GetBinEntries(prof->GetBin(ix, iy));
                if (entries < 10) continue;  // Skip bins with low statistics

                double x = prof->GetXaxis()->GetBinCenter(ix);
                double y = prof->GetYaxis()->GetBinCenter(iy);
                double rms = prof->GetBinError(ix, iy);  // RMS of relative resolution

                // Scale marker size based on RMS (resolution)
                double markerSize = TMath::Min(rms * 1000.0, 3.0);  // Scale factor and cap

                TMarker* marker = new TMarker(x, y, 20);
                marker->SetMarkerSize(markerSize);

                // Hollow if entries < 100, filled if >= 100
                if (entries < 100) {
                    marker->SetMarkerStyle(24);  // Hollow circle
                    marker->SetMarkerColor(kBlue);
                } else {
                    marker->SetMarkerStyle(20);  // Filled circle
                    marker->SetMarkerColor(kRed);
                }
                marker->Draw();
            }
        }

        c->SaveAs(saveName);
        delete c;
        delete frame;
    };

    // Q2 Resolution Circle Maps
    createCirclePlot("Q2_RelRes_vs_xy_EM", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_EM.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_DA", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_DA.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_Sigma", "figs/resolutions/2d_maps/Q2_RelRes_Q2x_Sigma.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    // x Resolution Circle Maps
    createCirclePlot("x_RelRes_vs_xQ2_EM", "figs/resolutions/2d_maps/x_RelRes_xQ2_EM.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_DA", "figs/resolutions/2d_maps/x_RelRes_xQ2_DA.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_Sigma", "figs/resolutions/2d_maps/x_RelRes_xQ2_Sigma.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    // y Resolution Circle Maps
    createCirclePlot("y_RelRes_vs_xQ2_EM", "figs/resolutions/2d_maps/y_RelRes_xQ2_EM.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_DA", "figs/resolutions/2d_maps/y_RelRes_xQ2_DA.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_Sigma", "figs/resolutions/2d_maps/y_RelRes_xQ2_Sigma.png",
                     "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    // Diffractive Variable Resolution Circle Maps
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

    // Best Method Comparison Circle Plots
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

        TH2D* frame = new TH2D("frame", Form(";%s;%s", xTitle, yTitle),
                                prof_EM->GetNbinsX(), prof_EM->GetXaxis()->GetXmin(), prof_EM->GetXaxis()->GetXmax(),
                                prof_EM->GetNbinsY(), prof_EM->GetYaxis()->GetXmin(), prof_EM->GetYaxis()->GetXmax());
        frame->Draw();

        for (int ix = 1; ix <= prof_EM->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= prof_EM->GetNbinsY(); iy++) {
                double rms_EM = prof_EM->GetBinError(ix, iy);
                double rms_DA = prof_DA->GetBinError(ix, iy);
                double rms_Sigma = prof_Sigma->GetBinError(ix, iy);
                double entries_EM = prof_EM->GetBinEntries(prof_EM->GetBin(ix, iy));

                if (entries_EM < 10) continue;

                // Find best method (smallest RMS)
                double best_rms = rms_EM;
                int best_method = 0; // 0=EM, 1=DA, 2=Sigma
                if (rms_DA < best_rms && rms_DA > 0) { best_rms = rms_DA; best_method = 1; }
                if (rms_Sigma < best_rms && rms_Sigma > 0) { best_rms = rms_Sigma; best_method = 2; }

                double x = prof_EM->GetXaxis()->GetBinCenter(ix);
                double y = prof_EM->GetYaxis()->GetBinCenter(iy);

                TMarker* marker = new TMarker(x, y, 20);
                marker->SetMarkerSize(1.5);

                // Color code: Red=EM, Blue=DA, Green=Sigma
                if (best_method == 0) marker->SetMarkerColor(kRed);
                else if (best_method == 1) marker->SetMarkerColor(kBlue);
                else marker->SetMarkerColor(kGreen+2);

                marker->Draw();
            }
        }

        // Add legend
        TLegend* leg = new TLegend(0.7, 0.85, 0.9, 0.95);
        TMarker* m1 = new TMarker(0, 0, 20); m1->SetMarkerColor(kRed);
        TMarker* m2 = new TMarker(0, 0, 20); m2->SetMarkerColor(kBlue);
        TMarker* m3 = new TMarker(0, 0, 20); m3->SetMarkerColor(kGreen+2);
        leg->AddEntry(m1, "EM Best", "p");
        leg->AddEntry(m2, "DA Best", "p");
        leg->AddEntry(m3, "Sigma Best", "p");
        leg->Draw();

        c->SaveAs(saveName);
        delete c;
        delete frame;
    };

    createBestMethodPlot("Q2_RelRes_vs_xy_EM", "Q2_RelRes_vs_xy_DA", "Q2_RelRes_vs_xy_Sigma",
                         "figs/resolutions/2d_maps/Q2_RelRes_Q2x_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("x_RelRes_vs_xQ2_EM", "x_RelRes_vs_xQ2_DA", "x_RelRes_vs_xQ2_Sigma",
                         "figs/resolutions/2d_maps/x_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("y_RelRes_vs_xQ2_EM", "y_RelRes_vs_xQ2_DA", "y_RelRes_vs_xQ2_Sigma",
                         "figs/resolutions/2d_maps/y_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    // MX2 Resolution Heatmap (uses COLZ TEXT instead of circles)
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

    // E-pz 2D Correlation Plot
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

    //=================================================================
    // SECTION 7: DIFFERENTIAL CROSS SECTION PLOTS
    //=================================================================

    std::cout << "Creating differential cross section plots..." << std::endl;

    gSystem->mkdir("figs/cross_sections", kTRUE);

    // Create sum of B0 and RP reconstructed dsigma/dt
    TH1D* h_dsigma_dt_B0_temp = (TH1D*)inputFile->Get("dsigma_dt_B0");
    TH1D* h_dsigma_dt_RP_temp = (TH1D*)inputFile->Get("dsigma_dt_RP");

    if (h_dsigma_dt_B0_temp && h_dsigma_dt_RP_temp) {
        TH1D* h_dsigma_dt_Sum = (TH1D*)h_dsigma_dt_B0_temp->Clone("dsigma_dt_Sum");
        h_dsigma_dt_Sum->SetTitle("B0+RP Sum d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]");
        h_dsigma_dt_Sum->Add(h_dsigma_dt_RP_temp);
    }

    // Create exponential fit function: f(t) = A * e^(-b*|t|)
    TH1D* h_dsigma_dt_MC_temp = (TH1D*)inputFile->Get("dsigma_dt_MC");
    TF1* fit_exp = nullptr;
    double b_value = 0.0;
    double b_error = 0.0;

    if (h_dsigma_dt_MC_temp) {
        fit_exp = new TF1("fit_exp", "[0]*TMath::Exp(-[1]*x)", 0.01, 2.0);
        fit_exp->SetParameters(1000, 5.0);
        fit_exp->SetParNames("A", "b");
        fit_exp->SetLineColor(kMagenta+2);
        fit_exp->SetLineWidth(3);
        fit_exp->SetLineStyle(2);

        TFitResultPtr fitResult = h_dsigma_dt_MC_temp->Fit(fit_exp, "RSQ");

        b_value = fit_exp->GetParameter(1);
        b_error = fit_exp->GetParError(1);

        std::cout << "Exponential fit: b = " << b_value << " +/- " << b_error << " GeV^{-2}" << std::endl;
    }

    // Create custom dsigma/dt plots with exponential fit
    TH1D* h_dsigma_dt_MC = (TH1D*)inputFile->Get("dsigma_dt_MC");
    TH1D* h_dsigma_dt_B0 = (TH1D*)inputFile->Get("dsigma_dt_B0");
    TH1D* h_dsigma_dt_RP = (TH1D*)inputFile->Get("dsigma_dt_RP");
    TH1D* h_dsigma_dt_Sum = (TH1D*)gDirectory->Get("dsigma_dt_Sum");

    if (h_dsigma_dt_MC && h_dsigma_dt_B0 && h_dsigma_dt_RP && h_dsigma_dt_Sum && fit_exp) {
        // dsigma/dt plot with linear y-axis
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
        h_dsigma_dt_Sum->SetMarkerColor(kOrange+7);
        h_dsigma_dt_Sum->SetLineColor(kOrange+7);
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

        // dsigma/dt plot with log y-axis
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

    // Triple differential cross section plots (d³σ/(dQ² dβ dx_pom))
    TH3D* h_d3sigma_MC = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_MC");
    TH3D* h_d3sigma_B0 = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_B0");
    TH3D* h_d3sigma_RP = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_RP");

    if (h_d3sigma_MC && h_d3sigma_B0 && h_d3sigma_RP) {
        int n_Q2_bins = h_d3sigma_MC->GetNbinsX();
        int n_beta_bins = h_d3sigma_MC->GetNbinsY();
        int n_xpom_bins = h_d3sigma_MC->GetNbinsZ();

        // Create grid canvas for d³σ vs beta (varying Q², x_pom)
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

                proj_MC->SetTitle(Form("Q^{2} bin %d, x_{pom} bin %d;#beta;d^{3}#sigma", i_Q2, i_xpom));
                proj_MC->Draw("HIST");
                proj_B0->Draw("PE SAME");
                proj_RP->Draw("PE SAME");
            }
        }

        c_beta->SaveAs("figs/cross_sections/d3sigma_vs_beta.png");
        delete c_beta;

        // Create grid canvas for d³σ vs x_pom (varying Q², beta)
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

                proj_MC->SetTitle(Form("Q^{2} bin %d, #beta bin %d;x_{pom};d^{3}#sigma", i_Q2, i_beta));
                proj_MC->Draw("HIST");
                proj_B0->Draw("PE SAME");
                proj_RP->Draw("PE SAME");
            }
        }

        c_xpom->SaveAs("figs/cross_sections/d3sigma_vs_xpom.png");
        delete c_xpom;
    }

    //=================================================================
    // SECTION 8: ACCEPTANCE, EFFICIENCY, AND PURITY PLOTS
    //=================================================================

    std::cout << "Creating acceptance, efficiency, and purity plots..." << std::endl;

    gSystem->mkdir("figs/performance", kTRUE);

    // Calculate acceptance, efficiency, and purity from tracking histograms
    TH1D* h_gen_Q2 = (TH1D*)inputFile->Get("h_gen_Q2");
    TH1D* h_gen_and_reco_after_cuts_Q2_EM = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_EM");
    TH1D* h_gen_and_reco_after_cuts_Q2_DA = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_DA");
    TH1D* h_gen_and_reco_after_cuts_Q2_Sigma = (TH1D*)inputFile->Get("h_gen_and_reco_after_cuts_Q2_Sigma");
    TH1D* h_reco_Q2_EM = (TH1D*)inputFile->Get("h_reco_Q2_EM");
    TH1D* h_reco_Q2_DA = (TH1D*)inputFile->Get("h_reco_Q2_DA");
    TH1D* h_reco_Q2_Sigma = (TH1D*)inputFile->Get("h_reco_Q2_Sigma");

    if (h_gen_Q2 && h_gen_and_reco_after_cuts_Q2_EM && h_reco_Q2_EM) {
        // Create acceptance histograms
        TH1D* h_acceptance_Q2_EM = (TH1D*)h_gen_and_reco_after_cuts_Q2_EM->Clone("h_acceptance_Q2_EM");
        h_acceptance_Q2_EM->Divide(h_gen_Q2);
        h_acceptance_Q2_EM->SetTitle("Acceptance vs Q^{2} (EM);Q^{2} [GeV^{2}];Acceptance");

        TH1D* h_acceptance_Q2_DA = (TH1D*)h_gen_and_reco_after_cuts_Q2_DA->Clone("h_acceptance_Q2_DA");
        h_acceptance_Q2_DA->Divide(h_gen_Q2);

        TH1D* h_acceptance_Q2_Sigma = (TH1D*)h_gen_and_reco_after_cuts_Q2_Sigma->Clone("h_acceptance_Q2_Sigma");
        h_acceptance_Q2_Sigma->Divide(h_gen_Q2);

        // Create purity histograms
        TH1D* h_purity_Q2_EM = (TH1D*)h_gen_and_reco_after_cuts_Q2_EM->Clone("h_purity_Q2_EM");
        h_purity_Q2_EM->Divide(h_reco_Q2_EM);
        h_purity_Q2_EM->SetTitle("Purity vs Q^{2} (EM);Q^{2} [GeV^{2}];Purity");

        TH1D* h_purity_Q2_DA = (TH1D*)h_gen_and_reco_after_cuts_Q2_DA->Clone("h_purity_Q2_DA");
        h_purity_Q2_DA->Divide(h_reco_Q2_DA);

        TH1D* h_purity_Q2_Sigma = (TH1D*)h_gen_and_reco_after_cuts_Q2_Sigma->Clone("h_purity_Q2_Sigma");
        h_purity_Q2_Sigma->Divide(h_reco_Q2_Sigma);

        // Plot acceptance, efficiency, and purity
        plot_ptr = new PlotOptions1D(
            {"h_acceptance_Q2_EM", "h_acceptance_Q2_DA", "h_acceptance_Q2_Sigma"},
            {"EM Method", "DA Method", "Sigma Method"},
            {"pe", "pe", "pe"},
            "Acceptance vs Q^{2}",
            "Q^{2} [GeV^{2}]",
            "Acceptance",
            "figs/performance/acceptance_vs_Q2.png",
            true,   // logX
            false   // logY
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
            true,   // logX
            false   // logY
        );
        plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
        plots.push_back(plot_ptr);
    }

    //=================================================================
    // EXECUTE ALL STANDARD PLOTS (using PlotOptions framework)
    //=================================================================

    std::cout << "\n========================================" << std::endl;
    std::cout << "Executing all standard plots..." << std::endl;
    std::cout << "Total number of plots: " << plots.size() << std::endl;
    std::cout << "========================================\n" << std::endl;

    for (size_t i = 0; i < plots.size(); i++) {
        std::cout << "Plot " << (i+1) << "/" << plots.size() << "..." << std::endl;
        plots[i]->Plot(inputFile);
    }

    // Clean up
    for (const auto& plot : plots) {
        delete plot;
    }

    inputFile->Close();
    delete inputFile;

    std::cout << "\n========================================" << std::endl;
    std::cout << "PLOTTING COMPLETE!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nAll plots saved to figs/ directory with subdirectories:" << std::endl;
    std::cout << "  - figs/distributions/          : Base distributions" << std::endl;
    std::cout << "  - figs/response_matrices/      : 2D correlation plots" << std::endl;
    std::cout << "  - figs/resolutions/simple/     : Overall 1D resolutions" << std::endl;
    std::cout << "  - figs/resolutions/binned/     : Binned resolutions with Gaussian fits" << std::endl;
    std::cout << "  - figs/resolutions/binned/bins/: Individual bin plots" << std::endl;
    std::cout << "  - figs/resolutions/2d_maps/    : 2D resolution circle maps" << std::endl;
    std::cout << "  - figs/cross_sections/         : Differential cross sections" << std::endl;
    std::cout << "  - figs/performance/            : Acceptance, efficiency, purity" << std::endl;
    std::cout << "\n========================================\n" << std::endl;

    return 0;
}
