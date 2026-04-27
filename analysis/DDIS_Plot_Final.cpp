// Inclusive DIS plotter
// g++ DDIS_Plot_Final.cpp -o DDIS_Plot_Final $(root-config --cflags --glibs)
// ./DDIS_Plot_Final <combined.root>

#include "Plotting.hpp"
#include "RecoMethods.hpp"
#include "YAMLBinning.hpp"
#include "PlotDrawing.hpp"
#include "GridDrawing.hpp"
#include "Utility.hpp"
#include "plots/DDISPlots.hpp"

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
#include <TParameter.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <map>
#include <limits>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>

static bool InferBeamEnergiesFromInputPath(const std::string& inputPath,
                                           double& eBeamGeV,
                                           double& pBeamGeV,
                                           std::string* firstMatchedToken = nullptr,
                                           std::string* mismatchToken = nullptr) {
    std::vector<std::string> tokens;
    tokens.reserve(16);
    tokens.push_back(inputPath);

    std::string token;
    token.reserve(inputPath.size());
    for (char c : inputPath) {
        if (c == '/' || c == '\\') {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            continue;
        }
        token.push_back(c);
    }
    if (!token.empty()) tokens.push_back(token);

    bool found = false;
    double eRef = 0.0;
    double pRef = 0.0;
    std::string firstToken;
    for (const auto& tok : tokens) {
        double eCand = 0.0;
        double pCand = 0.0;
        if (!ParseBeamEnergiesFromFilename(tok, eCand, pCand)) continue;

        if (!found) {
            found = true;
            eRef = eCand;
            pRef = pCand;
            firstToken = tok;
            continue;
        }

        const bool sameE = std::fabs(eCand - eRef) < 1e-9;
        const bool sameP = std::fabs(pCand - pRef) < 1e-9;
        if (!sameE || !sameP) {
            if (mismatchToken) *mismatchToken = tok;
            return false;
        }
    }

    if (!found) return false;
    eBeamGeV = eRef;
    pBeamGeV = pRef;
    if (firstMatchedToken) *firstMatchedToken = firstToken;
    return true;
}

static bool EstimateDISsFromQ2Tree(TFile* inputFile, double& sGeV2) {
    if (!inputFile) return false;
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Q2_tree"));
    if (!tree) return false;
    if (!tree->GetBranch("Q2_truth") || !tree->GetBranch("x_truth") || !tree->GetBranch("y_truth")) {
        return false;
    }

    float q2 = -1.0f;
    float x = -1.0f;
    float y = -1.0f;
    tree->SetBranchAddress("Q2_truth", &q2);
    tree->SetBranchAddress("x_truth", &x);
    tree->SetBranchAddress("y_truth", &y);

    const Long64_t nEntries = tree->GetEntries();
    if (nEntries <= 0) return false;
    const Long64_t maxSample = 200000;
    const Long64_t stride = std::max<Long64_t>(1, nEntries / maxSample);

    std::vector<double> sValues;
    sValues.reserve(static_cast<size_t>(nEntries / stride + 1));
    for (Long64_t i = 0; i < nEntries; i += stride) {
        tree->GetEntry(i);
        if (!std::isfinite(q2) || !std::isfinite(x) || !std::isfinite(y)) continue;
        if (q2 <= 0.0f || x <= 0.0f || y <= 0.0f) continue;
        const double sCand = static_cast<double>(q2) / (static_cast<double>(x) * static_cast<double>(y));
        if (!std::isfinite(sCand) || sCand <= 0.0) continue;
        sValues.push_back(sCand);
    }

    if (sValues.empty()) return false;
    auto mid = sValues.begin() + static_cast<long>(sValues.size() / 2);
    std::nth_element(sValues.begin(), mid, sValues.end());
    sGeV2 = *mid;
    return std::isfinite(sGeV2) && sGeV2 > 0.0;
}














// =====================================================================
// Task 1: EPz distribution with cut boundary lines
// =====================================================================

// =====================================================================
// Task 2: Generic resolution comparison (Fitted + RMS for all methods)
// =====================================================================

// =====================================================================
// Task 3B: Plot NxN 3D response matrix + diagonal fraction
// =====================================================================

// =====================================================================
// Task 5: Phase-space plot with iso-y and iso-W2 lines
// =====================================================================

// Compares the EICrecon ScatteredElectronsTruth_objIdx electron ("old")
// against the ElectronID §A.1 selection ("eid"). Consumes histograms written
// by DDIS_Skim_Final: Ep_e / pT_e / phi_e / EPz_reco_mc (old) and
// Ep_e_eid / pT_e_eid / phi_e_eid / EPz_eid (eid), plus 2D correlations
// and the category counter.

int main(int argc, char** argv) {
    if (argc < 3) {
        Logger::error(std::string("Usage: ") + argv[0] + " <combined.root> <bins.yaml>");
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
        Logger::error("Could not open file " + std::string(inputFileName));
        return 1;
    }

    double nGenFromFile = -1.0;
    if (auto* nGenParam = dynamic_cast<TParameter<double>*>(inputFile->Get("N_gen"))) {
        nGenFromFile = nGenParam->GetVal();
    } else if (auto* nEvtParam = dynamic_cast<TParameter<Long64_t>*>(inputFile->Get("nEventsProcessed"))) {
        nGenFromFile = static_cast<double>(nEvtParam->GetVal());
    }
    double sigmaTotalNbFromFile = -1.0;
    if (auto* sigmaParam = dynamic_cast<TParameter<double>*>(inputFile->Get("sigma_total_nb"))) {
        sigmaTotalNbFromFile = sigmaParam->GetVal();
    }
    double yMinCutFromFile = 0.01;
    double yMaxCutFromFile = 0.95;
    if (auto* yMinParam = dynamic_cast<TParameter<double>*>(inputFile->Get("DISCut_y_min"))) {
        yMinCutFromFile = yMinParam->GetVal();
    }
    if (auto* yMaxParam = dynamic_cast<TParameter<double>*>(inputFile->Get("DISCut_y_max"))) {
        yMaxCutFromFile = yMaxParam->GetVal();
    }
    if (nGenFromFile > 0.0) {
        Logger::info("Loaded N_gen from output file: " + std::to_string(nGenFromFile));
    } else {
        Logger::info("N_gen metadata not found in output file.");
    }
    if (sigmaTotalNbFromFile > 0.0) {
        Logger::info("Loaded sigma_total from output file: " + std::to_string(sigmaTotalNbFromFile) + " nb");
    }
    Logger::info("Using DIS y-window for reduced-cross-section plot: " +
                 std::to_string(yMinCutFromFile) + " < y < " + std::to_string(yMaxCutFromFile));

    double eBeamGeV = 0.0;
    double pBeamGeV = 0.0;
    double disSGeV2 = 0.0;
    std::string firstBeamToken;
    std::string mismatchBeamToken;
    const bool haveBeamFromPath = InferBeamEnergiesFromInputPath(
        inputFileName.Data(), eBeamGeV, pBeamGeV, &firstBeamToken, &mismatchBeamToken
    );
    if (haveBeamFromPath && eBeamGeV > 0.0 && pBeamGeV > 0.0) {
        disSGeV2 = 4.0 * eBeamGeV * pBeamGeV;
        Logger::info("Reduced-cross-section kinematics: parsed beam tag " +
                     std::to_string(eBeamGeV) + "x" + std::to_string(pBeamGeV) +
                     " GeV from token '" + firstBeamToken +
                     "', using s = " + std::to_string(disSGeV2) + " GeV^2.");
    } else {
        if (!mismatchBeamToken.empty()) {
            Logger::warning("Inconsistent beam-energy tags found in input path (first match '" +
                            firstBeamToken + "', mismatch '" + mismatchBeamToken +
                            "'). Falling back to kinematic s estimate from Q2_tree.");
        } else {
            Logger::info("Beam-energy tag not found in input path; estimating s from Q2_tree truth kinematics.");
        }
        if (EstimateDISsFromQ2Tree(inputFile, disSGeV2)) {
            Logger::info("Estimated s = " + std::to_string(disSGeV2) +
                         " GeV^2 from median(Q2_truth/(x_truth*y_truth)).");
        } else {
            Logger::warning("Could not determine s from beam tags or Q2_tree. "
                            "Reduced-cross-section stacked plots will be skipped.");
        }
    }

    std::string yamlPath = argv[2];
    std::vector<BinDef> yaml_bins = ReadBinsFromYAML(yamlPath);
    std::vector<double> q2_edges;
    std::vector<double> beta_edges;
    std::vector<double> xpom_edges;
    if (!yaml_bins.empty()) {
        CollectEdges(yaml_bins, q2_edges, beta_edges, xpom_edges);
    } else {
        Logger::error("No 3D_bins loaded from " + yamlPath);
        return 1;
    }

    std::vector<PlotOptions*> plots;
    PlotOptions1D* plot_ptr = nullptr;
    PlotOptionsBinnedRelRes* binned_plot_ptr = nullptr;

    // Create output directory structure
    gSystem->mkdir("figs/inclusive/distributions", kTRUE);
    gSystem->mkdir("figs/inclusive/resolutions/simple", kTRUE);
    gSystem->mkdir("figs/inclusive/resolutions/binned/bins", kTRUE);
    gSystem->mkdir("figs/inclusive/resolutions/2d_maps", kTRUE);
    gSystem->mkdir("figs/inclusive/response", kTRUE);
    gSystem->mkdir("figs/inclusive/control", kTRUE);
    gSystem->mkdir("figs/inclusive/efficiency", kTRUE);
    gSystem->mkdir("figs/inclusive/performance", kTRUE);
    gSystem->mkdir("figs/diffractive/distributions", kTRUE);
    gSystem->mkdir("figs/diffractive/resolutions/simple", kTRUE);
    gSystem->mkdir("figs/diffractive/resolutions/binned/bins", kTRUE);
    gSystem->mkdir("figs/diffractive/resolutions/2d_maps", kTRUE);
    gSystem->mkdir("figs/diffractive/response", kTRUE);
    gSystem->mkdir("figs/diffractive/control", kTRUE);
    gSystem->mkdir("figs/diffractive/efficiency", kTRUE);
    gSystem->mkdir("figs/diffractive/performance", kTRUE);
    gSystem->mkdir("figs/cross_sections/debug", kTRUE);
    gSystem->mkdir("figs/electron_id", kTRUE);

    // Scattered-electron finder comparison: old (ScatteredElectronsTruth_objIdx)
    // vs ElectronID paper A.1 selection. Consumes the _eid / _old_vs_eid /
    // _finder_category histograms written by DDIS_Skim_Final.
    PlotElectronIDComparison(inputFile, "figs/electron_id");

    // M_X^2 hadronic-sum vs kinematic comparison (truth and reco). Hand-styled
    // colors so the four curves are visually distinct.
    gSystem->mkdir("figs/diffractive/distributions", kTRUE);
    PlotMX2Comparison(inputFile, "figs/diffractive/distributions/MX2_comparison.png",      false);
    PlotMX2Comparison(inputFile, "figs/diffractive/distributions/MX2_comparison_logy.png", true);

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

    // ---- Uncorrected MC vs pseudo-data control plots ----
    // Inclusive variables (raw counts)
    PlotRecoSetUncorrected(inputFile, "Q2_reco_mc", "Q2_reco_pdata",
                           "Q^{2} [GeV^{2}]", "Q^{2} Uncorrected: MC vs Pseudo-data",
                           "figs/inclusive/control/q2_mc_vs_pdata.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "x_reco_mc", "x_reco_pdata",
                           "x_{Bj}", "x_{Bj} Uncorrected: MC vs Pseudo-data",
                           "figs/inclusive/control/xbj_mc_vs_pdata.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "y_reco_mc", "y_reco_pdata",
                           "y", "y Uncorrected: MC vs Pseudo-data",
                           "figs/inclusive/control/y_mc_vs_pdata.png",
                           false, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "W2_reco_mc", "W2_reco_pdata",
                           "W^{2} [GeV^{2}]", "W^{2} Uncorrected: MC vs Pseudo-data",
                           "figs/inclusive/control/w2_mc_vs_pdata.png",
                           true, true, controlSetALabel, controlSetBLabel);
    // Inclusive variables (PDF-normalized)
    PlotRecoSetUncorrected(inputFile, "Q2_reco_mc", "Q2_reco_pdata",
                           "Q^{2} [GeV^{2}]", "Q^{2} Uncorrected: MC vs Pseudo-data (PDF)",
                           "figs/inclusive/control/q2_mc_vs_pdata_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "x_reco_mc", "x_reco_pdata",
                           "x_{Bj}", "x_{Bj} Uncorrected: MC vs Pseudo-data (PDF)",
                           "figs/inclusive/control/xbj_mc_vs_pdata_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "y_reco_mc", "y_reco_pdata",
                           "y", "y Uncorrected: MC vs Pseudo-data (PDF)",
                           "figs/inclusive/control/y_mc_vs_pdata_pdf.png",
                           false, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "W2_reco_mc", "W2_reco_pdata",
                           "W^{2} [GeV^{2}]", "W^{2} Uncorrected: MC vs Pseudo-data (PDF)",
                           "figs/inclusive/control/w2_mc_vs_pdata_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    // Diffractive variables (raw counts)
    PlotRecoSetUncorrected(inputFile, "t_reco_mc_B0", "t_reco_pdata_B0",
                           "|t| [GeV^{2}]", "|t| Uncorrected: MC vs Pseudo-data (B0)",
                           "figs/diffractive/control/t_mc_vs_pdata_b0.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "t_reco_mc_RP", "t_reco_pdata_RP",
                           "|t| [GeV^{2}]", "|t| Uncorrected: MC vs Pseudo-data (RP)",
                           "figs/diffractive/control/t_mc_vs_pdata_rp.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "beta_reco_mc_B0", "beta_reco_pdata_B0",
                           "#beta", "#beta Uncorrected: MC vs Pseudo-data (B0)",
                           "figs/diffractive/control/beta_mc_vs_pdata_b0.png",
                           false, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "beta_reco_mc_RP", "beta_reco_pdata_RP",
                           "#beta", "#beta Uncorrected: MC vs Pseudo-data (RP)",
                           "figs/diffractive/control/beta_mc_vs_pdata_rp.png",
                           false, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "xpom_reco_mc_B0", "xpom_reco_pdata_B0",
                           "x_{pom}", "x_{pom} Uncorrected: MC vs Pseudo-data (B0)",
                           "figs/diffractive/control/xpom_mc_vs_pdata_b0.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "xpom_reco_mc_RP", "xpom_reco_pdata_RP",
                           "x_{pom}", "x_{pom} Uncorrected: MC vs Pseudo-data (RP)",
                           "figs/diffractive/control/xpom_mc_vs_pdata_rp.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "MX2_reco_mc", "MX2_reco_pdata",
                           "M_{X}^{2} [GeV^{2}]", "M_{X}^{2} Uncorrected: MC vs Pseudo-data",
                           "figs/diffractive/control/mx2_mc_vs_pdata.png",
                           true, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "theta_reco_mc_B0", "theta_reco_pdata_B0",
                           "#theta [mrad]", "#theta Uncorrected: MC vs Pseudo-data (B0)",
                           "figs/diffractive/control/theta_mc_vs_pdata_b0.png",
                           false, true, controlSetALabel, controlSetBLabel);
    PlotRecoSetUncorrected(inputFile, "theta_reco_mc_RP", "theta_reco_pdata_RP",
                           "#theta [mrad]", "#theta Uncorrected: MC vs Pseudo-data (RP)",
                           "figs/diffractive/control/theta_mc_vs_pdata_rp.png",
                           false, true, controlSetALabel, controlSetBLabel);
    // Diffractive variables (PDF-normalized)
    PlotRecoSetUncorrected(inputFile, "t_reco_mc_B0", "t_reco_pdata_B0",
                           "|t| [GeV^{2}]", "|t| Uncorrected: MC vs Pseudo-data (B0, PDF)",
                           "figs/diffractive/control/t_mc_vs_pdata_b0_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "t_reco_mc_RP", "t_reco_pdata_RP",
                           "|t| [GeV^{2}]", "|t| Uncorrected: MC vs Pseudo-data (RP, PDF)",
                           "figs/diffractive/control/t_mc_vs_pdata_rp_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "beta_reco_mc_B0", "beta_reco_pdata_B0",
                           "#beta", "#beta Uncorrected: MC vs Pseudo-data (B0, PDF)",
                           "figs/diffractive/control/beta_mc_vs_pdata_b0_pdf.png",
                           false, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "beta_reco_mc_RP", "beta_reco_pdata_RP",
                           "#beta", "#beta Uncorrected: MC vs Pseudo-data (RP, PDF)",
                           "figs/diffractive/control/beta_mc_vs_pdata_rp_pdf.png",
                           false, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "xpom_reco_mc_B0", "xpom_reco_pdata_B0",
                           "x_{pom}", "x_{pom} Uncorrected: MC vs Pseudo-data (B0, PDF)",
                           "figs/diffractive/control/xpom_mc_vs_pdata_b0_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "xpom_reco_mc_RP", "xpom_reco_pdata_RP",
                           "x_{pom}", "x_{pom} Uncorrected: MC vs Pseudo-data (RP, PDF)",
                           "figs/diffractive/control/xpom_mc_vs_pdata_rp_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "MX2_reco_mc", "MX2_reco_pdata",
                           "M_{X}^{2} [GeV^{2}]", "M_{X}^{2} Uncorrected: MC vs Pseudo-data (PDF)",
                           "figs/diffractive/control/mx2_mc_vs_pdata_pdf.png",
                           true, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "theta_reco_mc_B0", "theta_reco_pdata_B0",
                           "#theta [mrad]", "#theta Uncorrected: MC vs Pseudo-data (B0, PDF)",
                           "figs/diffractive/control/theta_mc_vs_pdata_b0_pdf.png",
                           false, false, controlSetALabel, controlSetBLabel, true);
    PlotRecoSetUncorrected(inputFile, "theta_reco_mc_RP", "theta_reco_pdata_RP",
                           "#theta [mrad]", "#theta Uncorrected: MC vs Pseudo-data (RP, PDF)",
                           "figs/diffractive/control/theta_mc_vs_pdata_rp_pdf.png",
                           false, false, controlSetALabel, controlSetBLabel, true);

    PlotDensityFromHistWithOverlay(inputFile, "beta_Q2_reco", "#beta", "Q^{2} [GeV^{2}]",
                                   "figs/diffractive/distributions/beta_q2_density_reco.png",
                                   false, true, beta_edges, q2_edges, false, true);
    PlotDensityFromHist(inputFile, "t_Q2_reco", "|t| [GeV^{2}]", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/distributions/t_q2_density_reco.png", true, true);
    PlotDensityFromHistWithOverlay(inputFile, "xpom_Q2_reco", "x_{pom}", "Q^{2} [GeV^{2}]",
                                   "figs/diffractive/distributions/xpom_q2_density_reco.png",
                                   true, true, xpom_edges, q2_edges, true, true);
    PlotDensityFromHist(inputFile, "beta_t_reco", "#beta", "|t| [GeV^{2}]",
                        "figs/diffractive/distributions/beta_t_density_reco.png", false, true);
    PlotDensityFromHist(inputFile, "xbj_t_reco", "x_{Bj}", "|t| [GeV^{2}]",
                        "figs/diffractive/distributions/xbj_t_density_reco.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_t_reco", "x_{pom}", "|t| [GeV^{2}]",
                        "figs/diffractive/distributions/xpom_t_density_reco.png", true, true);
    PlotDensityFromHistWithOverlay(inputFile, "xpom_beta_reco", "x_{pom}", "#beta",
                                   "figs/diffractive/distributions/xpom_beta_density_reco.png",
                                   true, false, xpom_edges, beta_edges, true, false);
    PlotDensityFromHist(inputFile, "xbj_beta_reco", "x_{Bj}", "#beta",
                        "figs/diffractive/distributions/xbj_beta_density_reco.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_xpom_reco", "x_{Bj}", "x_{pom}",
                        "figs/diffractive/distributions/xbj_xpom_density_reco.png", true, true);
    // Inclusive truth-level 2D density
    PlotDensityFromHist(inputFile, "yQ2_truth", "y", "Q^{2} [GeV^{2}]",
                        "figs/inclusive/distributions/y_q2_density_truth.png", false, true);

    // MC truth density versions
    PlotDensityFromHist(inputFile, "beta_Q2_truth", "#beta", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/distributions/beta_q2_density_truth.png", false, true);
    PlotDensityFromHist(inputFile, "t_Q2_truth", "|t| [GeV^{2}]", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/distributions/t_q2_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_Q2_truth", "x_{pom}", "Q^{2} [GeV^{2}]",
                        "figs/diffractive/distributions/xpom_q2_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "beta_t_truth", "#beta", "|t| [GeV^{2}]",
                        "figs/diffractive/distributions/beta_t_density_truth.png", false, true);
    PlotDensityFromHist(inputFile, "xbj_t_truth", "x_{Bj}", "|t| [GeV^{2}]",
                        "figs/diffractive/distributions/xbj_t_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_t_truth", "x_{pom}", "|t| [GeV^{2}]",
                        "figs/diffractive/distributions/xpom_t_density_truth.png", true, true);
    PlotDensityFromHist(inputFile, "xpom_beta_truth", "x_{pom}", "#beta",
                        "figs/diffractive/distributions/xpom_beta_density_truth.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_beta_truth", "x_{Bj}", "#beta",
                        "figs/diffractive/distributions/xbj_beta_density_truth.png", true, false);
    PlotDensityFromHist(inputFile, "xbj_xpom_truth", "x_{Bj}", "x_{pom}",
                        "figs/diffractive/distributions/xbj_xpom_density_truth.png", true, true);
    PlotPhaseSpaceSlices(inputFile, yamlPath);

    PlotGraphDensity(inputFile,
                     "g_W2_EM",
                     "W^{2} Correlation (EM, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/distributions/w2_corr_unbinned_em_density.png",
                     10.0, 1.0e4, 140, true, true);
    PlotGraphDensity(inputFile,
                     "g_W2_DA",
                     "W^{2} Correlation (DA, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/distributions/w2_corr_unbinned_da_density.png",
                     10.0, 1.0e4, 140, true, true);
    PlotGraphDensity(inputFile,
                     "g_W2_Best",
                     "W^{2} Correlation (Best, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/distributions/w2_corr_unbinned_best_density.png",
                     10.0, 1.0e4, 140, true, true);
    PlotGraphDensity(inputFile,
                     "g_W2_Sigma",
                     "W^{2} Correlation (#Sigma, Unbinned)",
                     "W^{2}_{truth} [GeV^{2}]",
                     "W^{2}_{reco} [GeV^{2}]",
                     "figs/inclusive/distributions/w2_corr_unbinned_sigma_density.png",
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
            "figs/inclusive/performance/acceptance_vs_Q2.png",
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
            "figs/inclusive/performance/purity_vs_Q2.png",
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
        "Number of events",
        "figs/inclusive/distributions/q2_methods_hist.png",
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
        "figs/inclusive/distributions/q2_methods_pdf.png",
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
        "Number of events",
        "figs/inclusive/distributions/xbj_methods_hist.png",
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
        "figs/inclusive/distributions/xbj_methods_pdf.png",
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
        "Number of events",
        "figs/inclusive/distributions/y_methods_hist.png",
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
        "figs/inclusive/distributions/y_methods_pdf.png",
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
        {"W2_truth", "W2_EM", "W2_DA", "W2_Best", "W2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Best", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe", "pe"},
        "W^{2} Reconstruction Methods",
        "W^{2} [GeV^{2}]",
        "Number of events",
        "figs/inclusive/distributions/w2_methods_hist.png",
        true,
        true
    );
    plot_ptr->SetRangeX(10.0, 1.0e4);
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.3, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"W2_truth", "W2_EM", "W2_DA", "W2_Best", "W2_Sigma"},
        {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Best", "Reco. Sigma"},
        {"hist", "pe", "pe", "pe", "pe"},
        "W^{2} PDF Comparison",
        "W^{2} [GeV^{2}]",
        "PDF",
        "figs/inclusive/distributions/w2_methods_pdf.png",
        true,
        true,
        true
    );
    plot_ptr->SetRangeX(10.0, 1.0e4);
    plot_ptr->SetLegendPosition(0.7, 0.7, 0.9, 0.9);
    plots.push_back(plot_ptr);

    // Fine-binned W^{2} (200 log bins, 1-1e4 GeV^{2}) — for informed choice
    // of the first analysis-bin edge. Use points only (no fills) so overlaid
    // curves are distinguishable.
    {
        auto* fine_ptr = new PlotOptions1D(
            {"W2_truth_fine", "W2_EM_fine", "W2_DA_fine", "W2_Best_fine", "W2_Sigma_fine", "W2_ESigma_fine"},
            {"MC: truth", "Reco. EM", "Reco. DA", "Reco. Best", "Reco. Sigma", "Reco. ESigma"},
            {"hist", "pe", "pe", "pe", "pe", "pe"},
            "W^{2} Distribution (Fine Binning)",
            "W^{2} [GeV^{2}]",
            "Number of events",
            "figs/inclusive/distributions/w2_methods_fine.png",
            true,
            true
        );
        fine_ptr->SetDisableFills(true);
        fine_ptr->SetRangeX(1.0, 1.0e4);
        fine_ptr->SetLegendPosition(0.15, 0.7, 0.35, 0.9);
        plots.push_back(fine_ptr);
    }

    // 2D: relative resolution of Best-blend W^{2} vs truth W^{2} (fine).
    // Use to see where reco starts diverging from truth on the low-W^{2} side.
    PlotDensityFromHist(inputFile,
                        "W2_RelRes_Best_fine",
                        "W^{2}_{MC} [GeV^{2}]",
                        "(W^{2}_{best}-W^{2}_{MC})/W^{2}_{MC}",
                        "figs/inclusive/resolutions/w2_relres_best_vs_truth_fine.png",
                        true,
                        false);

    // =================================================================
    // NEW: Scattered electron leptonic quantities
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"Ep_e_truth", "Ep_e"},
        {"MC Truth", "Reconstruction"},
        {"hist", "E1"},
        "Scattered Electron Energy",
        "E'_{e} [GeV]",
        "Number of events",
        "figs/inclusive/distributions/electron_energy.png",
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
        "Number of events",
        "figs/inclusive/distributions/electron_energy_logy.png",
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
        "Number of events",
        "figs/inclusive/distributions/electron_phi.png",
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
        "Number of events",
        "figs/inclusive/distributions/electron_pt.png",
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
        "Number of events",
        "figs/inclusive/distributions/electron_pt_logy.png",
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
        "figs/inclusive/distributions/q2_corr_unbinned_em.png",
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
        "figs/inclusive/distributions/q2_corr_unbinned_da.png",
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
        "figs/inclusive/distributions/q2_corr_unbinned_sigma.png",
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
        "figs/inclusive/distributions/xbj_corr_unbinned_em.png",
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
        "figs/inclusive/distributions/xbj_corr_unbinned_da.png",
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
        "figs/inclusive/distributions/xbj_corr_unbinned_sigma.png",
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
        "figs/inclusive/distributions/y_corr_unbinned_em.png",
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
        "figs/inclusive/distributions/y_corr_unbinned_da.png",
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
        "figs/inclusive/distributions/y_corr_unbinned_sigma.png",
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
        "figs/inclusive/distributions/w2_corr_unbinned_em.png",
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
        "figs/inclusive/distributions/w2_corr_unbinned_da.png",
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
        "figs/inclusive/distributions/w2_corr_unbinned_sigma.png",
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
        "figs/inclusive/distributions/w2_corr_unbinned_all.png",
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
        "figs/inclusive/distributions/electron_energy_corr_unbinned.png",
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
        "figs/inclusive/distributions/electron_phi_corr_unbinned.png",
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
        "figs/inclusive/distributions/electron_pt_corr_unbinned.png",
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
        "Response_W2_Best",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_best.png",
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
        "Response_W2_Best_rowNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_best_rowNorm.png",
        true,
        true,
        {10.0, 1.0e4},
        {10.0, 1.0e4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Response_W2_Best_colNorm",
        "W^{2}_{truth} [GeV^{2}]",
        "W^{2}_{reco} [GeV^{2}]",
        "figs/inclusive/response/response_w2_best_colNorm.png",
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
        "Number of events",
        -0.05, 0.05,
        "figs/inclusive/resolutions/q2_relres_em.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_DA",
        "#frac{Q^{2}_{DA} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Number of events",
        -0.1, 0.1,
        "figs/inclusive/resolutions/q2_relres_da.png",
        "double_sided_crystalball"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "Q2_RelRes_Sigma",
        "#frac{Q^{2}_{#Sigma} - Q^{2}_{MC}}{ Q^{2}_{MC}}",
        "Number of events",
        -0.2, 0.05,
        "figs/inclusive/resolutions/q2_relres_sigma.png",
        "crystalball"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_EM",
        "#frac{x_{EM} - x_{MC}}{ x_{MC}}",
        "Number of events",
        -0.1, 0.1,
        "figs/inclusive/resolutions/xbj_relres_em.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_DA",
        "#frac{x_{DA} - x_{MC}}{X_{MC}}",
        "Number of events",
        -0.4, 0.5,
        "figs/inclusive/resolutions/xbj_relres_da.png",
        "double_sided_crystalball"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "x_RelRes_Sigma",
        "#frac{x_{#Sigma} - x_{MC}}{X_{MC}}",
        "Number of events",
        -0.5, 0.5,
        "figs/inclusive/resolutions/xbj_relres_sigma.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_EM",
        "#frac{y_{EM} - y_{MC}}{ y_{MC}}",
        "Number of events",
        -0.05, 0.05,
        "figs/inclusive/resolutions/y_relres_em.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_DA",
        "#frac{y_{DA} - y_{MC}}{y_{MC}}",
        "Number of events",
        -0.22, 0.3,
        "figs/inclusive/resolutions/y_relres_da.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "y_RelRes_Sigma",
        "#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "Number of events",
        -0.5, 0.5,
        "figs/inclusive/resolutions/y_relres_sigma.png",
        "double_sided_crystalball"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_EM",
        "#frac{W^{2}_{EM} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Number of events",
        -0.1, 0.1,
        "figs/inclusive/resolutions/w2_relres_em.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_DA",
        "#frac{W^{2}_{DA} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Number of events",
        -0.25, 0.3,
        "figs/inclusive/resolutions/w2_relres_da.png",
        "double_sided_crystalball"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_Best",
        "#frac{W^{2}_{best} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Number of events",
        -0.12, 0.12,
        "figs/inclusive/resolutions/w2_relres_best.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "W2_RelRes_Sigma",
        "#frac{W^{2}_{#Sigma} - W^{2}_{MC}}{W^{2}_{MC}}",
        "Number of events",
        -0.8, 0.5,
        "figs/inclusive/resolutions/w2_relres_sigma.png",
        "dscb"
    ));

    // =================================================================
    // Binned resolution plots (Q2/x/y)
    // =================================================================
    binned_plot_ptr = new PlotOptionsBinnedRelRes(
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
        "figs/inclusive/resolutions/q2_relres_binned_em.png",
        "figs/inclusive/resolutions/binned/bins/q2_relres_binned_em",
        std::make_pair(5.0, 200),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.2, 0.25, 0.35, 0.4);
    binned_plot_ptr->SetRangeY(-0.15, 0.10);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
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
        "figs/inclusive/resolutions/q2_relres_binned_da.png",
        "figs/inclusive/resolutions/binned/bins/q2_relres_binned_da",
        std::make_pair(5.0, 200),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.2, 0.25, 0.35, 0.4);
    binned_plot_ptr->SetRangeY(-0.15, 0.10);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
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
        "figs/inclusive/resolutions/q2_relres_binned_sigma.png",
        "figs/inclusive/resolutions/binned/bins/q2_relres_binned_sigma",
        std::make_pair(5.0, 200),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.7, 0.75, 0.85, 0.9);
    binned_plot_ptr->SetRangeY(-0.15, 0.10);
    plots.push_back(binned_plot_ptr);

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
        "figs/inclusive/resolutions/xbj_relres_binned_em.png",
        "figs/inclusive/resolutions/binned/bins/xbj_relres_binned_em",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    binned_plot_ptr->SetRangeY(-0.05, 0.05);
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
        "figs/inclusive/resolutions/xbj_relres_binned_da.png",
        "figs/inclusive/resolutions/binned/bins/xbj_relres_binned_da",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    binned_plot_ptr->SetRangeY(-0.05, 0.05);
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
        "figs/inclusive/resolutions/xbj_relres_binned_sigma.png",
        "figs/inclusive/resolutions/binned/bins/xbj_relres_binned_sigma",
        std::make_pair(1e-3, 0.3),
        true
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    binned_plot_ptr->SetRangeY(-0.05, 0.05);
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
        "figs/inclusive/resolutions/y_relres_binned_em.png",
        "figs/inclusive/resolutions/binned/bins/y_relres_binned_em",
        std::make_pair(0.0, 1.0),
        false
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    binned_plot_ptr->SetRangeY(-0.6, 0.15);
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
        "figs/inclusive/resolutions/y_relres_binned_da.png",
        "figs/inclusive/resolutions/binned/bins/y_relres_binned_da",
        std::make_pair(0.0, 1.0),
        false
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    binned_plot_ptr->SetRangeY(-0.6, 0.15);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "y_RelRes_binned_Sigma",
        ";y_{MC};#frac{y_{#Sigma} - y_{MC}}{y_{MC}}",
        "y_{#Sigma}",
        "",
        {
            {-0.6,0.15},{-0.6,0.0},{-0.6,0.0},{-0.6,0.02},{-0.5,0.0},
            {-0.5,0.02},{-0.5,0},{-0.5,0},{-0.3,0},{-0.2,0}
        },
        "figs/inclusive/resolutions/y_relres_binned_sigma.png",
        "figs/inclusive/resolutions/binned/bins/y_relres_binned_sigma",
        std::make_pair(0.0, 1.0),
        false,
        "crystalball"
    );
    binned_plot_ptr->SetLegendPosition(0.15, 0.15, 0.3, 0.3);
    binned_plot_ptr->SetRangeY(-0.6, 0.15);
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
        "figs/inclusive/resolutions/w2_relres_binned_em.png",
        "figs/inclusive/resolutions/binned/bins/w2_relres_binned_em",
        std::make_pair(10.0, 1.0e4),
        true
    );
    binned_plot_ptr->SetRangeY(-0.4, 0.4);
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
        "figs/inclusive/resolutions/w2_relres_binned_da.png",
        "figs/inclusive/resolutions/binned/bins/w2_relres_binned_da",
        std::make_pair(10.0, 1.0e4),
        true
    );
    binned_plot_ptr->SetRangeY(-0.4, 0.4);
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "W2_RelRes_binned_Best",
        "W^{2} relative bin by bin resolution (Best);W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{best} - W^{2}_{MC}}{W^{2}_{MC}}",
        "W^{2}_{best}",
        "",
        {
            {0,0}, {-0.35,0.35}, {-0.2,0.2}, {-0.12,0.12}, {-0.08,0.08},
            {-0.05,0.05}, {-0.05,0.05}
        },
        "figs/inclusive/resolutions/w2_relres_binned_best.png",
        "figs/inclusive/resolutions/binned/bins/w2_relres_binned_best",
        std::make_pair(10.0, 1.0e4),
        true
    );
    binned_plot_ptr->SetRangeY(-0.4, 0.4);
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
        "figs/inclusive/resolutions/w2_relres_binned_sigma.png",
        "figs/inclusive/resolutions/binned/bins/w2_relres_binned_sigma",
        std::make_pair(10.0, 1.0e4),
        true
    );
    binned_plot_ptr->SetRangeY(-0.4, 0.4);
    plots.push_back(binned_plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xL_MC", "xL_B0", "xL_RP"},
        {"MC Truth", "B0 Reco", "RP Reco"},
        {"hist", "pe", "pe"},
        "x_{L} Distributions",
        "x_{L}",
        "Number of events",
        "figs/diffractive/distributions/xL_distributions.png",
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
        "figs/diffractive/distributions/xL_corr_unbinned.png",
        {0.75, 1.05},
        {0.75, 1.05}
    ));
#endif

    // =================================================================
    // Diffractive: x_{pom} (definition) - split by B0/RP
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"xpom_truth_all", "xpom_reco_W2Best_all", "xpom_reco_DA_all", "xpom_reco_Sigma_all"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{pom} Distributions (All)",
        "x_{pom}",
        "Number of events",
        "figs/diffractive/distributions/xpom_def_comparison_all.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-4, 0.3);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_truth_B0", "xpom_reco_W2Best_B0", "xpom_reco_DA_B0", "xpom_reco_Sigma_B0"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{pom} Distributions (B0)",
        "x_{pom}",
        "Number of events",
        "figs/diffractive/distributions/xpom_def_comparison_b0.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-4, 0.3);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"xpom_truth_RP", "xpom_reco_W2Best_RP", "xpom_reco_DA_RP", "xpom_reco_Sigma_RP"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "x_{pom} Distributions (RP)",
        "x_{pom}",
        "Number of events",
        "figs/diffractive/distributions/xpom_def_comparison_rp.png",
        true,
        false
    );
    plot_ptr->SetRangeX(1e-4, 0.3);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_W2Best_B0"},
        {"EM"},
        {kRed},
        {20},
        "x_{pom} Correlation (W^{2}_{best}, B0)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/distributions/xpom_corr_unbinned_em_b0.png",
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
        "figs/diffractive/distributions/xpom_corr_unbinned_da_b0.png",
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
        "figs/diffractive/distributions/xpom_corr_unbinned_sigma_b0.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_xpom_W2Best_RP"},
        {"EM"},
        {kRed},
        {20},
        "x_{pom} Correlation (W^{2}_{best}, RP)",
        "Truth x_{pom}",
        "Reco x_{pom}",
        "figs/diffractive/distributions/xpom_corr_unbinned_em_rp.png",
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
        "figs/diffractive/distributions/xpom_corr_unbinned_da_rp.png",
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
        "figs/diffractive/distributions/xpom_corr_unbinned_sigma_rp.png",
        {1e-4, 0.3},
        {1e-4, 0.3},
        true,
        true
    ));

    // Combined B0+RP response matrices for diffractive variables
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_xpom_W2Best_B0", "Truth x_{pom}", "Reco x_{pom}",
            "figs/diffractive/response/response_xpom_em.png", true, true, {1e-4, 0.3}, {1e-4, 0.3});
        p->SetSecondHistogram("Response_xpom_W2Best_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_xpom_EM_B0", "Truth x_{pom}", "Reco x_{pom}",
            "figs/diffractive/response/response_xpom_w2best.png", true, true, {1e-4, 0.3}, {1e-4, 0.3});
        p->SetSecondHistogram("Response_xpom_EM_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_xpom_DA_B0", "Truth x_{pom}", "Reco x_{pom}",
            "figs/diffractive/response/response_xpom_da.png", true, true, {1e-4, 0.3}, {1e-4, 0.3});
        p->SetSecondHistogram("Response_xpom_DA_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_xpom_Sigma_B0", "Truth x_{pom}", "Reco x_{pom}",
            "figs/diffractive/response/response_xpom_sigma.png", true, true, {1e-4, 0.3}, {1e-4, 0.3});
        p->SetSecondHistogram("Response_xpom_Sigma_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_beta_W2Best_B0", "Truth #beta", "Reco #beta",
            "figs/diffractive/response/response_beta_em.png", false, false, {0.0, 1.0}, {0.0, 1.0});
        p->SetSecondHistogram("Response_beta_W2Best_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_beta_EM_B0", "Truth #beta", "Reco #beta",
            "figs/diffractive/response/response_beta_w2best.png", false, false, {0.0, 1.0}, {0.0, 1.0});
        p->SetSecondHistogram("Response_beta_EM_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_beta_DA_B0", "Truth #beta", "Reco #beta",
            "figs/diffractive/response/response_beta_da.png", false, false, {0.0, 1.0}, {0.0, 1.0});
        p->SetSecondHistogram("Response_beta_DA_RP");
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsResponseMatrix(
            "Response_beta_Sigma_B0", "Truth #beta", "Reco #beta",
            "figs/diffractive/response/response_beta_sigma.png", false, false, {0.0, 1.0}, {0.0, 1.0});
        p->SetSecondHistogram("Response_beta_Sigma_RP");
        plots.push_back(p);
    }

    // =================================================================
    // Diffractive: beta = x_{Bj}/x_{pom} (B0/RP)
    // =================================================================
    plot_ptr = new PlotOptions1D(
        {"beta_truth_all", "beta_reco_W2Best_all", "beta_reco_DA_all", "beta_reco_Sigma_all"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "#beta Distributions (All)",
        "#beta",
        "Number of events",
        "figs/diffractive/distributions/beta_comparison_all.png",
        false,
        false
    );
    plot_ptr->SetRangeX(0.0, 1.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"beta_truth_B0", "beta_reco_W2Best_B0", "beta_reco_DA_B0", "beta_reco_Sigma_B0"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "#beta Distributions (B0)",
        "#beta",
        "Number of events",
        "figs/diffractive/distributions/beta_comparison_b0.png",
        false,
        false
    );
    plot_ptr->SetRangeX(0.0, 1.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plot_ptr = new PlotOptions1D(
        {"beta_truth_RP", "beta_reco_W2Best_RP", "beta_reco_DA_RP", "beta_reco_Sigma_RP"},
        {"MC Truth", "Reco EM", "Reco DA", "Reco Sigma"},
        {"hist", "pe", "pe", "pe"},
        "#beta Distributions (RP)",
        "#beta",
        "Number of events",
        "figs/diffractive/distributions/beta_comparison_rp.png",
        false,
        false
    );
    plot_ptr->SetRangeX(0.0, 1.0);
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_W2Best_B0"},
        {"EM"},
        {kRed},
        {20},
        "#beta Correlation (W^{2}_{best}, B0)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/distributions/beta_corr_unbinned_em_b0.png",
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
        "figs/diffractive/distributions/beta_corr_unbinned_da_b0.png",
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
        "figs/diffractive/distributions/beta_corr_unbinned_sigma_b0.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_beta_W2Best_RP"},
        {"EM"},
        {kRed},
        {20},
        "#beta Correlation (W^{2}_{best}, RP)",
        "#beta_{truth}",
        "#beta_{reco}",
        "figs/diffractive/distributions/beta_corr_unbinned_em_rp.png",
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
        "figs/diffractive/distributions/beta_corr_unbinned_da_rp.png",
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
        "figs/diffractive/distributions/beta_corr_unbinned_sigma_rp.png",
        {0.0, 1.0},
        {0.0, 1.0},
        false,
        false
    ));

    // =================================================================
    // Diffractive: M_X^2 plots
    // =================================================================
    // The MX^2 comparison overlay (hadronic-sum vs kinematic, truth vs reco)
    // is drawn by the dedicated PlotMX2Comparison() function below; see the
    // call site near the top of main(). The PlotOptions1D auto-styling
    // collapses all "*truth*" curves to the same color, so a hand-styled
    // canvas is used to keep the four curves visually distinct.

    plots.push_back(new PlotOptionsCombinedCorrelation(
        {"g_MX2"},
        {"Reco vs Truth"},
        {kBlack},
        {20},
        "M_{X}^{2} Correlation (Unbinned)",
        "M_{X,truth}^{2} [GeV^{2}]",
        "M_{X,reco}^{2} [GeV^{2}]",
        "figs/diffractive/distributions/MX2_corr_unbinned.png",
        {1e-3, 1000.0},
        {1e-3, 1000.0},
        true,
        true
    ));

    plots.push_back(new PlotOptionsRelRes(
        "MX2_RelRes",
        "#frac{M_{X,reco}^{2} - M_{X,truth}^{2}}{M_{X,truth}^{2}}",
        "Number of events",
        -0.6, 0.6,
        "figs/diffractive/resolutions/simple/MX2_resolution.png",
        "double_sided_crystalball"
    ));

#if 0
    plots.push_back(new PlotOptions2D(
        "MX2_corr",
        "M_{X,truth}^{2} [GeV^{2}]",
        "M_{X,reco}^{2} [GeV^{2}]",
        "figs/diffractive/distributions/MX2_truth_vs_reco.png",
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
        "figs/diffractive/distributions/MX2_vs_t_truth.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptions2D(
        "MX2_t_B0",
        "M_{X}^{2} [GeV^{2}]",
        "|t| [GeV^{2}]",
        "figs/diffractive/distributions/MX2_vs_t_b0.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptions2D(
        "MX2_t_RP",
        "M_{X}^{2} [GeV^{2}]",
        "|t| [GeV^{2}]",
        "figs/diffractive/distributions/MX2_vs_t_rp.png",
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
        "Number of events",
        "figs/diffractive/distributions/t_distributions.png",
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
        "Number of events",
        "figs/diffractive/distributions/t_distributions_logy.png",
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
        "figs/diffractive/distributions/dsigma_dt.png",
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
        "Number of events",
        "figs/diffractive/distributions/theta_distributions.png",
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
        "figs/diffractive/distributions/t_corr_unbinned.png",
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
        "Number of events",
        -0.1, 0.2,
        "figs/diffractive/resolutions/t_res_b0.png",
        "dscb"
    ));

    plots.push_back(new PlotOptionsRelRes(
        "t_res_RP",
        "(|t|_{reco} - |t|_{truth})/|t|_{truth}",
        "Number of events",
        -0.3, 0.4,
        "figs/diffractive/resolutions/t_res_rp.png",
        "dscb"
    ));

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "t_RelRes_binned_B0",
        ";|t|_{truth} [GeV^{2}];#frac{|t|_{reco}-|t|_{truth}}{|t|_{truth}}",
        "|t|_{B0}",
        "",
        {},
        "figs/diffractive/resolutions/t_relres_binned_b0.png",
        "figs/diffractive/resolutions/binned/bins/t_relres_binned_b0",
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
        "figs/diffractive/resolutions/t_relres_binned_rp.png",
        "figs/diffractive/resolutions/binned/bins/t_relres_binned_rp",
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
        "Number of events",
        "figs/inclusive/distributions/EPz_distribution.png",
        false,
        false
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    // EPz logY with cut lines is handled by PlotEPzWithCuts() below

    plot_ptr = new PlotOptions1D(
        {"h_eta_max_truth", "h_eta_max"},
        {"MC Truth", "Reconstruction"},
        {"hist", "pe"},
        "Maximum Pseudorapidity per Event",
        "#eta_{max}",
        "Number of events",
        "figs/inclusive/distributions/eta_max_distribution.png",
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
        "Number of events",
        "figs/inclusive/distributions/eta_max_distribution_logY.png",
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
        "Number of events",
        "figs/diffractive/distributions/MX2_distribution.png",
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
        "Number of events",
        "figs/diffractive/distributions/MX2_distribution_logY.png",
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
        "figs/diffractive/distributions/t_pdf_comparison.png",
        true,
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsRelRes(
        "xL_res_B0",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        "Number of events",
        -0.4, 0.3,
        "figs/diffractive/resolutions/simple/xL_resolution_B0.png",
        "dscb"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "xL_res_RP",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        "Number of events",
        -0.3, 0.3,
        "figs/diffractive/resolutions/simple/xL_resolution_RP.png",
        "dscb"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "xpom_res_B0",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        "Number of events",
        -1.0, 1.0,
        "figs/diffractive/resolutions/simple/xpom_resolution_B0.png",
        "dscb"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "xpom_res_RP",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        "Number of events",
        -1.0, 1.0,
        "figs/diffractive/resolutions/simple/xpom_resolution_RP.png",
        "dscb"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "beta_res_B0",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Number of events",
        -0.5, 0.5,
        "figs/diffractive/resolutions/simple/beta_resolution_B0.png",
        "dscb"
    ));
    plots.push_back(new PlotOptionsRelRes(
        "beta_res_RP",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        "Number of events",
        -0.5, 0.5,
        "figs/diffractive/resolutions/simple/beta_resolution_RP.png",
        "dscb"
    ));

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xL_RelRes_binned_B0",
        "B0 x_{L} Resolution vs Truth x_{L}",
        "x_{L,truth}",
        "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
        {{-0.3, 0.3}},
        "figs/diffractive/resolutions/binned/xL_resolution_binned_B0.png",
        "figs/diffractive/resolutions/binned/bins/xL_B0",
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
        "figs/diffractive/resolutions/binned/xL_resolution_binned_RP.png",
        "figs/diffractive/resolutions/binned/bins/xL_RP",
        std::make_pair(0.75, 1.05),
        false
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_B0",
        "B0 x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{0,0},{0,0},{0,0},{0,0},{0,0},{-1.0,1.0}},
        "figs/diffractive/resolutions/binned/xpom_resolution_binned_B0.png",
        "figs/diffractive/resolutions/binned/bins/xpom_B0",
        std::make_pair(1e-4, 0.4),
        true,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_RP",
        "RP x_{pom} Resolution vs Truth x_{pom}",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{0,0},{0,0},{0,0},{-0.8,0.9},{-1.0,1.0},{-1.0,1.0}},
        "figs/diffractive/resolutions/binned/xpom_resolution_binned_RP.png",
        "figs/diffractive/resolutions/binned/bins/xpom_RP",
        std::make_pair(1e-4, 0.4),
        true,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_EM_B0",
        "B0 x_{pom} Resolution vs Truth x_{pom} (W^{2}_{best})",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{0,0},{0,0},{0,0},{0,0},{0,0},{-1.0,1.0}},
        "figs/diffractive/resolutions/binned/xpom_resolution_binned_w2best_B0.png",
        "figs/diffractive/resolutions/binned/bins/xpom_w2best_B0",
        std::make_pair(1e-4, 0.4),
        true,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "xpom_RelRes_binned_EM_RP",
        "RP x_{pom} Resolution vs Truth x_{pom} (W^{2}_{best})",
        "x_{pom,truth}",
        "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
        {{0,0},{0,0},{0,0},{-0.8,0.9},{-1.0,1.0},{-1.0,1.0}},
        "figs/diffractive/resolutions/binned/xpom_resolution_binned_w2best_RP.png",
        "figs/diffractive/resolutions/binned/bins/xpom_w2best_RP",
        std::make_pair(1e-4, 0.4),
        true,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_B0",
        "B0 #beta Resolution vs Truth #beta",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.4, 0.8}, {-0.4,0.8},{-0.3,0.6},{-0.2,0.4},{-0.1,0.2}},
        "figs/diffractive/resolutions/binned/beta_resolution_binned_B0.png",
        "figs/diffractive/resolutions/binned/bins/beta_B0",
        std::make_pair(0.0, 1.0),
        false,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_RP",
        "RP #beta Resolution vs Truth #beta",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.4, 0.8}, {-0.4,0.8},{-0.3,0.6},{-0.2,0.4},{-0.1,0.2}},
        "figs/diffractive/resolutions/binned/beta_resolution_binned_RP.png",
        "figs/diffractive/resolutions/binned/bins/beta_RP",
        std::make_pair(0.0, 1.0),
        false,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_EM_B0",
        "B0 #beta Resolution vs Truth #beta (W^{2}_{best})",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.4, 0.8}, {-0.4,0.8},{-0.3,0.6},{-0.2,0.4},{-0.1,0.2}},
        "figs/diffractive/resolutions/binned/beta_resolution_binned_w2best_B0.png",
        "figs/diffractive/resolutions/binned/bins/beta_w2best_B0",
        std::make_pair(0.0, 1.0),
        false,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "beta_RelRes_binned_EM_RP",
        "RP #beta Resolution vs Truth #beta (W^{2}_{best})",
        "#beta_{truth}",
        "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
        {{-0.4, 0.8}, {-0.4,0.8},{-0.3,0.6},{-0.2,0.4},{-0.1,0.2}},
        "figs/diffractive/resolutions/binned/beta_resolution_binned_w2best_RP.png",
        "figs/diffractive/resolutions/binned/bins/beta_w2best_RP",
        std::make_pair(0.0, 1.0),
        false,
        "dscb"
    );
    plots.push_back(binned_plot_ptr);

    // Combined (stitched) RP/B0 binned resolution plots. Below the boundary
    // points are taken from the RP histogram; above from the B0 histogram.
    // Boundary choices follow detector acceptance (θ ≈ 5 mrad).
    {
        auto* p = new PlotOptionsBinnedRelRes(
            "t_RelRes_binned_RP",
            ";|t|_{truth} [GeV^{2}];#frac{|t|_{reco}-|t|_{truth}}{|t|_{truth}}",
            "|t|",
            "",
            {},
            "figs/diffractive/resolutions/t_relres_binned_combined.png",
            "figs/diffractive/resolutions/binned/bins/t_relres_binned_combined",
            std::make_pair(1e-3, 2.0),
            true);
        p->SetStitchedDetectors("t_RelRes_binned_RP", "t_RelRes_binned_B0", 0.15);
        p->SetLegendPosition(0.18, 0.72, 0.40, 0.88);
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsBinnedRelRes(
            "xpom_RelRes_binned_RP",
            "x_{pom} Resolution vs Truth x_{pom}",
            "x_{pom,truth}",
            "(x_{pom,reco} - x_{pom,truth})/x_{pom,truth}",
            {},
            "figs/diffractive/resolutions/binned/xpom_resolution_binned_combined.png",
            "figs/diffractive/resolutions/binned/bins/xpom_combined",
            std::make_pair(1e-4, 0.4),
            true);
        p->SetStitchedDetectors("xpom_RelRes_binned_RP", "xpom_RelRes_binned_B0", 0.03);
        p->SetLegendPosition(0.18, 0.72, 0.40, 0.88);
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsBinnedRelRes(
            "xL_RelRes_binned_RP",
            "x_{L} Resolution vs Truth x_{L}",
            "x_{L,truth}",
            "(x_{L,reco} - x_{L,truth})/x_{L,truth}",
            {},
            "figs/diffractive/resolutions/binned/xL_resolution_binned_combined.png",
            "figs/diffractive/resolutions/binned/bins/xL_combined",
            std::make_pair(0.75, 1.05),
            false);
        p->SetStitchedDetectors("xL_RelRes_binned_RP", "xL_RelRes_binned_B0", 0.97, "B0", "RP");
        p->SetLegendPosition(0.18, 0.72, 0.40, 0.88);
        plots.push_back(p);
    }
    {
        auto* p = new PlotOptionsBinnedRelRes(
            "beta_RelRes_binned_RP",
            "#beta Resolution vs Truth #beta",
            "#beta_{truth}",
            "(#beta_{reco} - #beta_{truth})/#beta_{truth}",
            {},
            "figs/diffractive/resolutions/binned/beta_resolution_binned_combined.png",
            "figs/diffractive/resolutions/binned/bins/beta_combined",
            std::make_pair(0.0, 1.0),
            false);
        // β has no sharp detector cut in β itself; boundary at β = 0 so all
        // points fall on the "B0" side by default — instead we overlay both
        // for β by pointing RP hist at itself up to β = 0.55 (below that RP
        // covers more of the proton phase space) and B0 above.
        p->SetStitchedDetectors("beta_RelRes_binned_RP", "beta_RelRes_binned_B0", 0.55);
        p->SetLegendPosition(0.18, 0.72, 0.40, 0.88);
        plots.push_back(p);
    }

    binned_plot_ptr = new PlotOptionsBinnedRelRes(
        "MX2_RelRes_binned",
        "M_{X}^{2} Resolution vs Truth M_{X}^{2}",
        "M_{X,truth}^{2} [GeV^{2}]",
        "(M_{X,reco}^{2} - M_{X,truth}^{2})/M_{X,truth}^{2}",
        {{-0.5, 0.5}},
        "figs/diffractive/resolutions/binned/MX2_resolution_binned.png",
        "figs/diffractive/resolutions/binned/bins/MX2",
        std::make_pair(1e-3, 1000.0),
        true
    );
    plots.push_back(binned_plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_EM",
        "Q^{2} (true) [GeV^{2}]",
        "Q^{2} (EM) [GeV^{2}]",
        "figs/inclusive/response/response_matrix_Q2_EM.png",
        true,
        true,
        {1.0, 300.0},
        {1.0, 300.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_DA",
        "Q^{2} (true) [GeV^{2}]",
        "Q^{2} (DA) [GeV^{2}]",
        "figs/inclusive/response/response_matrix_Q2_DA.png",
        true,
        true,
        {1.0, 300.0},
        {1.0, 300.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "Corr_Q2_Sigma",
        "Q^{2} (true) [GeV^{2}]",
        "Q^{2} (Sigma) [GeV^{2}]",
        "figs/inclusive/response/response_matrix_Q2_Esigma.png",
        true,
        true,
        {1.0, 300.0},
        {1.0, 300.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_EM",
        "x_{Bj} (true)",
        "x_{Bj} (EM)",
        "figs/inclusive/response/response_matrix_x_EM.png",
        true,
        true,
        {1e-3, 0.3},
        {1e-3, 0.3}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_DA",
        "x_{Bj} (true)",
        "x_{Bj} (DA)",
        "figs/inclusive/response/response_matrix_x_DA.png",
        true,
        true,
        {1e-3, 0.3},
        {1e-3, 0.3}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "x_Corr_Sigma",
        "x_{Bj} (true)",
        "x_{Bj} (Sigma)",
        "figs/inclusive/response/response_matrix_x_Sigma.png",
        true,
        true,
        {1e-3, 0.3},
        {1e-3, 0.3}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_EM",
        "y (true)",
        "y (EM)",
        "figs/inclusive/response/response_matrix_y_EM.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_DA",
        "y (true)",
        "y (DA)",
        "figs/inclusive/response/response_matrix_y_DA.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "y_Corr_Sigma",
        "y (true)",
        "y (Sigma)",
        "figs/inclusive/response/response_matrix_y_Sigma.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_B0",
        "Truth |t| [GeV^{2}]",
        "B0 Reco |t| [GeV^{2}]",
        "figs/diffractive/response/response_matrix_t_B0.png",
        true,
        true,
        {1e-3, 2.0},
        {1e-3, 2.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "t_corr_RP",
        "Truth |t| [GeV^{2}]",
        "RP Reco |t| [GeV^{2}]",
        "figs/diffractive/response/response_matrix_t_RP.png",
        true,
        true,
        {1e-3, 2.0},
        {1e-3, 2.0}
    ));

    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_B0",
        "Truth x_{L}",
        "B0 Reco x_{L}",
        "figs/diffractive/response/response_matrix_xL_B0.png",
        false,
        false,
        {0.75, 1.05},
        {0.75, 1.05}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xL_corr_RP",
        "Truth x_{L}",
        "RP Reco x_{L}",
        "figs/diffractive/response/response_matrix_xL_RP.png",
        false,
        false,
        {0.75, 1.05},
        {0.75, 1.05}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_B0",
        "Truth x_{pom} (1-x_{L})",
        "B0 Reco x_{pom} (1-x_{L})",
        "figs/diffractive/response/response_matrix_xpom_B0.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_corr_RP",
        "Truth x_{pom} (1-x_{L})",
        "RP Reco x_{pom} (1-x_{L})",
        "figs/diffractive/response/response_matrix_xpom_RP.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_B0",
        "Truth #beta",
        "B0 Reco #beta",
        "figs/diffractive/response/response_matrix_beta_B0.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_corr_RP",
        "Truth #beta",
        "RP Reco #beta",
        "figs/diffractive/response/response_matrix_beta_RP.png",
        false,
        false,
        {0.0, 1.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "MX2_corr",
        "M_{X,truth}^{2} [GeV^{2}]",
        "M_{X,reco}^{2} [GeV^{2}]",
        "figs/diffractive/response/response_matrix_MX2.png",
        true,
        true,
        {1e-3, 1000.0},
        {1e-3, 1000.0}
    ));

    // Combined (B0+RP) response matrices with a dashed RP/B0 boundary line.
    // Boundary values chosen from detector acceptance: θ ≈ 5 mrad.
    // |t|: RP dominates |t| ≲ 0.15 GeV²; xpom/xL: RP dominates xpom ≲ 0.03 (xL ≳ 0.97).
    {
        auto* p_t = new PlotOptionsResponseMatrix(
            "t_corr_RP",
            "Truth |t| [GeV^{2}]",
            "Reco |t| [GeV^{2}]",
            "figs/diffractive/response/response_matrix_t_combined.png",
            true, true, {1e-3, 2.0}, {1e-3, 2.0});
        p_t->SetSecondHistogram("t_corr_B0");
        p_t->SetDetectorBoundary(0.15, "RP", "B0");
        plots.push_back(p_t);

        auto* p_xL = new PlotOptionsResponseMatrix(
            "xL_corr_RP",
            "Truth x_{L}",
            "Reco x_{L}",
            "figs/diffractive/response/response_matrix_xL_combined.png",
            false, false, {0.75, 1.05}, {0.75, 1.05});
        p_xL->SetSecondHistogram("xL_corr_B0");
        p_xL->SetDetectorBoundary(0.97, "B0", "RP");
        plots.push_back(p_xL);

        auto* p_xpom = new PlotOptionsResponseMatrix(
            "xpom_corr_RP",
            "Truth x_{pom} (1-x_{L})",
            "Reco x_{pom} (1-x_{L})",
            "figs/diffractive/response/response_matrix_xpom_combined.png",
            true, true, {1e-4, 0.4}, {1e-4, 0.4});
        p_xpom->SetSecondHistogram("xpom_corr_B0");
        p_xpom->SetDetectorBoundary(0.03, "RP", "B0");
        plots.push_back(p_xpom);

        auto* p_beta = new PlotOptionsResponseMatrix(
            "beta_corr_RP",
            "Truth #beta",
            "Reco #beta",
            "figs/diffractive/response/response_matrix_beta_combined.png",
            false, false, {0.0, 1.0}, {0.0, 1.0});
        p_beta->SetSecondHistogram("beta_corr_B0");
        plots.push_back(p_beta);
    }

    plot_ptr = new PlotOptions1D(
        {"xpom_MC", "xpom_def_MC"},
        {"x_{pom} = 1 - x_{L}", "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})"},
        {"hist", "hist"},
        "MC Truth x_{pom} Comparison",
        "x_{pom}",
        "Number of events",
        "figs/diffractive/distributions/xpom_comparison_MC_logxy.png",
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
        "Number of events",
        "figs/diffractive/distributions/xpom_comparison_B0_logxy.png",
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
        "Number of events",
        "figs/diffractive/distributions/xpom_comparison_RP_logxy.png",
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
        "Number of events",
        "figs/diffractive/distributions/xpom_comparison_all_logxy.png",
        true,
        true
    );
    plot_ptr->SetLegendPosition(0.1, 0.7, 0.45, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_MC",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/diffractive/distributions/xpom_2D_comparison_MC.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_B0",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/diffractive/distributions/xpom_2D_comparison_B0.png",
        true,
        true,
        {1e-4, 0.4},
        {1e-4, 0.4}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "xpom_comp_RP",
        "x_{pom} = 1 - x_{L}",
        "x_{pom} = (M_{X}^{2}+Q^{2}+|t|)/(W^{2}+Q^{2}-m_{p}^{2})",
        "figs/diffractive/distributions/xpom_2D_comparison_RP.png",
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
        "Number of events",
        "figs/diffractive/distributions/theta_comparison_B0_acceptance.png",
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
        "Number of events",
        "figs/diffractive/distributions/theta_comparison_B0_acceptance_logxy.png",
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
        "Number of events",
        "figs/diffractive/distributions/beta_distributions_logy.png",
        false,
        true
    );
    plot_ptr->SetLegendPosition(0.6, 0.7, 0.85, 0.9);
    plots.push_back(plot_ptr);

    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_Q2",
        "Q^{2} [GeV^{2}]",
        "#beta",
        "figs/diffractive/distributions/beta_vs_Q2.png",
        true,
        false,
        {1.0, 300.0},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_xpom",
        "x_{pom}",
        "#beta",
        "figs/diffractive/distributions/beta_vs_xpom.png",
        true,
        false,
        {1e-4, 0.4},
        {0.0, 1.0}
    ));
    plots.push_back(new PlotOptionsResponseMatrix(
        "beta_vs_t",
        "|t| [GeV^{2}]",
        "#beta",
        "figs/diffractive/distributions/beta_vs_t.png",
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
    // Standalone plot functions (Tasks 1, 2, 3B, 5)
    // =================================================================
    PlotEPzWithCuts(inputFile);
    // Resolution comparison plots for all inclusive kinematic quantities
    PlotResolutionComparison(inputFile,
        {"W2_RelRes_binned_EM", "W2_RelRes_binned_DA", "W2_RelRes_binned_Best", "W2_RelRes_binned_Sigma", "W2_RelRes_binned_ESigma"},
        {"EM", "DA", "Best", "Sigma", "ESigma"},
        "W^{2}_{MC} [GeV^{2}]", "Relative resolution #mu #pm #sigma",
        10.0, 1.0e4, -0.4, 0.4, true,
        "figs/inclusive/resolutions/w2_relres_binned_comparison.png");

    PlotResolutionComparison(inputFile,
        {"Q2_RelRes_binned_EM", "Q2_RelRes_binned_DA", "Q2_RelRes_binned_Sigma"},
        {"EM", "DA", "Sigma"},
        "Q^{2}_{MC} [GeV^{2}]", "Relative resolution #mu #pm #sigma",
        5.0, 200.0, -0.15, 0.15, true,
        "figs/inclusive/resolutions/q2_relres_binned_comparison.png");

    PlotResolutionComparison(inputFile,
        {"x_RelRes_binned_EM", "x_RelRes_binned_DA", "x_RelRes_binned_Sigma"},
        {"EM", "DA", "Sigma"},
        "x_{Bj,MC}", "Relative resolution #mu #pm #sigma",
        1e-3, 0.3, -0.3, 0.3, true,
        "figs/inclusive/resolutions/xbj_relres_binned_comparison.png");

    PlotResolutionComparison(inputFile,
        {"y_RelRes_binned_EM", "y_RelRes_binned_DA", "y_RelRes_binned_Sigma"},
        {"EM", "DA", "Sigma"},
        "y_{MC}", "Relative resolution #mu #pm #sigma",
        0.0, 1.0, -0.5, 0.5, false,
        "figs/inclusive/resolutions/y_relres_binned_comparison.png");
    Plot3DResponseMatrix(inputFile);
    PlotXQ2PhaseSpaceWithLines(inputFile, disSGeV2);

    // =================================================================
    // Relative resolution vs global bin index (k)
    // =================================================================
    PlotRelResVsK(inputFile,
                  {"Q2_RelRes_vs_k_EM", "Q2_RelRes_vs_k_DA", "Q2_RelRes_vs_k_Sigma"},
                  {"EM", "DA", "Sigma"},
                  "Q^{2} relative resolution vs k",
                  "figs/inclusive/resolutions/binned/q2_relres_vs_k.png");

    PlotRelResVsK(inputFile,
                  {"xpom_RelRes_vs_k_W2Best_B0", "xpom_RelRes_vs_k_DA_B0", "xpom_RelRes_vs_k_Sigma_B0"},
                  {"EM", "DA", "Sigma"},
                  "x_{pom} relative resolution vs k (B0)",
                  "figs/diffractive/resolutions/binned/xpom_relres_vs_k_b0.png");
    PlotRelResVsK(inputFile,
                  {"xpom_RelRes_vs_k_W2Best_RP", "xpom_RelRes_vs_k_DA_RP", "xpom_RelRes_vs_k_Sigma_RP"},
                  {"EM", "DA", "Sigma"},
                  "x_{pom} relative resolution vs k (RP)",
                  "figs/diffractive/resolutions/binned/xpom_relres_vs_k_rp.png");

    PlotRelResVsK(inputFile,
                  {"beta_RelRes_vs_k_W2Best_B0", "beta_RelRes_vs_k_DA_B0", "beta_RelRes_vs_k_Sigma_B0"},
                  {"EM", "DA", "Sigma"},
                  "#beta relative resolution vs k (B0)",
                  "figs/diffractive/resolutions/binned/beta_relres_vs_k_b0.png");
    PlotRelResVsK(inputFile,
                  {"beta_RelRes_vs_k_W2Best_RP", "beta_RelRes_vs_k_DA_RP", "beta_RelRes_vs_k_Sigma_RP"},
                  {"EM", "DA", "Sigma"},
                  "#beta relative resolution vs k (RP)",
                  "figs/diffractive/resolutions/binned/beta_relres_vs_k_rp.png");

    // Combined RP+B0 overlays (one method — EM/W2Best — per variable, both
    // detectors on the same axis, distinguished by marker style & colour).
    PlotRelResVsK(inputFile,
                  {"xpom_RelRes_vs_k_W2Best_RP", "xpom_RelRes_vs_k_W2Best_B0"},
                  {"RP", "B0"},
                  "x_{pom} relative resolution vs k (RP + B0)",
                  "figs/diffractive/resolutions/binned/xpom_relres_vs_k_combined.png");
    PlotRelResVsK(inputFile,
                  {"beta_RelRes_vs_k_W2Best_RP", "beta_RelRes_vs_k_W2Best_B0"},
                  {"RP", "B0"},
                  "#beta relative resolution vs k (RP + B0)",
                  "figs/diffractive/resolutions/binned/beta_relres_vs_k_combined.png");

    // |t|-bin performance metrics
    PlotTBinMetric(inputFile,
                   TBinMetricKind::Efficiency,
                   "|t| Bin Efficiency (Set A)",
                   "Efficiency",
                   "figs/diffractive/performance/t_efficiency_binned.png");
    PlotTBinMetric(inputFile,
                   TBinMetricKind::Acceptance,
                   "|t| Bin Acceptance (Set A)",
                   "Acceptance",
                   "figs/diffractive/performance/t_acceptance_binned.png");
    PlotTBinMetric(inputFile,
                   TBinMetricKind::Purity,
                   "|t| Bin Purity (Set A)",
                   "Purity",
                   "figs/diffractive/performance/t_purity_binned.png");

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
                double markerSize = 0.3 + 30.0 * rms;
                if (markerSize > 4.0) markerSize = 4.0;
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

    // Stitched (B0+RP) circle plot: draws markers from RP on the low side of
    // the x-axis boundary and from B0 on the high side, with a dashed
    // vertical line at the detector boundary.
    auto createCirclePlotStitched = [&](const char* histNameRP,
                                        const char* histNameB0,
                                        const char* saveName,
                                        const char* xTitle, const char* yTitle,
                                        double boundary,
                                        bool boundaryRPIsLow = true,
                                        bool logX = true, bool logY = true) {
        TProfile2D* profRP = (TProfile2D*)inputFile->Get(histNameRP);
        TProfile2D* profB0 = (TProfile2D*)inputFile->Get(histNameB0);
        if (!profRP || !profB0) return;

        TCanvas* c = new TCanvas(Form("c_circle_stitched_%s", histNameRP), "Stitched Circle", 1200, 900);
        c->SetLogx(logX);
        c->SetLogy(logY);
        c->SetGrid();
        c->SetRightMargin(0.15);

        TH2D* frame = new TH2D(Form("frame_stitched_%s", histNameRP), Form(";%s;%s", xTitle, yTitle),
                               profRP->GetNbinsX(), profRP->GetXaxis()->GetXmin(), profRP->GetXaxis()->GetXmax(),
                               profRP->GetNbinsY(), profRP->GetYaxis()->GetXmin(), profRP->GetYaxis()->GetXmax());
        frame->Draw();

        auto drawFromProf = [&](TProfile2D* prof, bool useLowSide) {
            for (int ix = 1; ix <= prof->GetNbinsX(); ix++) {
                double xc = prof->GetXaxis()->GetBinCenter(ix);
                const bool onLowSide = (xc < boundary);
                if (useLowSide != onLowSide) continue;
                for (int iy = 1; iy <= prof->GetNbinsY(); iy++) {
                    double entries = prof->GetBinEntries(prof->GetBin(ix, iy));
                    if (entries < 10) continue;
                    double yc = prof->GetYaxis()->GetBinCenter(iy);
                    double rms = prof->GetBinError(ix, iy);
                    double markerSize = 0.3 + 30.0 * rms;
                    if (markerSize > 4.0) markerSize = 4.0;
                    const bool isRP = (useLowSide == boundaryRPIsLow);
                    TMarker* marker = new TMarker(xc, yc, isRP ? 20 : 21);
                    marker->SetMarkerSize(markerSize);
                    marker->SetMarkerColor(isRP ? kBlue + 1 : kRed + 1);
                    marker->Draw();
                }
            }
        };
        drawFromProf(boundaryRPIsLow ? profRP : profB0, true);
        drawFromProf(boundaryRPIsLow ? profB0 : profRP, false);

        double ymin = profRP->GetYaxis()->GetXmin();
        double ymax = profRP->GetYaxis()->GetXmax();
        TLine* vline = new TLine(boundary, ymin, boundary, ymax);
        vline->SetLineColor(kBlack);
        vline->SetLineStyle(7);
        vline->SetLineWidth(2);
        vline->Draw("SAME");

        TLegend* leg = new TLegend(0.18, 0.82, 0.44, 0.92);
        TMarker mRP(0, 0, 20); mRP.SetMarkerColor(kBlue + 1);
        TMarker mB0(0, 0, 21); mB0.SetMarkerColor(kRed + 1);
        leg->AddEntry(&mRP, "RP", "p");
        leg->AddEntry(&mB0, "B0", "p");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        SaveCanvas(c, saveName);
        delete leg;
        delete vline;
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

    createCirclePlot("Q2_RelRes_vs_xy_EM", "figs/inclusive/resolutions/2d_maps/Q2_RelRes_Q2x_EM.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_DA", "figs/inclusive/resolutions/2d_maps/Q2_RelRes_Q2x_DA.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("Q2_RelRes_vs_xy_Sigma", "figs/inclusive/resolutions/2d_maps/Q2_RelRes_Q2x_Sigma.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_EM", "figs/inclusive/resolutions/2d_maps/x_RelRes_xQ2_EM.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_DA", "figs/inclusive/resolutions/2d_maps/x_RelRes_xQ2_DA.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("x_RelRes_vs_xQ2_Sigma", "figs/inclusive/resolutions/2d_maps/x_RelRes_xQ2_Sigma.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_EM", "figs/inclusive/resolutions/2d_maps/y_RelRes_xQ2_EM.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_DA", "figs/inclusive/resolutions/2d_maps/y_RelRes_xQ2_DA.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("y_RelRes_vs_xQ2_Sigma", "figs/inclusive/resolutions/2d_maps/y_RelRes_xQ2_Sigma.png", "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("t_RelRes_vs_xpomQ2_B0", "figs/diffractive/resolutions/2d_maps/t_RelRes_xpomQ2_B0.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("t_RelRes_vs_xpomQ2_RP", "figs/diffractive/resolutions/2d_maps/t_RelRes_xpomQ2_RP.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("xpom_RelRes_vs_xpomQ2_B0", "figs/diffractive/resolutions/2d_maps/xpom_RelRes_xpomQ2_B0.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("xpom_RelRes_vs_xpomQ2_RP", "figs/diffractive/resolutions/2d_maps/xpom_RelRes_xpomQ2_RP.png", "x_{pom}", "Q^{2} [GeV^{2}]", true, true);
    createCirclePlot("beta_RelRes_vs_betaQ2_B0", "figs/diffractive/resolutions/2d_maps/beta_RelRes_betaQ2_B0.png", "#beta", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("beta_RelRes_vs_betaQ2_RP", "figs/diffractive/resolutions/2d_maps/beta_RelRes_betaQ2_RP.png", "#beta", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("xL_RelRes_vs_xLQ2_B0", "figs/diffractive/resolutions/2d_maps/xL_RelRes_xLQ2_B0.png", "x_{L}", "Q^{2} [GeV^{2}]", false, true);
    createCirclePlot("xL_RelRes_vs_xLQ2_RP", "figs/diffractive/resolutions/2d_maps/xL_RelRes_xLQ2_RP.png", "x_{L}", "Q^{2} [GeV^{2}]", false, true);

    // Combined (RP+B0) circle plots with dashed detector boundary.
    createCirclePlotStitched("t_RelRes_vs_xpomQ2_RP", "t_RelRes_vs_xpomQ2_B0",
                             "figs/diffractive/resolutions/2d_maps/t_RelRes_xpomQ2_combined.png",
                             "x_{pom}", "Q^{2} [GeV^{2}]", 0.03, /*RPlow*/true, true, true);
    createCirclePlotStitched("xpom_RelRes_vs_xpomQ2_RP", "xpom_RelRes_vs_xpomQ2_B0",
                             "figs/diffractive/resolutions/2d_maps/xpom_RelRes_xpomQ2_combined.png",
                             "x_{pom}", "Q^{2} [GeV^{2}]", 0.03, /*RPlow*/true, true, true);
    createCirclePlotStitched("beta_RelRes_vs_betaQ2_RP", "beta_RelRes_vs_betaQ2_B0",
                             "figs/diffractive/resolutions/2d_maps/beta_RelRes_betaQ2_combined.png",
                             "#beta", "Q^{2} [GeV^{2}]", 0.55, /*RPlow*/true, false, true);
    createCirclePlotStitched("xL_RelRes_vs_xLQ2_B0", "xL_RelRes_vs_xLQ2_RP",
                             "figs/diffractive/resolutions/2d_maps/xL_RelRes_xLQ2_combined.png",
                             "x_{L}", "Q^{2} [GeV^{2}]", 0.97, /*RPlow*/false, false, true);

    createBestMethodPlot("Q2_RelRes_vs_xy_EM", "Q2_RelRes_vs_xy_DA", "Q2_RelRes_vs_xy_Sigma",
                         "figs/inclusive/resolutions/2d_maps/Q2_RelRes_Q2x_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("x_RelRes_vs_xQ2_EM", "x_RelRes_vs_xQ2_DA", "x_RelRes_vs_xQ2_Sigma",
                         "figs/inclusive/resolutions/2d_maps/x_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);
    createBestMethodPlot("y_RelRes_vs_xQ2_EM", "y_RelRes_vs_xQ2_DA", "y_RelRes_vs_xQ2_Sigma",
                         "figs/inclusive/resolutions/2d_maps/y_RelRes_xQ2_BestMethod.png",
                         "x_{Bj}", "Q^{2} [GeV^{2}]", true, true);

    createCirclePlot("MX2_RelRes_vs_MX2Q2", "figs/diffractive/resolutions/2d_maps/MX2_RelRes_MX2Q2.png",
                     "M_{X}^{2} [GeV^{2}]", "Q^{2} [GeV^{2}]", true, true);

    if (auto* hEPz2D = (TH2D*)inputFile->Get("h_EPz_2D")) {
        TCanvas* cEPz = new TCanvas("c_EPz_2D", "E-pz Correlation", 1200, 900);
        cEPz->SetRightMargin(0.15);
        cEPz->SetGrid();
        hEPz2D->Draw("COLZ");
        SaveCanvas(cEPz, "figs/inclusive/distributions/EPz_2D.png");
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

    // --- Geometric acceptance correction A(t) ---
    TH1D* hTruthFull = (TH1D*)inputFile->Get("t_truth_mc_full");
    TH1D* hTruthB0   = (TH1D*)inputFile->Get("t_truth_mc_B0");
    TH1D* hTruthRP   = (TH1D*)inputFile->Get("t_truth_mc_RP");

    TH1D* hDSigSumCorr = nullptr;
    if (hDSigSum && hTruthFull && hTruthB0 && hTruthRP) {
        hDSigSumCorr = (TH1D*)hDSigSum->Clone("dsigma_dt_Sum_geomcorr");
        hDSigSumCorr->SetDirectory(nullptr);

        for (int i = 1; i <= hDSigSum->GetNbinsX(); ++i) {
            const double Nfull = hTruthFull->GetBinContent(i);
            const double NB0   = hTruthB0->GetBinContent(i);
            const double NRP   = hTruthRP->GetBinContent(i);
            const double Nacc  = NB0 + NRP;

            if (Nfull <= 0.0 || Nacc <= 0.0) continue;

            const double p     = Nacc / Nfull;
            const double A     = 1.0 / p;
            const double sig_p = std::sqrt(p * (1.0 - p) / Nfull);
            const double sig_A = sig_p / (p * p);

            const double y  = hDSigSum->GetBinContent(i);
            const double dy = hDSigSum->GetBinError(i);
            hDSigSumCorr->SetBinContent(i, y * A);
            hDSigSumCorr->SetBinError(i,
                std::sqrt((dy * A) * (dy * A) + (y * sig_A) * (y * sig_A)));
        }
    }

    // --- Exponential fit to the geom-corrected pseudo-data ---
    TF1* fitExp = nullptr;
    double bValue = 0.0;
    double bError = 0.0;
    TH1D* hDSigDraw = hDSigSumCorr ? hDSigSumCorr : hDSigSum;  // prefer corrected
    if (hDSigDraw) {
        fitExp = new TF1("fit_exp_plot_final", "[0]*TMath::Exp(-[1]*x)", 0.01, 2.0);
        fitExp->SetParameters(1000.0, 5.0);
        fitExp->SetLineColor(kRed);
        fitExp->SetLineWidth(3);
        fitExp->SetLineStyle(2);
        hDSigDraw->Fit(fitExp, "RSQ");
        bValue = fitExp->GetParameter(1);
        bError = fitExp->GetParError(1);
    }

    if (hDSigMC && hDSigDraw && fitExp) {
        // --- Linear canvas ---
        TCanvas* cLin = new TCanvas("c_dsigma_with_fit_linear", "dSigma/dt", 1200, 900);
        cLin->SetLogx();
        cLin->SetGrid();
        hDSigMC->SetLineColor(kBlack);
        hDSigMC->SetLineWidth(2);
        hDSigMC->SetFillColor(kYellow - 9);
        hDSigMC->SetFillStyle(1001);
        hDSigMC->Draw("HIST");
        hDSigDraw->SetMarkerStyle(20); hDSigDraw->SetMarkerColor(kBlue); hDSigDraw->SetLineColor(kBlue);
        hDSigDraw->Draw("PE SAME");
        fitExp->Draw("SAME");
        TLegend leg1(0.56, 0.6, 0.87, 0.9);
        leg1.SetBorderSize(0); leg1.SetFillStyle(0);
        leg1.AddEntry(hDSigMC, "MC", "lf");
        leg1.AddEntry(hDSigDraw, "FF det.s", "pe");
        leg1.AddEntry(fitExp, Form("e^{-b|t|}, b = %.2f #pm %.2f GeV^{-2}", bValue, bError), "l");
        leg1.Draw();
        SaveCanvas(cLin, "figs/cross_sections/dsigma_dt_with_fit.png");
        delete cLin;

        // --- Log-y canvas ---
        TCanvas* cLog = new TCanvas("c_dsigma_with_fit_logy", "dSigma/dt", 1200, 900);
        cLog->SetLogx();
        cLog->SetLogy();
        cLog->SetGrid();
        hDSigMC->Draw("HIST");
        hDSigDraw->Draw("PE SAME");
        fitExp->Draw("SAME");
        TLegend leg2(0.56, 0.6, 0.87, 0.9);
        leg2.SetBorderSize(0); leg2.SetFillStyle(0);
        leg2.AddEntry(hDSigMC, "MC", "lf");
        leg2.AddEntry(hDSigDraw, "FF det.s", "pe");
        leg2.AddEntry(fitExp, Form("e^{-b|t|}, b = %.2f #pm %.2f GeV^{-2}", bValue, bError), "l");
        leg2.Draw();
        SaveCanvas(cLog, "figs/cross_sections/dsigma_dt_logy_with_fit.png");
        delete cLog;
    }
    delete fitExp;
    delete hDSigSumCorr;
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
                pMC->SetFillColor(kYellow - 9);
                pMC->SetFillStyle(1001);
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
                pMC->SetFillColor(kYellow - 9);
                pMC->SetFillStyle(1001);
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

        // H1/ZEUS-like stacked layout in reduced-cross-section form:
        // fixed x_{pom}, Q^{2} on x-axis, and vertically offset #beta series by 3^{i}.
        if (disSGeV2 > 0.0) {
            TH3D* hD3Data = hD3Sum ? hD3Sum : hD3B0;
            const char* dataLegendLabel = "Pseudodata";
            constexpr double kAlphaEM = 1.0 / 137.035999084;
            constexpr double kOffsetBase = 3.0;
            constexpr double kMaxBandRelErr = 0.80;
            struct LegendBox {
                double x1;
                double y1;
                double x2;
                double y2;
            };
            const LegendBox kDefaultReducedLegendBox{0.52, 0.78, 0.90, 0.92};
            const LegendBox kDefaultRawLegendBox{0.52, 0.78, 0.90, 0.92};
            // Keys are xpom bin-center tags produced by makeXpomTag (e.g. "0p002").
            // Low-xpom tags go top-right; high-xpom tags go bottom-left.
            const std::map<std::string, LegendBox> kReducedLegendOverrides = {
                {"0p002",  {0.52, 0.78, 0.90, 0.88}},
                {"0p0065", {0.52, 0.78, 0.90, 0.88}},
                {"0p015",  {0.54, 0.13, 0.92, 0.23}},
                {"0p025",  {0.54, 0.13, 0.92, 0.23}},
                {"0p045",  {0.54, 0.13, 0.92, 0.23}},
                {"0p08",   {0.54, 0.13, 0.92, 0.23}}
            };
            const std::map<std::string, LegendBox> kRawLegendOverrides = {
                {"0p002",  {0.52, 0.78, 0.90, 0.88}},
                {"0p0065", {0.52, 0.78, 0.90, 0.88}},
                {"0p015",  {0.17, 0.13, 0.55, 0.23}},
                {"0p025",  {0.17, 0.13, 0.55, 0.23}},
                {"0p045",  {0.17, 0.13, 0.55, 0.23}},
                {"0p08",   {0.17, 0.13, 0.55, 0.23}}
            };
            auto getLegendBox = [](const std::map<std::string, LegendBox>& overrides,
                                   const LegendBox& defaultBox,
                                   const std::string& key) {
                const auto it = overrides.find(key);
                return (it != overrides.end()) ? it->second : defaultBox;
            };

            auto makeXpomTag = [](double xpomCenter) {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(4) << xpomCenter;
                std::string s = oss.str();
                while (!s.empty() && s.back() == '0') s.pop_back();
                if (!s.empty() && s.back() == '.') s.pop_back();
                for (char& c : s) {
                    if (c == '.') c = 'p';
                }
                return s;
            };

            std::vector<double> q2Edges;
            q2Edges.reserve(nQ2 + 1);
            for (int iQ2 = 1; iQ2 <= nQ2; ++iQ2) {
                q2Edges.push_back(hD3MC->GetXaxis()->GetBinLowEdge(iQ2));
            }
            q2Edges.push_back(hD3MC->GetXaxis()->GetBinUpEdge(nQ2));
            gSystem->mkdir("figs/cross_sections/debug", true);

            for (int iXpom = 1; iXpom <= nXpom; ++iXpom) {
                const double xpomLow = hD3MC->GetZaxis()->GetBinLowEdge(iXpom);
                const double xpomHigh = hD3MC->GetZaxis()->GetBinUpEdge(iXpom);
                const double xpomCenter = hD3MC->GetZaxis()->GetBinCenter(iXpom);
                if (!std::isfinite(xpomCenter) || xpomCenter <= 0.0) continue;

                struct SeriesRow {
                    double betaLow = -1.0;
                    double betaHigh = -1.0;
                    double betaCenter = -1.0;
                    std::vector<double> q2;
                    std::vector<double> theory;
                    std::vector<double> theoryErr;
                    std::vector<double> data;
                    std::vector<double> dataErr;
                };
                struct DebugPoint {
                    int iQ2 = -1;
                    int iBeta = -1;
                    double q2 = -1.0;
                    double beta = -1.0;
                    double y = -1.0;
                    double d3Theory = -1.0;
                    double d3TheoryErr = 0.0;
                    double d3Data = -1.0;
                    double d3DataErr = 0.0;
                    double redTheory = -1.0;
                    double redTheoryErr = 0.0;
                    double redData = -1.0;
                    double redDataErr = 0.0;
                };

                std::vector<SeriesRow> rows;
                rows.reserve(nBeta);
                std::vector<DebugPoint> debugPoints;
                debugPoints.reserve(static_cast<size_t>(nBeta * nQ2));
                for (int iBeta = 1; iBeta <= nBeta; ++iBeta) {
                    SeriesRow row;
                    row.betaLow = hD3MC->GetYaxis()->GetBinLowEdge(iBeta);
                    row.betaHigh = hD3MC->GetYaxis()->GetBinUpEdge(iBeta);
                    row.betaCenter = hD3MC->GetYaxis()->GetBinCenter(iBeta);
                    if (!std::isfinite(row.betaCenter) || row.betaCenter <= 0.0) continue;

                    for (int iQ2 = 1; iQ2 <= nQ2; ++iQ2) {
                        const double q2 = hD3MC->GetXaxis()->GetBinCenter(iQ2);
                        if (!std::isfinite(q2) || q2 <= 0.0) continue;

                        const double y = q2 / (disSGeV2 * xpomCenter * row.betaCenter);
                        if (!std::isfinite(y) || !(y > yMinCutFromFile && y < yMaxCutFromFile)) continue;

                        const double yTerm = 1.0 + (1.0 - y) * (1.0 - y);
                        if (!(yTerm > 0.0)) continue;
                        const double redFactor = row.betaCenter * q2 * q2 / (2.0 * TMath::Pi() * kAlphaEM * kAlphaEM * yTerm);

                        const double th = hD3MC->GetBinContent(iQ2, iBeta, iXpom);
                        const double thErr = hD3MC->GetBinError(iQ2, iBeta, iXpom);
                        const double dt = hD3Data->GetBinContent(iQ2, iBeta, iXpom);
                        const double dtErr = hD3Data->GetBinError(iQ2, iBeta, iXpom);
                        const bool hasTheory = std::isfinite(th) && th > 0.0;
                        const bool hasData = std::isfinite(dt) && dt > 0.0;
                        if (!hasTheory && !hasData) continue;

                        const double thRed = hasTheory ? (xpomCenter * th * redFactor) : -1.0;
                        const double thRedErr = hasTheory ? (xpomCenter * std::max(0.0, thErr) * redFactor) : 0.0;
                        const double dtRed = hasData ? (xpomCenter * dt * redFactor) : -1.0;
                        const double dtRedErr = hasData ? (xpomCenter * std::max(0.0, dtErr) * redFactor) : 0.0;

                        row.q2.push_back(q2);
                        row.theory.push_back(thRed);
                        row.theoryErr.push_back(thRedErr);
                        row.data.push_back(dtRed);
                        row.dataErr.push_back(dtRedErr);

                        DebugPoint dp;
                        dp.iQ2 = iQ2;
                        dp.iBeta = iBeta;
                        dp.q2 = q2;
                        dp.beta = row.betaCenter;
                        dp.y = y;
                        dp.d3Theory = hasTheory ? th : -1.0;
                        dp.d3TheoryErr = hasTheory ? std::max(0.0, thErr) : 0.0;
                        dp.d3Data = hasData ? dt : -1.0;
                        dp.d3DataErr = hasData ? std::max(0.0, dtErr) : 0.0;
                        dp.redTheory = thRed;
                        dp.redTheoryErr = thRedErr;
                        dp.redData = dtRed;
                        dp.redDataErr = dtRedErr;
                        debugPoints.push_back(dp);
                    }
                    if (row.q2.size() >= 2) rows.push_back(std::move(row));
                }
                if (rows.empty()) continue;

                std::sort(rows.begin(), rows.end(), [](const SeriesRow& a, const SeriesRow& b) {
                    return a.betaCenter > b.betaCenter;
                });

                double yMin = std::numeric_limits<double>::max();
                double yMax = 0.0;
                for (size_t iRow = 0; iRow < rows.size(); ++iRow) {
                    const double offset = std::pow(kOffsetBase, static_cast<double>(iRow));
                    for (size_t i = 0; i < rows[iRow].q2.size(); ++i) {
                        if (rows[iRow].theory[i] > 0.0) {
                            const double yVal = rows[iRow].theory[i] * offset;
                            yMin = std::min(yMin, yVal);
                            yMax = std::max(yMax, yVal);
                        }
                        if (rows[iRow].data[i] > 0.0) {
                            const double yVal = rows[iRow].data[i] * offset;
                            yMin = std::min(yMin, yVal);
                            yMax = std::max(yMax, yVal);
                        }
                    }
                }
                if (!(yMax > 0.0) || !(yMin > 0.0) || !std::isfinite(yMin) || !std::isfinite(yMax)) {
                    continue;
                }

                const std::string xpomTag = makeXpomTag(xpomCenter);
                const std::string saveName = "figs/cross_sections/sigma_r_q2_stacked_xpom_" + xpomTag + ".png";
                const std::string legacySaveName = "figs/cross_sections/d3sigma_q2_stacked_xpom_" + xpomTag + ".png";
                TCanvas* cStack = new TCanvas(Form("c_sigma_r_q2_stacked_%d", iXpom),
                                              "reduced cross section stacked layout", 1100, 1400);
                cStack->SetLeftMargin(0.16);
                cStack->SetBottomMargin(0.12);
                cStack->SetLogx();
                cStack->SetLogy();
                cStack->SetGridx();
                cStack->SetGridy();

                TH1D* frame = new TH1D(Form("h_sigma_r_stack_frame_%d", iXpom),
                                       "",
                                       nQ2, q2Edges.data());
                frame->SetDirectory(nullptr);
                frame->SetStats(false);
                frame->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
                frame->GetYaxis()->SetTitle("3^{i} #times x_{pom} #sigma_{r}^{D(3)}");
                frame->GetXaxis()->SetTitleSize(0.048);
                frame->GetYaxis()->SetTitleSize(0.048);
                frame->GetXaxis()->SetLabelSize(0.040);
                frame->GetYaxis()->SetLabelSize(0.040);
                frame->GetXaxis()->SetTitleOffset(1.05);
                frame->GetYaxis()->SetTitleOffset(1.20);
                frame->SetMinimum(std::max(yMin / 2.0, 1e-12));
                frame->SetMaximum(yMax * 1.35);
                frame->Draw("AXIS");

                std::vector<TObject*> transientObjects;
                transientObjects.reserve(rows.size() * 3);
                TGraph* legendTheory = nullptr;
                TGraphErrors* legendData = nullptr;

                for (size_t iRow = 0; iRow < rows.size(); ++iRow) {
                    const double offset = std::pow(kOffsetBase, static_cast<double>(iRow));
                    std::vector<double> xTheory, yTheory, yTheoryErr;
                    std::vector<double> xTheoryBand, yTheoryBand, yTheoryBandErr;
                    std::vector<double> xData, yData, yDataErr;
                    xTheory.reserve(rows[iRow].q2.size());
                    yTheory.reserve(rows[iRow].q2.size());
                    yTheoryErr.reserve(rows[iRow].q2.size());
                    xTheoryBand.reserve(rows[iRow].q2.size());
                    yTheoryBand.reserve(rows[iRow].q2.size());
                    yTheoryBandErr.reserve(rows[iRow].q2.size());
                    xData.reserve(rows[iRow].q2.size());
                    yData.reserve(rows[iRow].q2.size());
                    yDataErr.reserve(rows[iRow].q2.size());

                    for (size_t i = 0; i < rows[iRow].q2.size(); ++i) {
                        if (rows[iRow].theory[i] > 0.0) {
                            xTheory.push_back(rows[iRow].q2[i]);
                            yTheory.push_back(rows[iRow].theory[i] * offset);
                            yTheoryErr.push_back(rows[iRow].theoryErr[i] * offset);
                            const double relErr = (rows[iRow].theory[i] > 0.0)
                                                  ? (rows[iRow].theoryErr[i] / rows[iRow].theory[i])
                                                  : std::numeric_limits<double>::infinity();
                            if (std::isfinite(relErr) && relErr <= kMaxBandRelErr) {
                                xTheoryBand.push_back(rows[iRow].q2[i]);
                                yTheoryBand.push_back(rows[iRow].theory[i] * offset);
                                yTheoryBandErr.push_back(rows[iRow].theoryErr[i] * offset);
                            }
                        }
                        if (rows[iRow].data[i] > 0.0) {
                            xData.push_back(rows[iRow].q2[i]);
                            yData.push_back(rows[iRow].data[i] * offset);
                            yDataErr.push_back(rows[iRow].dataErr[i] * offset);
                        }
                    }
                    if (xTheory.empty() && xData.empty()) continue;

                    if (!xTheoryBand.empty()) {
                        auto* band = new TGraphErrors(static_cast<int>(xTheoryBand.size()));
                        for (size_t i = 0; i < xTheoryBand.size(); ++i) {
                            band->SetPoint(static_cast<int>(i), xTheoryBand[i], yTheoryBand[i]);
                            band->SetPointError(static_cast<int>(i), 0.0, yTheoryBandErr[i]);
                        }
                        band->SetFillColorAlpha(kGreen + 1, 0.35);
                        band->SetLineColor(kGreen + 1);
                        band->SetLineWidth(1);
                        band->Draw("3 SAME");
                        transientObjects.push_back(band);
                    }
                    if (!xTheory.empty()) {
                        auto* line = new TGraph(static_cast<int>(xTheory.size()));
                        for (size_t i = 0; i < xTheory.size(); ++i) {
                            line->SetPoint(static_cast<int>(i), xTheory[i], yTheory[i]);
                        }
                        line->SetLineColor(kGreen + 2);
                        line->SetLineWidth(2);
                        line->SetMarkerStyle(20);
                        line->SetMarkerSize(0.8);
                        line->SetMarkerColor(kGreen + 2);
                        line->Draw("LP SAME");

                        transientObjects.push_back(line);
                        if (!legendTheory) legendTheory = line;
                    }

                    if (!xData.empty()) {
                        auto* data = new TGraphErrors(static_cast<int>(xData.size()));
                        for (size_t i = 0; i < xData.size(); ++i) {
                            data->SetPoint(static_cast<int>(i), xData[i], yData[i]);
                            data->SetPointError(static_cast<int>(i), 0.0, yDataErr[i]);
                        }
                        data->SetLineColor(kBlack);
                        data->SetMarkerColor(kBlack);
                        data->SetMarkerStyle(20);
                        data->SetMarkerSize(0.9);
                        data->SetLineWidth(1);
                        data->Draw("P SAME");
                        transientObjects.push_back(data);
                        if (!legendData) legendData = data;
                    }

                    double yText = -1.0;
                    if (!yTheory.empty()) {
                        yText = yTheory.front();
                    } else if (!yData.empty()) {
                        yText = yData.front();
                    }
                    if (yText > 0.0) {
                        const double xMin = frame->GetXaxis()->GetXmin();
                        const double xMax = frame->GetXaxis()->GetXmax();
                        const double xText = std::exp(std::log(xMin) + 0.06 * (std::log(xMax) - std::log(xMin)));
                        TLatex betaText;
                        betaText.SetTextSize(0.028);
                        betaText.SetTextFont(42);
                        betaText.SetTextColor(kBlack);
                        betaText.SetTextAlign(12);
                        betaText.DrawLatex(
                            xText,
                            yText,
                            Form("%.4g < #beta < %.4g (i = %zu)",
                                 rows[iRow].betaLow,
                                 rows[iRow].betaHigh,
                                 iRow)
                        );
                    }
                }

                const LegendBox reducedLegendBox = getLegendBox(kReducedLegendOverrides,
                                                                kDefaultReducedLegendBox,
                                                                xpomTag);
                TLegend* leg = new TLegend(reducedLegendBox.x1, reducedLegendBox.y1,
                                           reducedLegendBox.x2, reducedLegendBox.y2);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.032);
                leg->SetHeader(Form("%.4g < x_{pom} < %.4g", xpomLow, xpomHigh), "C");
                if (legendTheory) leg->AddEntry(legendTheory, "MC", "lp");
                if (legendData) leg->AddEntry(legendData, dataLegendLabel, "p");
                leg->Draw();

                {
                    const std::string debugName = "figs/cross_sections/debug/sigma_r_debug_xpom_" + xpomTag + ".tsv";
                    std::ofstream dbg(debugName);
                    if (dbg.is_open()) {
                        dbg << "# iQ2\tiBeta\tQ2\tbeta\txpom\ty\td3_theory\td3_theory_err\td3_data\td3_data_err\tred_theory\tred_theory_err\tred_data\tred_data_err\n";
                        dbg << std::setprecision(8);
                        for (const auto& dp : debugPoints) {
                            dbg << dp.iQ2 << '\t'
                                << dp.iBeta << '\t'
                                << dp.q2 << '\t'
                                << dp.beta << '\t'
                                << xpomCenter << '\t'
                                << dp.y << '\t'
                                << dp.d3Theory << '\t'
                                << dp.d3TheoryErr << '\t'
                                << dp.d3Data << '\t'
                                << dp.d3DataErr << '\t'
                                << dp.redTheory << '\t'
                                << dp.redTheoryErr << '\t'
                                << dp.redData << '\t'
                                << dp.redDataErr << '\n';
                        }
                    }
                }

                DrawSimLabels(inputFile);
                SaveCanvas(cStack, saveName.c_str());

                delete leg;
                for (TObject* obj : transientObjects) {
                    delete obj;
                }
                delete frame;
                delete cStack;

                // Also write a true (non-reduced) stacked x_{pom} d^{3}#sigma plot
                // to the legacy filename for direct visual comparison.
                std::vector<SeriesRow> rowsRaw;
                rowsRaw.reserve(nBeta);
                for (int iBeta = 1; iBeta <= nBeta; ++iBeta) {
                    SeriesRow rowRaw;
                    rowRaw.betaLow = hD3MC->GetYaxis()->GetBinLowEdge(iBeta);
                    rowRaw.betaHigh = hD3MC->GetYaxis()->GetBinUpEdge(iBeta);
                    rowRaw.betaCenter = hD3MC->GetYaxis()->GetBinCenter(iBeta);
                    if (!std::isfinite(rowRaw.betaCenter) || rowRaw.betaCenter <= 0.0) continue;
                    for (const auto& dp : debugPoints) {
                        if (dp.iBeta != iBeta) continue;
                        const bool hasTheory = std::isfinite(dp.d3Theory) && dp.d3Theory > 0.0;
                        const bool hasData = std::isfinite(dp.d3Data) && dp.d3Data > 0.0;
                        if (!hasTheory && !hasData) continue;
                        rowRaw.q2.push_back(dp.q2);
                        rowRaw.theory.push_back(hasTheory ? (xpomCenter * dp.d3Theory) : -1.0);
                        rowRaw.theoryErr.push_back(hasTheory ? (xpomCenter * std::max(0.0, dp.d3TheoryErr)) : 0.0);
                        rowRaw.data.push_back(hasData ? (xpomCenter * dp.d3Data) : -1.0);
                        rowRaw.dataErr.push_back(hasData ? (xpomCenter * std::max(0.0, dp.d3DataErr)) : 0.0);
                    }
                    if (rowRaw.q2.size() >= 2) rowsRaw.push_back(std::move(rowRaw));
                }
                if (!rowsRaw.empty()) {
                    std::sort(rowsRaw.begin(), rowsRaw.end(), [](const SeriesRow& a, const SeriesRow& b) {
                        return a.betaCenter > b.betaCenter;
                    });
                    double yMinRaw = std::numeric_limits<double>::max();
                    double yMaxRaw = 0.0;
                    for (size_t iRow = 0; iRow < rowsRaw.size(); ++iRow) {
                        const double offset = std::pow(kOffsetBase, static_cast<double>(iRow));
                        for (size_t i = 0; i < rowsRaw[iRow].q2.size(); ++i) {
                            if (rowsRaw[iRow].theory[i] > 0.0) {
                                const double yVal = rowsRaw[iRow].theory[i] * offset;
                                yMinRaw = std::min(yMinRaw, yVal);
                                yMaxRaw = std::max(yMaxRaw, yVal);
                            }
                            if (rowsRaw[iRow].data[i] > 0.0) {
                                const double yVal = rowsRaw[iRow].data[i] * offset;
                                yMinRaw = std::min(yMinRaw, yVal);
                                yMaxRaw = std::max(yMaxRaw, yVal);
                            }
                        }
                    }
                    if (yMaxRaw > 0.0 && yMinRaw > 0.0 && std::isfinite(yMinRaw) && std::isfinite(yMaxRaw)) {
                        TCanvas* cStackRaw = new TCanvas(Form("c_d3sigma_q2_stacked_%d", iXpom),
                                                         "raw d3sigma stacked layout", 1100, 1400);
                        cStackRaw->SetLeftMargin(0.16);
                        cStackRaw->SetBottomMargin(0.12);
                        cStackRaw->SetLogx();
                        cStackRaw->SetLogy();
                        cStackRaw->SetGridx();
                        cStackRaw->SetGridy();
                        TH1D* frameRaw = new TH1D(Form("h_d3sigma_stack_frame_%d", iXpom),
                                                  "",
                                                  nQ2, q2Edges.data());
                        frameRaw->SetDirectory(nullptr);
                        frameRaw->SetStats(false);
                        frameRaw->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
                        frameRaw->GetYaxis()->SetTitle("3^{i} #times x_{pom} d^{3}#sigma/(dQ^{2} d#beta dx_{pom})");
                        frameRaw->GetXaxis()->SetTitleSize(0.048);
                        frameRaw->GetYaxis()->SetTitleSize(0.048);
                        frameRaw->GetXaxis()->SetLabelSize(0.040);
                        frameRaw->GetYaxis()->SetLabelSize(0.040);
                        frameRaw->GetXaxis()->SetTitleOffset(1.05);
                        frameRaw->GetYaxis()->SetTitleOffset(1.30);
                        frameRaw->SetMinimum(std::max(yMinRaw / 2.0, 1e-12));
                        frameRaw->SetMaximum(yMaxRaw * 1.35);
                        frameRaw->Draw("AXIS");

                        std::vector<TObject*> rawObjects;
                        rawObjects.reserve(rowsRaw.size() * 3);
                        TGraph* legendTheoryRaw = nullptr;
                        TGraphErrors* legendDataRaw = nullptr;

                        for (size_t iRow = 0; iRow < rowsRaw.size(); ++iRow) {
                            const double offset = std::pow(kOffsetBase, static_cast<double>(iRow));
                            std::vector<double> xTheory, yTheory, yTheoryErr;
                            std::vector<double> xTheoryBand, yTheoryBand, yTheoryBandErr;
                            std::vector<double> xData, yData, yDataErr;
                            xTheory.reserve(rowsRaw[iRow].q2.size());
                            yTheory.reserve(rowsRaw[iRow].q2.size());
                            yTheoryErr.reserve(rowsRaw[iRow].q2.size());
                            xTheoryBand.reserve(rowsRaw[iRow].q2.size());
                            yTheoryBand.reserve(rowsRaw[iRow].q2.size());
                            yTheoryBandErr.reserve(rowsRaw[iRow].q2.size());
                            xData.reserve(rowsRaw[iRow].q2.size());
                            yData.reserve(rowsRaw[iRow].q2.size());
                            yDataErr.reserve(rowsRaw[iRow].q2.size());
                            for (size_t i = 0; i < rowsRaw[iRow].q2.size(); ++i) {
                                if (rowsRaw[iRow].theory[i] > 0.0) {
                                    xTheory.push_back(rowsRaw[iRow].q2[i]);
                                    yTheory.push_back(rowsRaw[iRow].theory[i] * offset);
                                    yTheoryErr.push_back(rowsRaw[iRow].theoryErr[i] * offset);
                                    const double relErr = (rowsRaw[iRow].theory[i] > 0.0)
                                                          ? (rowsRaw[iRow].theoryErr[i] / rowsRaw[iRow].theory[i])
                                                          : std::numeric_limits<double>::infinity();
                                    if (std::isfinite(relErr) && relErr <= kMaxBandRelErr) {
                                        xTheoryBand.push_back(rowsRaw[iRow].q2[i]);
                                        yTheoryBand.push_back(rowsRaw[iRow].theory[i] * offset);
                                        yTheoryBandErr.push_back(rowsRaw[iRow].theoryErr[i] * offset);
                                    }
                                }
                                if (rowsRaw[iRow].data[i] > 0.0) {
                                    xData.push_back(rowsRaw[iRow].q2[i]);
                                    yData.push_back(rowsRaw[iRow].data[i] * offset);
                                    yDataErr.push_back(rowsRaw[iRow].dataErr[i] * offset);
                                }
                            }
                            if (xTheory.empty() && xData.empty()) continue;

                            if (!xTheoryBand.empty()) {
                                auto* band = new TGraphErrors(static_cast<int>(xTheoryBand.size()));
                                for (size_t i = 0; i < xTheoryBand.size(); ++i) {
                                    band->SetPoint(static_cast<int>(i), xTheoryBand[i], yTheoryBand[i]);
                                    band->SetPointError(static_cast<int>(i), 0.0, yTheoryBandErr[i]);
                                }
                                band->SetFillColorAlpha(kGreen + 1, 0.35);
                                band->SetLineColor(kGreen + 1);
                                band->SetLineWidth(1);
                                band->Draw("3 SAME");
                                rawObjects.push_back(band);
                            }
                            if (!xTheory.empty()) {
                                auto* line = new TGraph(static_cast<int>(xTheory.size()));
                                for (size_t i = 0; i < xTheory.size(); ++i) {
                                    line->SetPoint(static_cast<int>(i), xTheory[i], yTheory[i]);
                                }
                                line->SetLineColor(kGreen + 2);
                                line->SetLineWidth(2);
                                line->SetMarkerStyle(20);
                                line->SetMarkerSize(0.8);
                                line->SetMarkerColor(kGreen + 2);
                                line->Draw("LP SAME");
                                rawObjects.push_back(line);
                                if (!legendTheoryRaw) legendTheoryRaw = line;
                            }
                            if (!xData.empty()) {
                                auto* data = new TGraphErrors(static_cast<int>(xData.size()));
                                for (size_t i = 0; i < xData.size(); ++i) {
                                    data->SetPoint(static_cast<int>(i), xData[i], yData[i]);
                                    data->SetPointError(static_cast<int>(i), 0.0, yDataErr[i]);
                                }
                                data->SetLineColor(kBlack);
                                data->SetMarkerColor(kBlack);
                                data->SetMarkerStyle(20);
                                data->SetMarkerSize(0.9);
                                data->SetLineWidth(1);
                                data->Draw("P SAME");
                                rawObjects.push_back(data);
                                if (!legendDataRaw) legendDataRaw = data;
                            }

                            double yText = -1.0;
                            if (!yTheory.empty()) {
                                yText = yTheory.front();
                            } else if (!yData.empty()) {
                                yText = yData.front();
                            }
                            if (yText > 0.0) {
                                const double xMin = frameRaw->GetXaxis()->GetXmin();
                                const double xMax = frameRaw->GetXaxis()->GetXmax();
                                const double xText = std::exp(std::log(xMin) + 0.06 * (std::log(xMax) - std::log(xMin)));
                                TLatex betaText;
                                betaText.SetTextSize(0.028);
                                betaText.SetTextFont(42);
                                betaText.SetTextColor(kBlack);
                                betaText.SetTextAlign(12);
                                betaText.DrawLatex(
                                    xText,
                                    yText,
                                    Form("%.4g < #beta < %.4g (i = %zu)",
                                         rowsRaw[iRow].betaLow,
                                         rowsRaw[iRow].betaHigh,
                                         iRow)
                                );
                            }
                        }

                        const LegendBox rawLegendBox = getLegendBox(kRawLegendOverrides,
                                                                    kDefaultRawLegendBox,
                                                                    xpomTag);
                        TLegend* legRaw = new TLegend(rawLegendBox.x1, rawLegendBox.y1,
                                                      rawLegendBox.x2, rawLegendBox.y2);
                        legRaw->SetBorderSize(0);
                        legRaw->SetFillStyle(0);
                        legRaw->SetTextFont(42);
                        legRaw->SetTextSize(0.032);
                        legRaw->SetHeader(Form("%.4g < x_{pom} < %.4g", xpomLow, xpomHigh), "C");
                        if (legendTheoryRaw) legRaw->AddEntry(legendTheoryRaw, "MC", "lp");
                        if (legendDataRaw) legRaw->AddEntry(legendDataRaw, dataLegendLabel, "p");
                        legRaw->Draw();

                        DrawSimLabels(inputFile);
                        SaveCanvas(cStackRaw, legacySaveName.c_str());

                        delete legRaw;
                        for (TObject* obj : rawObjects) delete obj;
                        delete frameRaw;
                        delete cStackRaw;
                    }
                }
            }
        }
    }
    delete hD3SumLocal;


    inputFile->Close();
    delete inputFile;

    return 0;
}
