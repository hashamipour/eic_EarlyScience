#include "plots/DDISPlots.hpp"

#include "Plotting.hpp"
#include "Utility.hpp"
#include "YAMLBinning.hpp"
#include "PlotDrawing.hpp"
#include "GridDrawing.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

void PlotPhaseSpaceSlices(TFile* inputFile, const std::string& yamlPath) {
    if (!inputFile) return;
    TH3D* h3 = (TH3D*)inputFile->Get("phase3D_reco");
    if (!h3) {
        Logger::warning("phase3D_reco not found; skipping phase-space slices.");
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
            Logger::info("Using YAML bins for phase slices: " + yamlPath);
        } else {
            Logger::warning("YAML bins invalid; using default slice edges.");
            yaml_bins.clear();
        }
    }

    DrawSliceGrid(h3, 3, beta_edges,
                  "xy",
                  "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                  "figs/diffractive/distributions/phase_slices_beta.png",
                  true, true,
                  (yaml_requested ? nullptr : &xpom_overlay_bins),
                  (yaml_requested ? nullptr : &q2_overlay_bins),
                  true, true,
                  yaml_requested ? &yaml_bins : nullptr,
                  false);

    DrawSliceGrid(h3, 3, beta_edges,
                  "xy",
                  "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                  "figs/diffractive/distributions/phase_slices_beta_bin_number.png",
                  true, true,
                  (yaml_requested ? nullptr : &xpom_overlay_bins),
                  (yaml_requested ? nullptr : &q2_overlay_bins),
                  true, true,
                  yaml_requested ? &yaml_bins : nullptr,
                  true);

    // Bin-wise acceptance/purity/efficiency overlays using explicit YAML bins.
    // Definitions (user-requested):
    //   A_i   = N_meas_i / N_gen_i
    //   P_i   = N_gen&meas_same_i / N_meas_i
    //   eff_i = N_gen&meas_same_i / N_gen_i
    if (yaml_requested && !yaml_bins.empty()) {
        TH1* h_gen = (TH1*)inputFile->Get("phase_bin_gen");
        TH1* h_meas = (TH1*)inputFile->Get("phase_bin_meas");
        TH1* h_same = (TH1*)inputFile->Get("phase_bin_gen_meas_same");
        if (!h_gen || !h_meas || !h_same) {
            Logger::warning("phase-bin metric histograms not found; skipping "
                            "phase_slices_beta_{eff,acceptance,purity}.png");
        } else {
            const int nGlobalBins =
                static_cast<int>((q2_edges.size() - 1) * (xpom_edges.size() - 1) * (beta_edges.size() - 1));
            if (nGlobalBins > 0) {
                std::vector<double> acceptance_values(nGlobalBins, -1.0);
                std::vector<double> purity_values(nGlobalBins, -1.0);
                std::vector<double> efficiency_values(nGlobalBins, -1.0);

                const int nCommon = std::min({nGlobalBins, h_gen->GetNbinsX(), h_meas->GetNbinsX(), h_same->GetNbinsX()});
                for (int k = 0; k < nCommon; ++k) {
                    const double n_gen = h_gen->GetBinContent(k + 1);
                    const double n_meas = h_meas->GetBinContent(k + 1);
                    const double n_same = h_same->GetBinContent(k + 1);

                    if (n_gen > 0.0) {
                        acceptance_values[k] = n_meas / n_gen;
                        efficiency_values[k] = n_same / n_gen;
                    }
                    if (n_meas > 0.0) {
                        purity_values[k] = n_same / n_meas;
                    }
                }

                // Default: no bin IDs
                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_eff.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        efficiency_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 1.0,
                                        nullptr, 1.0, false);

                // Variant with bin numbers
                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_eff_binnum.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        efficiency_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 1.0);

                // Default: no bin IDs
                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_acceptance.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        acceptance_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 1.2,
                                        nullptr, 1.0, false);

                // Variant with bin numbers
                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_acceptance_binnum.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        acceptance_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 1.2);

                // Default: no bin IDs
                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_purity.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        purity_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 1.0,
                                        nullptr, 1.0, false);

                // Variant with bin numbers
                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_purity_binnum.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        purity_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 1.0);
            }
        }

        // Closure test with user-requested formula:
        // N_i^(B,corr) = N_i^(B,meas) * P_i^A / eff_i^A
        // where:
        //   P_i^A   = N_i^(A,gen&meas_same) / N_i^(A,meas)
        //   eff_i^A = N_i^(A,gen&meas_same) / N_i^(A,gen)
        // and closure ratio:
        //   C_i = N_i^(B,corr) / N_i^(B,gen)
        TH1* h_gen_A = (TH1*)inputFile->Get("phase_bin_gen_setA");
        TH1* h_meas_A = (TH1*)inputFile->Get("phase_bin_meas_setA");
        TH1* h_same_A = (TH1*)inputFile->Get("phase_bin_gen_meas_same_setA");
        TH1* h_gen_B = (TH1*)inputFile->Get("phase_bin_gen_setB");
        TH1* h_meas_B = (TH1*)inputFile->Get("phase_bin_meas_setB");
        if (h_gen_A && h_meas_A && h_same_A && h_gen_B && h_meas_B) {
            const int nGlobalBins =
                static_cast<int>((q2_edges.size() - 1) * (xpom_edges.size() - 1) * (beta_edges.size() - 1));
            if (nGlobalBins > 0) {
                std::vector<double> closure_values(nGlobalBins, -1.0);
                std::vector<double> closure_unc_values(nGlobalBins, -1.0);
                const int nCommon = std::min({nGlobalBins,
                                              h_gen_A->GetNbinsX(),
                                              h_meas_A->GetNbinsX(),
                                              h_same_A->GetNbinsX(),
                                              h_gen_B->GetNbinsX(),
                                              h_meas_B->GetNbinsX()});
                for (int k = 0; k < nCommon; ++k) {
                    const double n_gen_A = h_gen_A->GetBinContent(k + 1);
                    const double n_meas_A = h_meas_A->GetBinContent(k + 1);
                    const double n_same_A = h_same_A->GetBinContent(k + 1);
                    const double n_gen_B = h_gen_B->GetBinContent(k + 1);
                    const double n_meas_B = h_meas_B->GetBinContent(k + 1);

                    const double purity_A = (n_meas_A > 0.0) ? (n_same_A / n_meas_A) : -1.0;
                    const double eff_A = (n_gen_A > 0.0) ? (n_same_A / n_gen_A) : -1.0;
                    if (eff_A > 0.0 && purity_A >= 0.0 && n_gen_B > 0.0) {
                        const double n_corr_B = n_meas_B * purity_A / eff_A;
                        const double closure = n_corr_B / n_gen_B;
                        closure_values[k] = closure;

                        // Approximate statistical uncertainty (Poisson, uncorrelated counts)
                        // using algebraically equivalent form:
                        // C = (N_meas^B * N_gen^A) / (N_meas^A * N_gen^B)
                        if (n_meas_B > 0.0 && n_gen_A > 0.0 && n_meas_A > 0.0 && n_gen_B > 0.0) {
                            const double rel2 = (1.0 / n_meas_B) +
                                                (1.0 / n_gen_A) +
                                                (1.0 / n_meas_A) +
                                                (1.0 / n_gen_B);
                            closure_unc_values[k] = std::abs(closure) * std::sqrt(rel2);
                        }
                    }
                }

                DrawSliceGridWithMetric(h3, 3, beta_edges,
                                        "xy",
                                        "x_{pom}", "Q^{2} [GeV^{2}]", "#beta",
                                        "figs/diffractive/distributions/phase_slices_beta_closure_ratio.png",
                                        true, true,
                                        yaml_bins,
                                        true, true,
                                        closure_values,
                                        q2_edges, xpom_edges, beta_edges,
                                        0.0, 2.0,
                                        &closure_unc_values, 0.85);
            }
        } else {
            Logger::warning("Set A/B phase-bin histograms for closure test are missing; "
                            "skipping phase_slices_beta_closure_ratio.png");
        }
    }
}
