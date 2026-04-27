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

void PlotRecoSetComparisonBR(TFile* inputFile,
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
        Logger::warning(std::string("Missing hist(s) for common-eff B0/RP reco comparison: ") +
                        effTruthHistNameSetA_B0 + ", " + effRecoHistNameSetA_B0 + ", " +
                        recoHistNameSetA_B0 + ", " + recoHistNameSetB_B0 + ", " +
                        effTruthHistNameSetA_RP + ", " + effRecoHistNameSetA_RP + ", " +
                        recoHistNameSetA_RP + ", " + recoHistNameSetB_RP);
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
