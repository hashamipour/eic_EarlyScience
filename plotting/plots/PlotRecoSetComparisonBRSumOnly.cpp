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

void PlotRecoSetComparisonBRSumOnly(TFile* inputFile,
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
        Logger::warning("Missing hist(s) for sum-only B0/RP reco comparison.");
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
    h_sum_A_corr->GetYaxis()->SetTitle("Number of events");
    h_sum_A_corr->GetXaxis()->SetTitleOffset(1.1);
    h_sum_A_corr->GetYaxis()->SetTitleOffset(1.2);

    h_sum_A_corr->SetLineColor(kBlack);
    h_sum_A_corr->SetLineStyle(1);
    h_sum_A_corr->SetLineWidth(2);
    h_sum_A_corr->SetFillColor(kYellow - 9);
    h_sum_A_corr->SetFillStyle(1001);

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
