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

void PlotRecoSetComparisonBRWithSetBUncorrected(TFile* inputFile,
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
        Logger::warning(std::string("Missing hist(s) for B0/RP set-B corrected/uncorrected comparison: ") +
                        effTruthHistNameSetA_B0 + ", " + effRecoHistNameSetA_B0 + ", " +
                        recoHistNameSetA_B0 + ", " + recoHistNameSetB_B0 + ", " +
                        effTruthHistNameSetA_RP + ", " + effRecoHistNameSetA_RP + ", " +
                        recoHistNameSetA_RP + ", " + recoHistNameSetB_RP);
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
    hA_sum_reco->GetYaxis()->SetTitle("Number of events");
    hA_sum_reco->GetXaxis()->SetTitleOffset(1.1);
    hA_sum_reco->GetYaxis()->SetTitleOffset(1.2);
    hA_sum_reco->SetLineColor(kBlack);
    hA_sum_reco->SetLineWidth(2);
    hA_sum_reco->SetLineStyle(1);
    hA_sum_reco->SetFillColor(kYellow - 9);
    hA_sum_reco->SetFillStyle(1001);

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
