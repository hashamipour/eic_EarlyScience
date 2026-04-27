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

void PlotRecoSetComparison(TFile* inputFile,
                                  const char* effTruthHistNameSetA,
                                  const char* effRecoHistNameSetA,
                                  const char* recoHistNameSetA,
                                  const char* recoHistNameSetB,
                                  const char* xLabel,
                                  const char* title,
                                  const char* saveName,
                                  const bool logX,
                                  const bool logY,
                                  const char* setALabel,
                                  const char* setBLabel) {
    if (!inputFile) return;
    TH1* h_truth_setA = (TH1*)inputFile->Get(effTruthHistNameSetA);
    TH1* h_reco_setA_for_eff = (TH1*)inputFile->Get(effRecoHistNameSetA);
    TH1* h_setA = (TH1*)inputFile->Get(recoHistNameSetA);
    TH1* h_setB = (TH1*)inputFile->Get(recoHistNameSetB);
    if (!h_truth_setA || !h_reco_setA_for_eff || !h_setA || !h_setB) {
        Logger::warning(std::string("Missing hist(s) for common-eff reco comparison: ") +
                        effTruthHistNameSetA + ", " + effRecoHistNameSetA + ", " +
                        recoHistNameSetA + ", " + recoHistNameSetB);
        return;
    }

    TH1D* h_eff = (TH1D*)h_reco_setA_for_eff->Clone(Form("eff_ref_%s", effRecoHistNameSetA));
    h_eff->SetDirectory(nullptr);
    h_eff->Divide(h_reco_setA_for_eff, h_truth_setA, 1.0, 1.0, "B");

    TH1D* hA_corr = (TH1D*)h_setA->Clone(Form("setA_corr_%s", recoHistNameSetA));
    TH1D* hB_corr = (TH1D*)h_setB->Clone(Form("setB_corr_%s", recoHistNameSetB));
    hA_corr->SetDirectory(nullptr);
    hB_corr->SetDirectory(nullptr);

    const int nbins = hA_corr->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        const double eff = h_eff->GetBinContent(i);
        const double a = hA_corr->GetBinContent(i);
        const double b = hB_corr->GetBinContent(i);
        const double a_err = hA_corr->GetBinError(i);
        const double b_err = hB_corr->GetBinError(i);
        if (eff > 0.0) {
            hA_corr->SetBinContent(i, a / eff);
            hA_corr->SetBinError(i, a_err / eff);
            hB_corr->SetBinContent(i, b / eff);
            hB_corr->SetBinError(i, b_err / eff);
        } else {
            hA_corr->SetBinContent(i, 0.0);
            hA_corr->SetBinError(i, 0.0);
            hB_corr->SetBinContent(i, 0.0);
            hB_corr->SetBinError(i, 0.0);
        }
    }

    TCanvas* c = new TCanvas(Form("c_cmp_%s", recoHistNameSetA), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    hA_corr->SetStats(false);
    hA_corr->SetTitle("");
    hA_corr->GetXaxis()->SetTitle(xLabel);
    hA_corr->GetYaxis()->SetTitle("Efficiency-corrected counts");
    hA_corr->GetXaxis()->SetTitleOffset(1.1);
    hA_corr->GetYaxis()->SetTitleOffset(1.2);

    hA_corr->SetLineWidth(2);
    hA_corr->SetLineColor(kBlack);
    hA_corr->SetFillColor(kYellow - 9);
    hA_corr->SetFillStyle(1001);
    hB_corr->SetMarkerStyle(20);
    hB_corr->SetMarkerSize(1.1);
    hB_corr->SetMarkerColor(kBlack);
    hB_corr->SetLineColor(kBlack);

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= hA_corr->GetNbinsX(); ++i) {
        const double a = hA_corr->GetBinContent(i);
        const double b = hB_corr->GetBinContent(i);
        max_val = std::max(max_val, a + hA_corr->GetBinError(i));
        max_val = std::max(max_val, b + hB_corr->GetBinError(i));
        if (a > 0.0 && a < min_pos) min_pos = a;
        if (b > 0.0 && b < min_pos) min_pos = b;
    }
    if (max_val > 0.0) {
        hA_corr->SetMaximum(max_val * 1.25);
        if (logY) {
            const double min_y = (min_pos < std::numeric_limits<double>::max()) ? (min_pos * 0.5) : 1e-6;
            hA_corr->SetMinimum(std::max(min_y, 1e-6));
        } else {
            hA_corr->SetMinimum(0.0);
        }
    }

    hA_corr->Draw("HIST");
    hB_corr->Draw("PE SAME");

    TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hA_corr, setALabel, "l");
    legend->AddEntry(hB_corr, setBLabel, "p");
    legend->AddEntry((TObject*)nullptr, "Common efficiency from Set A", "");
    legend->Draw();

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
    delete h_eff;
    delete hA_corr;
    delete hB_corr;
}
