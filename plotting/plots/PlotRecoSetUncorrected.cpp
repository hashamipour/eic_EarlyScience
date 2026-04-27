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

void PlotRecoSetUncorrected(TFile* inputFile,
                                   const char* recoHistNameSetA,
                                   const char* recoHistNameSetB,
                                   const char* xLabel,
                                   const char* title,
                                   const char* saveName,
                                   bool logX,
                                   bool logY,
                                   const char* setALabel,
                                   const char* setBLabel,
                                   bool normalizeToPDF) {
    if (!inputFile) return;
    TH1* h_setA = (TH1*)inputFile->Get(recoHistNameSetA);
    TH1* h_setB = (TH1*)inputFile->Get(recoHistNameSetB);
    if (!h_setA || !h_setB) {
        Logger::warning(std::string("Missing hist(s) for uncorrected reco comparison: ") +
                        recoHistNameSetA + ", " + recoHistNameSetB);
        return;
    }

    TH1D* hA = (TH1D*)h_setA->Clone(Form("uncorr_setA_%s", recoHistNameSetA));
    TH1D* hB = (TH1D*)h_setB->Clone(Form("uncorr_setB_%s", recoHistNameSetB));
    hA->SetDirectory(nullptr);
    hB->SetDirectory(nullptr);

    if (normalizeToPDF) {
        const double intA = hA->Integral();
        const double intB = hB->Integral();
        if (intA > 0.0) hA->Scale(1.0 / intA);
        if (intB > 0.0) hB->Scale(1.0 / intB);
    }

    TCanvas* c = new TCanvas(Form("c_uncorr_%s", recoHistNameSetA), title, 1200, 900);
    gStyle->SetOptStat(0);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    hA->SetStats(false);
    hA->SetTitle("");
    hA->GetXaxis()->SetTitle(xLabel);
    hA->GetYaxis()->SetTitle(normalizeToPDF ? "Normalized" : "Number of events");
    hA->GetXaxis()->SetTitleOffset(1.1);
    hA->GetYaxis()->SetTitleOffset(1.2);

    hA->SetLineWidth(2);
    hA->SetLineColor(kBlack);
    hA->SetFillColor(kYellow - 9);
    hA->SetFillStyle(1001);
    hB->SetMarkerStyle(20);
    hB->SetMarkerSize(1.4);
    hB->SetMarkerColor(kRed + 1);
    hB->SetLineColor(kRed + 1);

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= hA->GetNbinsX(); ++i) {
        double v = std::max(hA->GetBinContent(i), hB->GetBinContent(i));
        if (v > max_val) max_val = v;
        if (v > 0 && v < min_pos) min_pos = v;
    }
    if (logY && min_pos < std::numeric_limits<double>::max()) {
        hA->SetMinimum(min_pos * 0.5);
        hA->SetMaximum(max_val * 5.0);
    } else {
        hA->SetMaximum(max_val * 1.3);
    }

    hA->Draw("hist");
    hB->Draw("pe same");

    TLegend* leg = new TLegend(0.65, 0.75, 0.92, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hA, setALabel, "l");
    leg->AddEntry(hB, setBLabel, "p");
    leg->Draw();

    SaveCanvas(c, saveName);
    delete c;
}
