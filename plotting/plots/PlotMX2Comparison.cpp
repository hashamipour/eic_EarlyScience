// Hand-styled MX^2 comparison: overlays hadronic-sum and kinematic
// reconstructions on truth and reco with explicit, distinct colors. The
// PlotOptions1D auto-styler collapses all "*truth*" curves to the same color,
// so we draw this one by hand.

#include "plots/DDISPlots.hpp"

#include "Plotting.hpp"
#include "Utility.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>

void PlotMX2Comparison(TFile* inputFile,
                       const std::string& outpath,
                       bool logY) {
    if (!inputFile) return;

    auto get1d = [&](const char* name) -> TH1D* {
        return dynamic_cast<TH1D*>(inputFile->Get(name));
    };

    TH1D* h_truth_had = get1d("MX2_truth");
    TH1D* h_reco_had  = get1d("MX2_reco");
    TH1D* h_truth_kin = get1d("MX2_truth_kin");
    TH1D* h_reco_kin  = get1d("MX2_reco_kin");

    if (!h_truth_had || !h_reco_had || !h_truth_kin || !h_reco_kin) {
        Logger::warning("PlotMX2Comparison: missing one of MX2_truth, MX2_reco, "
                        "MX2_truth_kin, MX2_reco_kin -- skipping " + outpath);
        return;
    }

    TCanvas* c = new TCanvas(Form("c_mx2_cmp_%s", logY ? "logy" : "liny"),
                             "MX2 comparison", 1200, 900);
    gStyle->SetOptTitle(0);
    c->SetLogx();
    if (logY) c->SetLogy();
    c->SetGridx();
    c->SetGridy();
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.12);

    auto style_truth = [](TH1D* h, Color_t col) {
        h->SetTitle("");
        h->SetLineColor(col);
        h->SetLineWidth(3);
        h->SetMarkerStyle(0);
        h->SetMarkerSize(0.0);
        h->SetFillStyle(0);
        h->SetStats(false);
    };
    auto style_reco = [](TH1D* h, Color_t col, int marker) {
        h->SetTitle("");
        h->SetLineColor(col);
        h->SetMarkerColor(col);
        h->SetLineWidth(2);
        h->SetMarkerStyle(marker);
        h->SetMarkerSize(0.9);
        h->SetFillStyle(0);
        h->SetStats(false);
    };

    style_truth(h_truth_had, kBlack);
    style_truth(h_truth_kin, kGreen + 2);
    style_reco (h_reco_had,  kBlue + 1, 20);
    style_reco (h_reco_kin,  kRed + 1,  21);

    double ymax = 0.0;
    for (TH1D* h : {h_truth_had, h_reco_had, h_truth_kin, h_reco_kin}) {
        ymax = std::max(ymax, h->GetMaximum());
    }
    if (ymax <= 0.0) ymax = 1.0;

    h_truth_had->GetXaxis()->SetTitle("M_{X}^{2} [GeV^{2}]");
    h_truth_had->GetYaxis()->SetTitle("Number of events");
    h_truth_had->GetXaxis()->SetTitleSize(0.045);
    h_truth_had->GetYaxis()->SetTitleSize(0.045);
    h_truth_had->GetXaxis()->SetLabelSize(0.040);
    h_truth_had->GetYaxis()->SetLabelSize(0.040);
    h_truth_had->GetXaxis()->SetRangeUser(1.0e-3, 1000.0);
    if (logY) {
        h_truth_had->SetMinimum(0.5);
        h_truth_had->SetMaximum(ymax * 5.0);
    } else {
        h_truth_had->SetMinimum(0.0);
        h_truth_had->SetMaximum(ymax * 1.25);
    }

    h_truth_had->Draw("HIST");
    h_truth_kin->Draw("HIST SAME");
    h_reco_had ->Draw("PE SAME");
    h_reco_kin ->Draw("PE SAME");

    TLegend* leg = new TLegend(logY ? 0.16 : 0.55, 0.66,
                               logY ? 0.46 : 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(h_truth_had, "MC Truth (hadronic sum)", "l");
    leg->AddEntry(h_truth_kin, "MC Truth (kinematic)",    "l");
    leg->AddEntry(h_reco_had,  "Reco (hadronic sum)",     "lep");
    leg->AddEntry(h_reco_kin,  "Reco (kinematic, RP #cup B0)", "lep");
    leg->Draw();

    DrawSimLabels(inputFile);
    SaveCanvas(c, outpath.c_str());

    delete leg;
    delete c;
}
