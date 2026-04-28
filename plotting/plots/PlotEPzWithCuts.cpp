#include "plots/DDISPlots.hpp"

#include "Plotting.hpp"
#include "Utility.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TParameter.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <limits>

void PlotEPzWithCuts(TFile* inputFile) {
    if (!inputFile) return;
    // Compare Set A (MC) vs Set B (pseudo-data) using histograms filled BEFORE the
    // E-p_z cut, so the accepted Sigma(E-p_z) window can be drawn on top to show
    // how many events would be removed by the cut.
    //
    // Sigma(E-pz) is computed by ElectronID::ComputeEventDeltaH() over
    // ReconstructedParticles (the EID stream), so the cut shown here matches
    // the paper A.1 event-level requirement. Falls back to the old
    // TTreeReader-matched EPz_reco_{mc,pdata} histograms if the EID
    // versions aren't present (e.g. skim ran with --no-eid).
    TH1D* h_setA = (TH1D*)inputFile->Get("EPz_eid_mc");
    TH1D* h_setB = (TH1D*)inputFile->Get("EPz_eid_pdata");
    if (!h_setA || !h_setB) {
        Logger::warning("EPz_eid_{mc,pdata} not found; falling back to EPz_reco_{mc,pdata}.");
        h_setA = (TH1D*)inputFile->Get("EPz_reco_mc");
        h_setB = (TH1D*)inputFile->Get("EPz_reco_pdata");
    }
    if (!h_setA || !h_setB) {
        Logger::warning("EPz histograms not found; skipping EPz cut plot.");
        return;
    }

    double epz_lo = 16.0, epz_hi = 24.0;
    if (auto* p = dynamic_cast<TParameter<double>*>(inputFile->Get("DISCut_epz_min")))
        epz_lo = p->GetVal();
    if (auto* p = dynamic_cast<TParameter<double>*>(inputFile->Get("DISCut_epz_max")))
        epz_hi = p->GetVal();

    TH1D* hA = (TH1D*)h_setA->Clone("epzcut_setA");
    TH1D* hB = (TH1D*)h_setB->Clone("epzcut_setB");
    hA->SetDirectory(nullptr);
    hB->SetDirectory(nullptr);

    TCanvas* c = new TCanvas("c_epz_cuts", "", 1200, 900);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);
    c->SetLogy();

    hA->SetStats(false);
    hA->SetTitle("");
    hA->GetXaxis()->SetTitle("#Sigma(E-p_{z}) [GeV]");
    hA->GetYaxis()->SetTitle("Number of events");
    hA->GetXaxis()->SetTitleOffset(1.1);
    hA->GetYaxis()->SetTitleOffset(1.2);

    hA->SetLineWidth(2);
    hA->SetLineColor(kBlack);
    hA->SetFillColor(kYellow - 9);
    hA->SetFillStyle(1001);
    hB->SetMarkerStyle(20);
    hB->SetMarkerSize(1.1);
    hB->SetMarkerColor(kRed + 1);
    hB->SetLineColor(kBlack);

    double max_val = 0.0;
    double min_pos = std::numeric_limits<double>::max();
    for (int i = 1; i <= hA->GetNbinsX(); ++i) {
        const double v = std::max(hA->GetBinContent(i), hB->GetBinContent(i));
        if (v > max_val) max_val = v;
        if (v > 0 && v < min_pos) min_pos = v;
    }
    const double ymin_draw = (min_pos < std::numeric_limits<double>::max()) ? std::max(min_pos * 0.5, 0.5) : 0.5;
    const double ymax_draw = max_val * 5.0;
    hA->SetMinimum(ymin_draw);
    hA->SetMaximum(ymax_draw);

    hA->Draw("hist");
    hB->Draw("pe same");

    TBox* box = new TBox(epz_lo, ymin_draw, epz_hi, ymax_draw);
    box->SetFillColorAlpha(kGreen + 1, 0.15);
    box->SetLineWidth(0);
    box->Draw("same");

    TLine* line_lo = new TLine(epz_lo, ymin_draw, epz_lo, ymax_draw);
    line_lo->SetLineColor(kGreen + 2);
    line_lo->SetLineStyle(2);
    line_lo->SetLineWidth(2);
    line_lo->Draw("same");

    TLine* line_hi = new TLine(epz_hi, ymin_draw, epz_hi, ymax_draw);
    line_hi->SetLineColor(kGreen + 2);
    line_hi->SetLineStyle(2);
    line_hi->SetLineWidth(2);
    line_hi->Draw("same");

    hA->Draw("hist same");
    hB->Draw("pe same");

    TLegend* leg = new TLegend(0.58, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hA, "MC", "l");
    leg->AddEntry(hB, "pseudo-data", "pe");
    leg->AddEntry(line_lo, Form("Cuts: %.0f < #Sigma(E-p_{z}) < %.0f GeV", epz_lo, epz_hi), "l");
    leg->Draw();

    DrawSimLabels(inputFile);
    c->Update();
    SaveCanvas(c, "figs/inclusive/distributions/EPz_distribution_logY.png");

    delete leg; delete box; delete line_lo; delete line_hi; delete c;
}
