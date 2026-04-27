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

void Plot3DResponseMatrix(TFile* inputFile) {
    if (!inputFile) return;

    // --- (a) Row-normalized response matrix ---
    TH2D* h_raw = (TH2D*)inputFile->Get("Response_3D");
    if (!h_raw) {
        Logger::warning("Response_3D not found; skipping 3D response matrix plots.");
        return;
    }

    const int nBins = h_raw->GetNbinsX();
    TH2D* h_rowNorm = (TH2D*)h_raw->Clone("Response_3D_rowNorm_plot");
    h_rowNorm->SetDirectory(nullptr);
    for (int i = 1; i <= nBins; ++i) {
        double rowSum = 0.0;
        for (int j = 1; j <= nBins; ++j) rowSum += h_raw->GetBinContent(i, j);
        if (rowSum > 0.0) {
            for (int j = 1; j <= nBins; ++j)
                h_rowNorm->SetBinContent(i, j, h_raw->GetBinContent(i, j) / rowSum);
        }
    }

    TCanvas* c1 = new TCanvas("c_resp3d", "", 1200, 1100);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.12);
    c1->SetTopMargin(0.08);
    c1->SetBottomMargin(0.12);
    gStyle->SetPaintTextFormat(".0f");

    h_rowNorm->SetTitle("");
    h_rowNorm->GetXaxis()->SetTitle("Truth bin index");
    h_rowNorm->GetYaxis()->SetTitle("Reco bin index");
    h_rowNorm->GetZaxis()->SetTitle("Migration probability");
    h_rowNorm->GetXaxis()->SetTitleOffset(1.2);
    h_rowNorm->GetYaxis()->SetTitleOffset(1.3);
    h_rowNorm->Draw("COLZ");

    TLine* diag = new TLine(0.5, 0.5, nBins + 0.5, nBins + 0.5);
    diag->SetLineColor(kRed);
    diag->SetLineStyle(2);
    diag->SetLineWidth(2);
    diag->Draw("same");

    DrawSimLabels(inputFile);
    c1->Update();
    SaveCanvas(c1, "figs/diffractive/response/response_3d_rowNorm.png");
    delete diag; delete c1;

    // --- (b) Diagonal fraction per bin ---
    TH1D* h_diagFrac = new TH1D("h_diag_frac",
        ";Truth bin index;Diagonal fraction",
        nBins, 0.5, nBins + 0.5);
    for (int i = 1; i <= nBins; ++i) {
        double rowSum = 0.0;
        for (int j = 1; j <= nBins; ++j) rowSum += h_raw->GetBinContent(i, j);
        double diagVal = h_raw->GetBinContent(i, i);
        if (rowSum > 0.0) {
            h_diagFrac->SetBinContent(i, diagVal / rowSum);
        }
    }

    TCanvas* c2 = new TCanvas("c_diagfrac", "", 1200, 600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c2->SetLeftMargin(0.10);
    c2->SetBottomMargin(0.13);
    c2->SetTopMargin(0.08);
    c2->SetRightMargin(0.05);

    h_diagFrac->SetFillColor(kAzure + 1);
    h_diagFrac->SetLineColor(kAzure + 3);
    h_diagFrac->GetYaxis()->SetRangeUser(0.0, 1.05);
    h_diagFrac->GetXaxis()->SetTitleSize(0.045);
    h_diagFrac->GetYaxis()->SetTitleSize(0.045);
    h_diagFrac->Draw("bar");

    TLine* half = new TLine(0.5, 0.5, nBins + 0.5, 0.5);
    half->SetLineColor(kRed);
    half->SetLineStyle(2);
    half->SetLineWidth(2);
    half->Draw("same");

    DrawSimLabels(inputFile);
    c2->Update();
    SaveCanvas(c2, "figs/diffractive/response/response_3d_diagonal_fraction.png");

    delete h_diagFrac; delete h_rowNorm; delete half; delete c2;
}
