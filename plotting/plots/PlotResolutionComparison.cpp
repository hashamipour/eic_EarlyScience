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

void PlotResolutionComparison(TFile* inputFile,
                                     const std::vector<const char*>& histNames,
                                     const std::vector<const char*>& labels,
                                     const char* xTitle,
                                     const char* yTitle,
                                     double xLo, double xHi,
                                     double yLo, double yHi,
                                     bool logX,
                                     const char* saveName,
                                     bool disableFit) {
    if (!inputFile) return;
    const int nMethods = static_cast<int>(histNames.size());

    // Distinct marker styles per method
    const Style_t markerStyles[] = {20, 21, 22, 23, 33, 34};
    // Horizontal offsets (multiplicative for logX, additive otherwise)
    const double hFrac[] = {-0.07, +0.07, 0.0, -0.05, +0.05, 0.0};

    std::vector<TGraphErrors*> gFit(nMethods, nullptr);
    std::vector<TGraphErrors*> gRMS(nMethods, nullptr);
    bool anyValid = false;

    for (int m = 0; m < nMethods; ++m) {
        TH2D* h2 = (TH2D*)inputFile->Get(histNames[m]);
        if (!h2) {
            Logger::warning(std::string(histNames[m]) + " not found; skipping.");
            continue;
        }
        const int nbinsX = h2->GetNbinsX();
        gFit[m] = new TGraphErrors();
        gRMS[m] = new TGraphErrors();
        int npFit = 0, npRMS = 0;
        for (int j = 1; j <= nbinsX; ++j) {
            TH1D* proj = h2->ProjectionY(Form("_comp_%s_%d", histNames[m], j), j, j);
            if (proj->GetEntries() < 20) { delete proj; continue; }

            double xc = h2->GetXaxis()->GetBinCenter(j);
            double histMean = proj->GetMean();
            double histRMS  = proj->GetRMS();

            // RMS point (always)
            double xRMS = logX ? xc * (1.0 + hFrac[m % 6] + 0.03)
                               : xc + ((m % 6)-3)*0.005  ;
            gRMS[m]->SetPoint(npRMS, xRMS, histMean);
            gRMS[m]->SetPointError(npRMS, 0.0, histRMS);
            ++npRMS;

            if (!disableFit) {
                // Fitted point: use Gaussian fit with moment-based sigma
                // (consistent with PlotOptionsBinnedRelRes::Plot)
                double fitLo = histMean - 1.5 * histRMS;
                double fitHi = histMean + 1.5 * histRMS;
                TF1* fit = new TF1(Form("gfit_comp_%s_%d", histNames[m], j),
                                   "gaus", fitLo, fitHi);
                fit->SetParameters(proj->GetMaximum(), histMean, histRMS);
                proj->Fit(fit, "RQ");

                double fitMean = fit->Mean(fitLo, fitHi);
                double fitVar  = fit->Variance(fitLo, fitHi);
                double fitSigma = (fitVar > 0.0) ? std::sqrt(fitVar) : -1.0;
                if (!std::isfinite(fitMean))  fitMean  = fit->GetParameter(1);
                if (!std::isfinite(fitSigma) || fitSigma < 0.0)
                    fitSigma = std::abs(fit->GetParameter(2));

                double xFit = logX ? xc * (1.0 + hFrac[m % 6])
                                   : xc + ((m % 6)-3)*0.01;
                gFit[m]->SetPoint(npFit, xFit, fitMean);
                gFit[m]->SetPointError(npFit, 0.0, fitSigma);
                ++npFit;

                delete fit;
            }
            delete proj;
        }
        if (npFit > 0 || npRMS > 0) anyValid = true;
    }

    if (!anyValid) {
        for (int m = 0; m < nMethods; ++m) { delete gFit[m]; delete gRMS[m]; }
        return;
    }

    TCanvas* c = new TCanvas(Form("c_comp_%s", saveName), "", 1600, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.15);
    if (logX) c->SetLogx();

    // Draw fitted graphs first (black, filled markers) — skip if fits disabled
    bool first = true;
    if (!disableFit) {
        for (int m = 0; m < nMethods; ++m) {
            if (!gFit[m] || gFit[m]->GetN() == 0) continue;
            gFit[m]->SetMarkerStyle(markerStyles[m % 6]);
            gFit[m]->SetMarkerColor(kBlack);
            gFit[m]->SetLineColor(kBlack);
            gFit[m]->SetLineWidth(2);
            if (first) {
                gFit[m]->SetTitle("");
                gFit[m]->GetXaxis()->SetTitle(xTitle);
                gFit[m]->GetYaxis()->SetTitle(yTitle);
                gFit[m]->GetXaxis()->SetLimits(xLo, xHi);
                gFit[m]->GetYaxis()->SetRangeUser(yLo, yHi);
                gFit[m]->GetXaxis()->SetTitleSize(0.045);
                gFit[m]->GetYaxis()->SetTitleSize(0.045);
                gFit[m]->GetXaxis()->SetTitleOffset(1.2);
                gFit[m]->GetYaxis()->SetTitleOffset(1.2);
                gFit[m]->Draw("AP");
                first = false;
            } else {
                gFit[m]->Draw("P SAME");
            }
        }
    }

    // Draw RMS graphs
    const Style_t openMarkers[] = {24, 25, 26, 32, 27, 28};
    const Color_t methodColors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
    for (int m = 0; m < nMethods; ++m) {
        if (!gRMS[m] || gRMS[m]->GetN() == 0) continue;
        gRMS[m]->SetMarkerStyle(disableFit ? markerStyles[m % 6] : openMarkers[m % 6]);
        gRMS[m]->SetMarkerColor(disableFit ? methodColors[m % 6] : kBlue);
        gRMS[m]->SetLineColor(disableFit ? methodColors[m % 6] : kBlue);
        gRMS[m]->SetLineWidth(2);
        if (first) {
            gRMS[m]->SetTitle("");
            gRMS[m]->GetXaxis()->SetTitle(xTitle);
            gRMS[m]->GetYaxis()->SetTitle(yTitle);
            gRMS[m]->GetXaxis()->SetLimits(xLo, xHi);
            gRMS[m]->GetYaxis()->SetRangeUser(yLo, yHi);
            gRMS[m]->GetXaxis()->SetTitleSize(0.045);
            gRMS[m]->GetYaxis()->SetTitleSize(0.045);
            gRMS[m]->GetXaxis()->SetTitleOffset(1.2);
            gRMS[m]->GetYaxis()->SetTitleOffset(1.2);
            gRMS[m]->Draw("AP");
            first = false;
        } else {
            gRMS[m]->Draw("P SAME");
        }
    }

    TLine* zero = new TLine(xLo, 0.0, xHi, 0.0);
    zero->SetLineColor(kGray + 1);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw();

    // Legend
    TLegend* leg = new TLegend(0.15, 0.60, 0.50, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    if (!disableFit) {
        leg->AddEntry((TObject*)nullptr, "#bf{Black: Fitted (#mu #pm #sigma_{fit})}", "");
        leg->AddEntry((TObject*)nullptr, "#bf{#color[4]{Blue: RMS (#mu #pm #sigma_{RMS})}}", "");
    }
    for (int m = 0; m < nMethods; ++m) {
        if (!disableFit && gFit[m] && gFit[m]->GetN() > 0)
            leg->AddEntry(gFit[m], Form("%s (Fitted)", labels[m]), "ep");
        if (gRMS[m] && gRMS[m]->GetN() > 0)
            leg->AddEntry(gRMS[m], disableFit ? labels[m] : Form("%s (RMS)", labels[m]), "ep");
    }
    leg->Draw();

    DrawSimLabels(inputFile);
    c->Update();
    SaveCanvas(c, saveName);

    for (int m = 0; m < nMethods; ++m) { delete gFit[m]; delete gRMS[m]; }
    delete zero; delete leg; delete c;
}
