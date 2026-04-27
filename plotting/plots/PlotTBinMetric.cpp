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

static std::pair<double, double> ComputeMetricFromCounts(const double nGen,
                                                         const double nMeas,
                                                         const double nSame,
                                                         const TBinMetricKind kind) {
    double value = 0.0;
    double error = 0.0;

    if (kind == TBinMetricKind::Acceptance) {
        if (nGen > 0.0) {
            value = nMeas / nGen;
            if (nMeas > 0.0) {
                error = std::abs(value) * std::sqrt((1.0 / nMeas) + (1.0 / nGen));
            }
        }
        return {value, error};
    }

    if (kind == TBinMetricKind::Purity) {
        if (nMeas > 0.0) {
            value = nSame / nMeas;
            if (nSame > 0.0) {
                error = std::abs(value) * std::sqrt((1.0 / nSame) + (1.0 / nMeas));
            }
        }
        return {value, error};
    }

    // Efficiency
    if (nGen > 0.0) {
        value = nSame / nGen;
        if (nSame > 0.0) {
            error = std::abs(value) * std::sqrt((1.0 / nSame) + (1.0 / nGen));
        }
    }
    return {value, error};
}

static bool BuildTMetricHistograms(TFile* inputFile,
                                   const TBinMetricKind kind,
                                   TH1D*& hB0Out,
                                   TH1D*& hRPOut,
                                   TH1D*& hSumOut) {
    hB0Out = nullptr;
    hRPOut = nullptr;
    hSumOut = nullptr;
    if (!inputFile) return false;

    TH1D* hTruthB0 = (TH1D*)inputFile->Get("t_truth_mc_B0");
    TH1D* hRecoB0 = (TH1D*)inputFile->Get("t_reco_mc_B0");
    TH2D* hCorrB0 = (TH2D*)inputFile->Get("t_corr_B0");
    TH1D* hTruthRP = (TH1D*)inputFile->Get("t_truth_mc_RP");
    TH1D* hRecoRP = (TH1D*)inputFile->Get("t_reco_mc_RP");
    TH2D* hCorrRP = (TH2D*)inputFile->Get("t_corr_RP");
    if (!hTruthB0 || !hRecoB0 || !hCorrB0 || !hTruthRP || !hRecoRP || !hCorrRP) {
        Logger::warning("Missing |t| histograms for metric plots. "
                        "(need t_truth_mc_{B0,RP}, t_reco_mc_{B0,RP}, t_corr_{B0,RP})");
        return false;
    }

    hB0Out = (TH1D*)hTruthB0->Clone("t_metric_b0_tmp");
    hRPOut = (TH1D*)hTruthRP->Clone("t_metric_rp_tmp");
    hSumOut = (TH1D*)hTruthB0->Clone("t_metric_sum_tmp");
    hB0Out->SetDirectory(nullptr);
    hRPOut->SetDirectory(nullptr);
    hSumOut->SetDirectory(nullptr);
    hB0Out->Reset("ICES");
    hRPOut->Reset("ICES");
    hSumOut->Reset("ICES");

    const int nbins = hTruthB0->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        const double genB0 = hTruthB0->GetBinContent(i);
        const double measB0 = hRecoB0->GetBinContent(i);
        const double sameB0 = hCorrB0->GetBinContent(i, i);

        const double genRP = hTruthRP->GetBinContent(i);
        const double measRP = hRecoRP->GetBinContent(i);
        const double sameRP = hCorrRP->GetBinContent(i, i);

        const auto b0Metric = ComputeMetricFromCounts(genB0, measB0, sameB0, kind);
        const auto rpMetric = ComputeMetricFromCounts(genRP, measRP, sameRP, kind);
        const auto sumMetric = ComputeMetricFromCounts(genB0 + genRP, measB0 + measRP, sameB0 + sameRP, kind);

        hB0Out->SetBinContent(i, b0Metric.first);
        hB0Out->SetBinError(i, b0Metric.second);
        hRPOut->SetBinContent(i, rpMetric.first);
        hRPOut->SetBinError(i, rpMetric.second);
        hSumOut->SetBinContent(i, sumMetric.first);
        hSumOut->SetBinError(i, sumMetric.second);
    }

    return true;
}

void PlotTBinMetric(TFile* inputFile,
                           const TBinMetricKind kind,
                           const char* title,
                           const char* yLabel,
                           const char* saveName) {
    TH1D* hB0 = nullptr;
    TH1D* hRP = nullptr;
    TH1D* hSum = nullptr;
    if (!BuildTMetricHistograms(inputFile, kind, hB0, hRP, hSum)) {
        return;
    }

    TCanvas* c = new TCanvas(Form("c_%s", saveName), title, 1200, 900);
    c->SetLogx();
    c->SetGridx();
    c->SetGridy();

    hSum->SetTitle("");
    hSum->GetXaxis()->SetTitle("|t| [GeV^{2}]");
    hSum->GetYaxis()->SetTitle(yLabel);
    hSum->GetXaxis()->SetTitleOffset(1.1);
    hSum->GetYaxis()->SetTitleOffset(1.2);
    hSum->GetXaxis()->SetMoreLogLabels();
    hSum->GetXaxis()->SetNoExponent();

    hB0->SetMarkerStyle(20);
    hB0->SetMarkerSize(1.0);
    hB0->SetMarkerColor(kRed + 1);
    hB0->SetLineColor(kRed + 1);
    hB0->SetLineWidth(2);

    hRP->SetMarkerStyle(21);
    hRP->SetMarkerSize(1.0);
    hRP->SetMarkerColor(kBlue + 1);
    hRP->SetLineColor(kBlue + 1);
    hRP->SetLineWidth(2);

    hSum->SetMarkerStyle(24);
    hSum->SetMarkerSize(1.0);
    hSum->SetMarkerColor(kBlack);
    hSum->SetLineColor(kBlack);
    hSum->SetLineWidth(2);

    double ymax = 0.0;
    for (int i = 1; i <= hSum->GetNbinsX(); ++i) {
        ymax = std::max(ymax, hB0->GetBinContent(i) + hB0->GetBinError(i));
        ymax = std::max(ymax, hRP->GetBinContent(i) + hRP->GetBinError(i));
        ymax = std::max(ymax, hSum->GetBinContent(i) + hSum->GetBinError(i));
    }
    if (ymax <= 0.0) ymax = 1.0;
    const double minY = 0.0;
    const double floorMax = (kind == TBinMetricKind::Acceptance) ? 1.2 : 1.05;
    hSum->GetYaxis()->SetRangeUser(minY, std::max(ymax * 1.25, floorMax));

    hSum->Draw("E1");
    hB0->Draw("E1 SAME");
    hRP->Draw("E1 SAME");

    const double xMin = hSum->GetXaxis()->GetXmin();
    const double xMax = hSum->GetXaxis()->GetXmax();
    TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
    unity->SetLineColor(kGray + 2);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->Draw("SAME");

    TLegend* leg = new TLegend(0.66, 0.72, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hB0, "B0", "lep");
    leg->AddEntry(hRP, "RP", "lep");
    leg->AddEntry(hSum, "B0+RP", "lep");
    leg->Draw();

    DrawSimLabels(inputFile);
    SaveCanvas(c, saveName);

    delete leg;
    delete unity;
    delete c;
    delete hB0;
    delete hRP;
    delete hSum;
}
