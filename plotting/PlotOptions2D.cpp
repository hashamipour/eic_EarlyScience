#include "Plotting.hpp"
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>

PlotOptions2D::PlotOptions2D(const TString& histName,
                             const char* xLabel,
                             const char* yLabel,
                             const char* saveName,
                             const bool isLogX,
                             const bool isLogY,
                             const std::pair<double, double>& xRange,
                             const std::pair<double, double>& yRange)
    : m_histName(histName),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName),
      m_isLogX(isLogX),
      m_isLogY(isLogY),
      m_xRange(xRange),
      m_yRange(yRange) {}

void PlotOptions2D::Plot(TFile* inputFile) {
    TH2D* hist = (TH2D*)inputFile->Get(m_histName);
    if (!hist) {
        std::cerr << "Error: 2D Histogram " << m_histName << " not found." << std::endl;
        return;
    }

    TCanvas* c = new TCanvas("c_density", "2D Density", 1200, 1000);
    gStyle->SetOptStat(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);

    if (m_isLogX) c->SetLogx();
    if (m_isLogY) c->SetLogy();

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(m_xLabel);
    hist->GetYaxis()->SetTitle(m_yLabel);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3);

    auto xRange = m_rangeX ? *m_rangeX : m_xRange;
    auto yRange = m_rangeY ? *m_rangeY : m_yRange;
    if (xRange.first != -999. && xRange.second != -999.) {
        hist->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
    }
    if (yRange.first != -999. && yRange.second != -999.) {
        hist->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
    }

    hist->Draw("COLZ");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.2, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

    c->Update();
    SaveCanvas(c, m_saveName);

    delete c;
}
