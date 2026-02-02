#include "Plotting.hpp"
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <iostream>
#include <algorithm>
#include "Utility.hpp"

PlotOptionsCombinedCorrelation::PlotOptionsCombinedCorrelation(const std::vector<TString>& histNames,
                                                               const std::vector<const char*>& legendEntries,
                                                               const std::vector<Color_t>& colors,
                                                               const std::vector<Style_t>& markerStyles,
                                                               const char* canvasTitle,
                                                               const char* xLabel,
                                                               const char* yLabel,
                                                               const char* saveName,
                                                               const std::pair<double, double>& xRange,
                                                               const std::pair<double, double>& yRange,
                                                               const bool isLogX,
                                                               const bool isLogY)
    : m_histNames(histNames),
      m_legendEntries(legendEntries),
      m_colors(colors),
      m_markerStyles(markerStyles),
      m_canvasTitle(canvasTitle),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName),
      m_xRange(xRange),
      m_yRange(yRange),
      m_isLogX(isLogX),
      m_isLogY(isLogY) {}

void PlotOptionsCombinedCorrelation::Plot(TFile* inputFile) {
    TCanvas* c = new TCanvas("c_combined_corr", m_canvasTitle, 1200, 1000);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    c->SetLogx(m_isLogX);
    c->SetLogy(m_isLogY);

    gStyle->SetOptStat(0);

    auto xRange = m_rangeX ? *m_rangeX : m_xRange;
    auto yRange = m_rangeY ? *m_rangeY : m_yRange;

    double xmin = xRange.first;
    double xmax = xRange.second;
    double ymin = yRange.first;
    double ymax = yRange.second;

    if (xmin == -999. || xmax == -999. || ymin == -999. || ymax == -999.) {
        xmin = 1.0e9; xmax = -1.0;
        ymin = 1.0e9; ymax = -1.0;
        for (size_t i = 0; i < m_histNames.size(); ++i) {
            TGraph* graph = (TGraph*)inputFile->Get(m_histNames[i]);
            if (!graph) continue;
            const int n = graph->GetN();
            for (int j = 0; j < n; ++j) {
                double x = 0.0, y = 0.0;
                graph->GetPoint(j, x, y);
                if (m_isLogX && x <= 0.0) continue;
                if (m_isLogY && y <= 0.0) continue;
                xmin = std::min(xmin, x);
                xmax = std::max(xmax, x);
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
            }
        }
        if (xmin >= xmax || ymin >= ymax) {
            xmin = 1.0; xmax = 10.0;
            ymin = 1.0; ymax = 10.0;
        }
    }

    const int nBinsX = 100;
    const int nBinsY = 100;
    std::vector<Double_t> xBins = m_isLogX ? GetLogBins(xmin, xmax, nBinsX) : std::vector<Double_t>();
    std::vector<Double_t> yBins = m_isLogY ? GetLogBins(ymin, ymax, nBinsY) : std::vector<Double_t>();

    TH2D* density = nullptr;
    if (m_isLogX && m_isLogY) {
        density = new TH2D("corr_density", "", xBins.size()-1, xBins.data(), yBins.size()-1, yBins.data());
    } else if (m_isLogX && !m_isLogY) {
        density = new TH2D("corr_density", "", xBins.size()-1, xBins.data(), nBinsY, ymin, ymax);
    } else if (!m_isLogX && m_isLogY) {
        density = new TH2D("corr_density", "", nBinsX, xmin, xmax, yBins.size()-1, yBins.data());
    } else {
        density = new TH2D("corr_density", "", nBinsX, xmin, xmax, nBinsY, ymin, ymax);
    }

    for (size_t i = 0; i < m_histNames.size(); ++i) {
        TGraph* graph = (TGraph*)inputFile->Get(m_histNames[i]);
        if (!graph) {
            std::cerr << "Error: TGraph " << m_histNames[i] << " not found." << std::endl;
            continue;
        }
        const int n = graph->GetN();
        for (int j = 0; j < n; ++j) {
            double x = 0.0, y = 0.0;
            graph->GetPoint(j, x, y);
            if (m_isLogX && x <= 0.0) continue;
            if (m_isLogY && y <= 0.0) continue;
            density->Fill(x, y);
        }
    }

    density->SetTitle("");
    density->GetXaxis()->SetTitle(m_xLabel);
    density->GetYaxis()->SetTitle(m_yLabel);
    density->GetXaxis()->SetTitleOffset(1.3);
    density->GetYaxis()->SetTitleOffset(1.3);
    density->Draw("COLZ");

    const double line_min = std::max(xmin, ymin);
    const double line_max = std::min(xmax, ymax);
    if (line_max > line_min) {
        TLine* diagLine = new TLine(line_min, line_min, line_max, line_max);
        diagLine->SetLineColor(kBlack);
        diagLine->SetLineWidth(2);
        diagLine->SetLineStyle(2);
        diagLine->Draw("same");
    }

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.17, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");
    
    c->Update();
    SaveCanvas(c, m_saveName);
    
    delete density;
    delete c;
}
