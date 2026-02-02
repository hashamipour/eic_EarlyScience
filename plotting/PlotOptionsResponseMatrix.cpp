#include "Plotting.hpp"
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TText.h>
#include <iostream>
#include "Utility.hpp"

PlotOptionsResponseMatrix::PlotOptionsResponseMatrix(const TString& histName,
                                                     const char* xLabel,
                                                     const char* yLabel,
                                                     const char* saveName,
                                                     const bool  isLogX,
                                                     const bool  isLogY,
                                                     const std::pair<double, double>& xRange,
                                                     const std::pair<double, double>& yRange)
    : m_histName(histName),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName) {
        this->m_isLogX = isLogX;
        this->m_isLogY = isLogY;
        this->m_xRange = xRange;
        this->m_yRange = yRange;
      }

void PlotOptionsResponseMatrix::Plot(TFile* inputFile) {
    TH2D* h_matrix_orig = (TH2D*)inputFile->Get(m_histName);
    if (!h_matrix_orig) {
        std::cerr << "Error: 2D Histogram " << m_histName << " not found." << std::endl;
        return;
    }

    TH2D* h_matrix_perc = (TH2D*)h_matrix_orig->Clone(Form("%s_percentages", h_matrix_orig->GetName()));
    h_matrix_perc->Reset();
    h_matrix_perc->SetTitle("");
    h_matrix_perc->GetZaxis()->SetTitle("Percentage (%)");
    h_matrix_perc->GetZaxis()->SetRangeUser(0, 100);

    TCanvas* c = new TCanvas("c_response_matrix", "Response Matrix", 1200, 1000);
    gStyle->SetPaintTextFormat(".0f");
    gStyle->SetOptStat(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);

    if (m_isLogX) {c->SetLogx();}
    if (m_isLogY) {c->SetLogy();}

    int nbinsX = h_matrix_orig->GetNbinsX();
    int nbinsY = h_matrix_orig->GetNbinsY();
    
    for (int ix = 1; ix <= nbinsX; ix++) {
        double total_true = 0.0;
        for (int iy = 1; iy <= nbinsY; iy++) {
            total_true += h_matrix_orig->GetBinContent(ix, iy);
        }

        if (total_true > 0) {
            for (int iy = 1; iy <= nbinsY; iy++) {
                double content = h_matrix_orig->GetBinContent(ix, iy);
                double percentage = (content / total_true) * 100.0;

                if (percentage >= 1.0) {
                    h_matrix_perc->SetBinContent(ix, iy, percentage);
                }
            }
        }
    }

    h_matrix_perc->GetXaxis()->SetTitle(m_xLabel);
    h_matrix_perc->GetYaxis()->SetTitle(m_yLabel);
    h_matrix_perc->GetXaxis()->SetTitleOffset(1.3);
    h_matrix_perc->GetYaxis()->SetTitleOffset(1.3);
    gPad->Clear();
    h_matrix_perc->Draw("COLZ");

    // Use base class ranges if set, otherwise fall back to constructor parameters
    auto xRange = m_rangeX ? *m_rangeX : m_xRange;
    auto yRange = m_rangeY ? *m_rangeY : m_yRange;

    if ((xRange.first != -999.) && (xRange.second != -999.)) {
        // h_matrix_perc->GetXaxis()->SetLimits(xRange.first, xRange.second);
        // Logger::debug(Form("Set X axis range to [%.2f, %.2f]", xRange.first, xRange.second));
        h_matrix_perc->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
    }
    if ((yRange.first != -999.) && (yRange.second != -999.)) {
        // h_matrix_perc->GetYaxis()->SetLimits(yRange.first, yRange.second);
        h_matrix_perc->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
    }

    double xmin = h_matrix_perc->GetXaxis()->GetBinLowEdge(h_matrix_perc->GetXaxis()->GetFirst());
    double xmax = h_matrix_perc->GetXaxis()->GetBinUpEdge(h_matrix_perc->GetXaxis()->GetLast());
    double ymin = h_matrix_perc->GetYaxis()->GetBinLowEdge(h_matrix_perc->GetYaxis()->GetFirst());
    double ymax = h_matrix_perc->GetYaxis()->GetBinUpEdge(h_matrix_perc->GetYaxis()->GetLast());

    // std::cout << "Drawing diagonal from (" << xmin << ", " << ymin 
        //   << ") to (" << xmax << ", " << ymax << ")" << std::endl;
    TLine* diagonal_line = new TLine(xmin, ymin, xmax, ymax);

    diagonal_line->SetLineColor(kBlue);
    diagonal_line->SetLineStyle(2);
    diagonal_line->SetLineWidth(2);
    diagonal_line->Draw("SAME");

    TText t;
    t.SetTextSize(0.02);
    t.SetTextAlign(22);
    for (int iy = 1; iy <= nbinsY; iy++) {
        for (int ix = 1; ix <= nbinsX; ix++) {
            double bin_content = h_matrix_perc->GetBinContent(ix, iy);
            double bin_x = h_matrix_perc->GetXaxis()->GetBinCenter(ix);
            double bin_y = h_matrix_perc->GetYaxis()->GetBinCenter(iy);
            
            if (bin_content > 1.0) {
                t.SetTextColor(kBlack);
                t.DrawText(bin_x, bin_y, Form("%.0f", bin_content));
            } else if (bin_content > 0.1) {
                t.SetTextColor(kBlue);
                t.DrawText(bin_x, bin_y, "<1%");
            }
        }
    }

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.2, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} #10x100 GeV");

    SetCustomPalette("SolarBloom");
    TColor::InvertPalette(); // Invert the palette for better visibility
    c->Update();
    SaveCanvas(c, m_saveName);

    delete h_matrix_perc;
    delete c;
}
