#include "Plotting.hpp"

#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TPad.h>

#include <TF1.h>
#include <TLatex.h>

#include <cstring>   // std::strchr
#include <iostream>
#include <utility>

// Constructor must match your header
PlotOptions1D::PlotOptions1D(const std::vector<TString>& histNames,
                             const std::vector<const char*>& legendEntries,
                             const std::vector<const char*>& drawOptions,
                             const char* canvasTitle,
                             const char* xLabel,
                             const char* yLabel,
                             const char* saveName,
                             const bool  isLogX,
                             const bool  isLogY,
                             const bool  normalizeToPDF)
    : m_histNames(histNames),
      m_legendEntries(legendEntries),
      m_drawOptions(drawOptions),
      m_canvasTitle(canvasTitle),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_saveName(saveName),
      m_isLogX(isLogX),
      m_isLogY(isLogY),
      m_normalizeToPDF(normalizeToPDF)
{}

void PlotOptions1D::Plot(TFile* inputFile) {
    // unique canvas name to avoid clashes
    TCanvas* c = new TCanvas(Form("c_%p", this), m_canvasTitle, 1200, 900);

    if (m_isLogX) c->SetLogx();
    if (m_isLogY) c->SetLogy();

    
    // ---------- Draw histograms ----------
    TH1* firstHist = nullptr;

    for (size_t i = 0; i < m_histNames.size(); ++i) {
        TH1* hist = static_cast<TH1*>(inputFile->Get(m_histNames[i]));
        if (!hist) {
            std::cerr << "Error: Histogram " << m_histNames[i] << " not found." << std::endl;
            continue;
        }

        hist->SetStats(false);

        if (m_normalizeToPDF) {
            hist->Sumw2();
            const double integral_w = hist->Integral("width");
            if (integral_w > 0.0) hist->Scale(1.0 / integral_w);
        }

        // simple styling by name convention (adapt as you like)
        hist->SetLineWidth((i == 0) ? 2 : 1);
        if (m_histNames[i].Contains("truth")|| m_histNames[i].Contains("MC")) {
            hist->SetLineColor(kBlack);
        } else if (m_histNames[i].Contains("EM")|| m_histNames[i].Contains("B0")) {
            hist->SetLineColor(kRed);
            hist->SetMarkerColor(kRed);
            hist->SetMarkerStyle(20);
        } else if (m_histNames[i].Contains("DA")|| m_histNames[i].Contains("RP")) {
            hist->SetLineColor(kBlue);
            hist->SetMarkerColor(kBlue);
            hist->SetMarkerStyle(20);
        } else if (m_histNames[i].Contains("Sigma")) {
            hist->SetLineColor(kGreen + 2);
            hist->SetMarkerColor(kGreen + 2);
            hist->SetMarkerStyle(20);
        } else if (m_histNames[i].Contains("Sum")) {
            hist->SetLineColor(kOrange+7);
            hist->SetMarkerColor(kOrange+7);
            hist->SetMarkerStyle(24);  // Hollow circle
        }

        // Draw with user draw options; ensure SAME for i>0
        TString drawOption = m_drawOptions[i];
        if (i == 0) {
            hist->Draw(drawOption);
            firstHist = hist;
        } else {
            if (!drawOption.Contains("SAME")) drawOption += " SAME";
            hist->Draw(drawOption);
        }

        // Legend entry: use legend styles ("l","p","lp"), not draw options
        // const char* opt = m_drawOptions[i];
        // const bool hasMarkers = (hist->GetMarkerStyle() > 0) ||
        //                         (opt && (std::strchr(opt, 'P') || std::strchr(opt, 'p')));
        // legend->AddEntry(hist, m_legendEntries[i], hasMarkers ? "lp" : "l");
    }

    // ---------- Legend placement (plot-area [0..1] -> pad NDC) ----------
    c->cd();
    gPad->Update(); // margins are known after Update()

    // Log the requested (plot-area) legend box
    // Logger::debug("Legend position (plot-area fractions): (" +
    //               std::to_string(m_legendLB.first) + ", " +
    //               std::to_string(m_legendLB.second) + ") - (" +
    //               std::to_string(m_legendRT.first) + ", " +
    //               std::to_string(m_legendRT.second) + ")");

    // Create legend with DIRECT NDC coordinates - no mapping
    TLegend* legend = nullptr;
    if (m_legendLB && m_legendRT) {
        // Logger::debug("Using user-defined legend position.");
        legend = new TLegend(m_legendLB->first, m_legendLB->second,
                                  m_legendRT->first, m_legendRT->second);
        // Logger::debug("Legend position set to: (" +
        //               std::to_string(m_legendLB->first) + ", " +
        //               std::to_string(m_legendLB->second) + ") - (" +
        //               std::to_string(m_legendRT->first) + ", " +
        //               std::to_string(m_legendRT->second) + ")");
    } else {
        // Default legend position if not set by user
        legend = new TLegend(0.7, 0.7, 0.98, 0.9); 
        // Logger::debug("Using explicit default legend position (0.7, 0.7) - (0.98, 0.9).");
    }
    
    
    for (size_t i = 0; i < m_histNames.size(); ++i) {
        TH1* hist = static_cast<TH1*>(inputFile->Get(m_histNames[i]));
        if (hist) legend->AddEntry(hist, m_legendEntries[i], "lp");
    }
    
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);

    // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    legend->Draw();
    gPad->Update(); // Update the pad after drawing the legend to finalize its position
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    // ---------- Axes labels and ranges ----------
    if (firstHist) {
        firstHist->GetXaxis()->SetTitle(m_xLabel);
        firstHist->GetYaxis()->SetTitle(m_yLabel);
        firstHist->GetXaxis()->SetTitleSize(0.04);
        firstHist->GetYaxis()->SetTitleSize(0.04);
        firstHist->GetXaxis()->SetTitleOffset(1.0);
        firstHist->GetYaxis()->SetTitleOffset(1.0);

        // Apply custom ranges if set
        if (m_rangeX) {
            firstHist->GetXaxis()->SetRangeUser(m_rangeX->first, m_rangeX->second);
        }
        if (m_rangeY) {
            firstHist->GetYaxis()->SetRangeUser(m_rangeY->first, m_rangeY->second);
        }
    }

    // ---------- Add ePIC simulation labels ----------
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.2, 0.92, "#bf{ePIC} Simulation (100k events)");
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");

    // ---------- Draw legend and save ----------
    // legend->Draw();
    // gPad->Update();
    c->Update();
    SaveCanvas(c, m_saveName);

    delete c;
}
