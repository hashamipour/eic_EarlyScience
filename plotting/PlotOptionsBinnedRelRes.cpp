#include "Plotting.hpp"
#include <TF1.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <iostream>

PlotOptionsBinnedRelRes::PlotOptionsBinnedRelRes(const TString& histName,
                                                 const char* title,
                                                 const char* xLabel,
                                                 const char* yLabel,
                                                 const std::vector<std::pair<double, double>>& fitRanges,
                                                 const char* saveName,
                                                 const char* binSavePrefix,
                                                 const std::pair<double, double>& x_axis_range,
                                                 const bool isLogX
                                                )
    : m_histName(histName),
      m_title(title),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_fitRanges(fitRanges),
      m_xMinFit(0.0),
      m_xMaxFit(0.0),
      m_saveName(saveName),
      m_binSavePrefix(binSavePrefix),
      m_xAxisRange(x_axis_range),
      m_isLogX(isLogX) {}

void PlotOptionsBinnedRelRes::SetFitRangeByBins(TH1D* hist) {
    int peakBin = hist->GetMaximumBin();
    double peakCenter = hist->GetBinCenter(peakBin);

    struct FitResult {
        int nBinsLeft;
        int nBinsRight;
        int totalBins;
        double chi2_ndf;
        double yValueDiffPercent;
        double xMinFit;
        double xMaxFit;
    };

    std::vector<FitResult> validFits;

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "   nBins_L | nBins_R | Total | Chi2/NDF | |f(x_L) - f(x_R)|% " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    for (int nBinsLeft = 2; nBinsLeft <= 3; ++nBinsLeft) {
        for (int nBinsRight = 2; nBinsRight <= 2; ++nBinsRight) {
            
            int minBinFit = peakBin - nBinsLeft;
            int maxBinFit = peakBin + nBinsRight;
            
            if (minBinFit < 1 || maxBinFit > hist->GetNbinsX()) {
                continue;
            }
            
            double x_min_fit = hist->GetBinLowEdge(minBinFit);
            double x_max_fit = hist->GetBinLowEdge(maxBinFit + 1);
            
            double x_left_bin_center = hist->GetBinCenter(minBinFit);
            double x_right_bin_center = hist->GetBinCenter(maxBinFit);
            
            TF1* currentGaus = new TF1("currentGaus", "gaus", x_min_fit, x_max_fit);
            currentGaus->SetParameters(hist->GetMaximum(), hist->GetBinCenter(peakBin), hist->GetRMS());
            
            hist->Fit(currentGaus, "RQN", "", x_min_fit, x_max_fit);
            
            double currentChi2 = currentGaus->GetChisquare();
            int ndf = currentGaus->GetNDF();
            
            if (ndf > 0) {
                double chi2_ndf = currentChi2 / ndf;
                
                if (chi2_ndf > 10.0) {
                    delete currentGaus;
                    continue;
                }
                
                double y_left = currentGaus->Eval(x_left_bin_center);
                double y_right = currentGaus->Eval(x_right_bin_center);
                double yValueDiff = TMath::Abs(y_left - y_right);
                double y_max = TMath::Max(y_left, y_right);
                double yValueDiffPercent = (y_max > 0) ? (yValueDiff / y_max * 100.0) : 0.0;
                
                if (yValueDiffPercent > 10.0) {
                    delete currentGaus;
                    continue;
                }
                
                FitResult result;
                result.nBinsLeft = nBinsLeft;
                result.nBinsRight = nBinsRight;
                result.totalBins = nBinsLeft + nBinsRight;
                result.chi2_ndf = chi2_ndf;
                result.yValueDiffPercent = yValueDiffPercent;
                result.xMinFit = x_min_fit;
                result.xMaxFit = x_max_fit;
                
                validFits.push_back(result);
                
                std::cout << TString::Format("   %7d | %7d | %5d | %8.4f | %8.2f%%", 
                    nBinsLeft, nBinsRight, result.totalBins, chi2_ndf, yValueDiffPercent) << std::endl;
            }
            
            delete currentGaus;
        }
    }
    
    if (validFits.empty()) {
        std::cerr << "Error: No valid fits found! (All fits had chi2/ndf > 10 or y-difference > 10%)" << std::endl;
        m_xMinFit = hist->GetBinCenter(peakBin) - 0.55 * hist->GetRMS();
        m_xMaxFit = hist->GetBinCenter(peakBin) + 0.45 * hist->GetRMS();
        return;
    }
    
    std::sort(validFits.begin(), validFits.end(), [](const FitResult& a, const FitResult& b) {
        if (TMath::Abs(a.yValueDiffPercent - b.yValueDiffPercent) > 0.01) 
            return a.yValueDiffPercent < b.yValueDiffPercent;
        if (a.totalBins != b.totalBins) 
            return a.totalBins > b.totalBins;
        return a.chi2_ndf < b.chi2_ndf;
    });
    
    double bestYSymmetry = validFits[0].yValueDiffPercent;
    
    std::vector<FitResult> bestYSymmetryFits;
    for (const auto& fit : validFits) {
        if (fit.yValueDiffPercent <= bestYSymmetry + 0.1) {
            bestYSymmetryFits.push_back(fit);
        }
    }
    
    int maxSizeAmongBestY = 0;
    for (const auto& fit : bestYSymmetryFits) {
        if (fit.totalBins > maxSizeAmongBestY) {
            maxSizeAmongBestY = fit.totalBins;
        }
    }
    
    std::vector<FitResult> bestYAndSizeFits;
    for (const auto& fit : bestYSymmetryFits) {
        if (fit.totalBins == maxSizeAmongBestY) {
            bestYAndSizeFits.push_back(fit);
        }
    }
    
    FitResult bestFit = *std::min_element(bestYAndSizeFits.begin(), bestYAndSizeFits.end(),
        [](const FitResult& a, const FitResult& b) {
            return a.chi2_ndf < b.chi2_ndf;
        });
    
    m_xMinFit = bestFit.xMinFit;
    m_xMaxFit = bestFit.xMaxFit;
    
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "SELECTION SUMMARY:" << std::endl;
    std::cout << "Total valid fits (chi2<10, y-diff<10%): " << validFits.size() << std::endl;
    std::cout << "Best y-symmetry: " << bestYSymmetry << "%" << std::endl;
    std::cout << "Fits with best y-symmetry: " << bestYSymmetryFits.size() << std::endl;
    std::cout << "Max interval size among best y-symmetry: " << maxSizeAmongBestY << " bins" << std::endl;
    std::cout << "Fits with best y-symmetry AND max size: " << bestYAndSizeFits.size() << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "OPTIMAL FIT SELECTED:" << std::endl;
    std::cout << "Bins: " << bestFit.nBinsLeft << " left + " << bestFit.nBinsRight << " right = " 
              << bestFit.totalBins << " total" << std::endl;
    std::cout << "Fit range: [" << m_xMinFit << ", " << m_xMaxFit << "]" << std::endl;
    std::cout << "Chi2/NDF: " << bestFit.chi2_ndf << std::endl;
    std::cout << "Y-symmetry: " << bestFit.yValueDiffPercent << "%" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}

void PlotOptionsBinnedRelRes::Plot(TFile* inputFile) {
    TH2D* h_RelRes_binned = (TH2D*)inputFile->Get(m_histName);
    if (!h_RelRes_binned) {
        std::cerr << "Error: 2D Histogram " << m_histName << " not found." << std::endl;
        return;
    }
    const std::string simLabel = BuildSimLabel(inputFile);

    TCanvas* c = new TCanvas("c_binned", "", 1200, 800);
    c->SetLeftMargin(0.15);

    if (m_isLogX){c->SetLogx();}

    int nbinsX = h_RelRes_binned->GetNbinsX();
    TGraphErrors* g = new TGraphErrors(nbinsX);
    TGraphErrors* g_RMS = new TGraphErrors(nbinsX);

    for (int j = 1; j <= nbinsX; ++j) {
        TH1D* projY = h_RelRes_binned->ProjectionY("", j, j);
        
        if (projY->GetEntries() == 0) {
            delete projY;
            continue;
        }

        double _center = h_RelRes_binned->GetXaxis()->GetBinCenter(j);
        
        // Check if we should skip fitting for this bin
        bool skipFit = false;
        if (!m_fitRanges.empty() && j - 1 < static_cast<int>(m_fitRanges.size())) {
            if (TMath::AreEqualRel(m_fitRanges[j - 1].first, 0., 1e-6) && 
                TMath::AreEqualRel(m_fitRanges[j - 1].second, 0., 1e-6)) {
                skipFit = true;
            }
        }
        
        TF1* gaus = nullptr;
        double mean = 0.0;
        double sigma = 0.0;
        
        if (!skipFit) {
            // Determine fit range: use provided ranges or automatic detection
            if (m_fitRanges.empty()) {
                SetFitRangeByBins(projY);
            } else {
                if (j - 1 < static_cast<int>(m_fitRanges.size())) {
                    m_xMinFit = m_fitRanges[j - 1].first;
                    m_xMaxFit = m_fitRanges[j - 1].second;
                } else {
                    std::cerr << "Warning: Not enough fit ranges provided for bin " << j 
                              << ". Using automatic range detection." << std::endl;
                    SetFitRangeByBins(projY);
                }
            }
            
            gaus = new TF1("gaus", "gaus", m_xMinFit, m_xMaxFit);
            gaus->SetParameters(projY->GetMaximum(), projY->GetMean(), projY->GetRMS());
            projY->Fit(gaus, "RQ");
            
            mean = gaus->GetParameter(1);
            sigma = gaus->GetParameter(2);
        }
        
        // Always create and save the projection plot
        TCanvas* c_proj = new TCanvas(Form("c_proj_%s_%d", m_binSavePrefix, j),
                                      Form("Bin Projection %d", j), 800, 600);
        projY->SetTitle(Form("%s Bin: %.1e-%.1e", m_xLabel, 
                     h_RelRes_binned->GetXaxis()->GetBinLowEdge(j), 
                     h_RelRes_binned->GetXaxis()->GetBinUpEdge(j)));
        
        projY->SetStats(0);
        projY->Draw();
        
        double ymax = 1.1 * projY->GetMaximum();
        if (gaus) {
            ymax = 1.1 * std::fmax(projY->GetMaximum(), gaus->GetMaximum());
            gaus->SetLineColor(kRed);
            gaus->Draw("same");
        }
        projY->GetYaxis()->SetRangeUser(0, ymax);
        
        // Create stats box
        TPaveText* statsBox = new TPaveText(0.6, 0.6, 0.9, 0.9, "NDC");
        statsBox->SetFillColor(kWhite);
        statsBox->SetBorderSize(1);
        statsBox->SetTextAlign(12);
        statsBox->SetTextSize(0.035);
        statsBox->SetTextFont(42);

        statsBox->AddText(Form("Entries: %d", (int)projY->GetEntries()));
        statsBox->AddText(Form("#mu: %.5f", projY->GetMean()));
        statsBox->AddText(Form("RMS: %.5f", projY->GetRMS()));

        if (gaus) {
            statsBox->AddText(Form("#chi^{2}/ndf: %.1f/%d", gaus->GetChisquare(), gaus->GetNDF()));
            statsBox->AddText(Form("Fit #mu: %.5f #pm %.5f", mean, gaus->GetParError(1)));
            statsBox->AddText(Form("Fit #sigma: %.4f #pm %.5f", sigma, gaus->GetParError(2)));
            statsBox->AddText(Form("Fit range: [%.3f, %.3f]", m_xMinFit, m_xMaxFit));
        } else {
            statsBox->AddText("(Fit skipped)");
        }
        statsBox->Draw();

        // Add ePIC simulation labels
        TLatex latex_proj;
        latex_proj.SetTextSize(0.035);
        latex_proj.SetNDC();
        latex_proj.SetTextColor(kBlack);
        latex_proj.DrawLatex(0.15, 0.85, simLabel.c_str());
        latex_proj.DrawLatex(0.15, 0.80, "#bf{Diff. DIS} 10x100 GeV");

        c_proj->Update();
        SaveCanvas(c_proj, Form("figs/%s_bin_%d.png", m_binSavePrefix, j));
        
        // Add to fit graph only if fit was performed
        if (!skipFit) {
            g->SetPoint(j - 1, _center, mean);
            g->SetPointError(j - 1, 0.0, sigma);
        }
        
        // Always add RMS points
        if (m_isLogX){
            g_RMS->SetPoint(j - 1, _center * 1.05, projY->GetMean());
        } else {
            g_RMS->SetPoint(j - 1, _center + 0.005, projY->GetMean());
        }
        g_RMS->SetPointError(j - 1, 0.0, projY->GetRMS());

        delete c_proj;
        delete statsBox;
        if (gaus) delete gaus;
        delete projY;
    }
    
    // Rest of the plotting code remains the same...
    g_RMS->SetTitle(m_title);
    g_RMS->SetMarkerStyle(20);
    g_RMS->SetMarkerColor(kBlack);
    g_RMS->SetLineColor(kBlack);
    g_RMS->SetLineWidth(2);

    // Use base class range if set, otherwise fall back to constructor parameter
    auto xRange = m_rangeX ? *m_rangeX : m_xAxisRange;
    if (xRange.first != -999. && xRange.second != -999.) {
        g_RMS->GetXaxis()->SetLimits(xRange.first, xRange.second);
    }
    g_RMS->Draw("AP");

    g->SetTitle(m_title);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue + 2);
    g->SetLineColor(kBlue + 2);
    g->SetLineWidth(2);
    g->Draw("PSAME");

    g_RMS->SetTitle(m_title);

    TLine* line = new TLine(0, 0, g_RMS->GetXaxis()->GetXmax(), 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    TLegend* legend = new TLegend(m_legendLB->first, m_legendLB->second, m_legendRT->first, m_legendRT->second);
    legend->AddEntry(g_RMS, "Histograms", "ep");
    legend->AddEntry(g, "Gaussian Fit", "ep");
    legend->Draw();

    // Add ePIC simulation labels
    TLatex latex;
    latex.SetTextSize(0.035);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.15, 0.85, simLabel.c_str());
    latex.DrawLatex(0.15, 0.80, "#bf{Diff. DIS} 10x100 GeV");

    c->Update();
    SaveCanvas(c, m_saveName);
    delete c;
    delete g;
    delete g_RMS;
    delete line;
}
