#include "Plotting.hpp"
#include <TF1.h>
#include <TLine.h>
#include <TLatex.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

namespace {

TString NormalizeFitFunction(const char* fitFunction) {
    TString normalized = fitFunction ? fitFunction : "gaus";
    normalized.ToLower();
    return normalized;
}

bool UseCrystalBall(const TString& fitFunction) {
    return fitFunction == "crystalball" || fitFunction == "crystal_ball" || fitFunction == "cb";
}

bool UseDSCB(const TString& fitFunction) {
    return fitFunction == "dscb" ||
           fitFunction == "doublecrystalball" ||
           fitFunction == "double_crystal_ball" ||
           fitFunction == "double_sided_crystalball" ||
           fitFunction == "double_sided_crystal_ball";
}

double CrystalBallCore(double* x, double* p) {
    const double amplitude = p[0];
    const double mean = p[1];
    const double sigma = std::max(std::abs(p[2]), 1e-9);
    const double alpha = p[3];
    const double n = std::max(p[4], 1.01);

    double t = (x[0] - mean) / sigma;
    if (alpha < 0.) {
        t = -t;
    }

    const double absAlpha = std::max(std::abs(alpha), 1e-9);
    if (t > -absAlpha) {
        return amplitude * std::exp(-0.5 * t * t);
    }

    const double A = std::pow(n / absAlpha, n) * std::exp(-0.5 * absAlpha * absAlpha);
    const double B = (n / absAlpha) - absAlpha;
    return amplitude * A * std::pow(B - t, -n);
}

double DoubleSidedCrystalBallCore(double* x, double* p) {
    const double amplitude = p[0];
    const double mean = p[1];
    const double sigma = std::max(std::abs(p[2]), 1e-9);
    const double alphaL = std::max(std::abs(p[3]), 1e-9);
    const double nL = std::max(p[4], 1.01);
    const double alphaR = std::max(std::abs(p[5]), 1e-9);
    const double nR = std::max(p[6], 1.01);

    const double t = (x[0] - mean) / sigma;

    if (t < -alphaL) {
        const double AL = std::pow(nL / alphaL, nL) * std::exp(-0.5 * alphaL * alphaL);
        const double BL = (nL / alphaL) - alphaL;
        return amplitude * AL * std::pow(BL - t, -nL);
    }

    if (t > alphaR) {
        const double AR = std::pow(nR / alphaR, nR) * std::exp(-0.5 * alphaR * alphaR);
        const double BR = (nR / alphaR) - alphaR;
        return amplitude * AR * std::pow(BR + t, -nR);
    }

    return amplitude * std::exp(-0.5 * t * t);
}

TF1* BuildFitFunction(const TString& fitFunction,
                      const TString& functionName,
                      const double xMinFit,
                      const double xMaxFit,
                      TH1D* hist) {
    const double amplitude = hist->GetMaximum();
    const double mean = hist->GetMean();
    const double sigma = std::max(hist->GetRMS(), 1e-6);

    if (UseCrystalBall(fitFunction)) {
        TF1* fit = new TF1(functionName, CrystalBallCore, xMinFit, xMaxFit, 5);
        fit->SetParNames("A", "Mean", "Sigma", "Alpha", "n");
        fit->SetParameters(amplitude, mean, sigma, 1.5, 3.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 0.01, 25.0);
        fit->SetParLimits(4, 1.01, 100.0);
        return fit;
    }

    if (UseDSCB(fitFunction)) {
        TF1* fit = new TF1(functionName, DoubleSidedCrystalBallCore, xMinFit, xMaxFit, 7);
        fit->SetParNames("A", "Mean", "Sigma", "AlphaL", "nL", "AlphaR", "nR");
        fit->SetParameters(amplitude, mean, sigma, 1.5, 3.0, 1.5, 3.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 0.01, 25.0);
        fit->SetParLimits(4, 1.01, 100.0);
        fit->SetParLimits(5, 0.01, 25.0);
        fit->SetParLimits(6, 1.01, 100.0);
        return fit;
    }

    TF1* fit = new TF1(functionName, "gaus", xMinFit, xMaxFit);
    fit->SetParameters(amplitude, mean, sigma);
    return fit;
}

const char* FitLegendLabel(const TString& fitFunction) {
    if (UseDSCB(fitFunction)) return "Double-sided CrystalBall Fit";
    if (UseCrystalBall(fitFunction)) return "CrystalBall Fit";
    return "Gaussian Fit";
}

const char* FitFunctionLabel(const TString& fitFunction) {
    if (UseDSCB(fitFunction)) return "Double-sided CrystalBall";
    if (UseCrystalBall(fitFunction)) return "CrystalBall";
    return "Gaussian";
}

}  // namespace

PlotOptionsBinnedRelRes::PlotOptionsBinnedRelRes(const TString& histName,
                                                 const char* title,
                                                 const char* xLabel,
                                                 const char* yLabel,
                                                 const std::vector<std::pair<double, double>>& fitRanges,
                                                 const char* saveName,
                                                 const char* binSavePrefix,
                                                 const std::pair<double, double>& x_axis_range,
                                                 const bool isLogX,
                                                 const char* fitFunction
                                                )
    : m_histName(histName),
      m_title(title),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_fitFunction(NormalizeFitFunction(fitFunction)),
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
        Logger::error("No valid fits found! (All fits had chi2/ndf > 10 or y-difference > 10%)");
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
    const bool stitched = (m_histNameRP.Length() > 0 && m_histNameB0.Length() > 0 &&
                           std::isfinite(m_boundaryValue));
    TH2D* h_rp = stitched ? (TH2D*)inputFile->Get(m_histNameRP) : nullptr;
    TH2D* h_b0 = stitched ? (TH2D*)inputFile->Get(m_histNameB0) : nullptr;
    TH2D* h_RelRes_binned = stitched ? h_rp : (TH2D*)inputFile->Get(m_histName);
    if (!h_RelRes_binned || (stitched && !h_b0)) {
        std::cerr << "Error: 2D Histogram " << (stitched ? m_histNameRP : m_histName)
                  << " not found." << std::endl;
        return;
    }
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TCanvas* c = new TCanvas("c_binned", "", 1200, 800);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.18);

    if (m_isLogX){c->SetLogx();}

    int nbinsX = h_RelRes_binned->GetNbinsX();
    TGraphErrors* g = new TGraphErrors(nbinsX);
    TGraphErrors* g_RMS = new TGraphErrors(nbinsX);
    TGraphErrors* g_RMS_RP = stitched ? new TGraphErrors(nbinsX) : nullptr;
    TGraphErrors* g_RMS_B0 = stitched ? new TGraphErrors(nbinsX) : nullptr;
    int nRP = 0, nB0 = 0;

    for (int j = 1; j <= nbinsX; ++j) {
        TH2D* h_src = h_RelRes_binned;
        bool isRPBin = true;
        if (stitched) {
            const double xCenter = h_RelRes_binned->GetXaxis()->GetBinCenter(j);
            isRPBin = (xCenter < m_boundaryValue);
            h_src = isRPBin ? h_rp : h_b0;
        }
        TH1D* projY = h_src->ProjectionY("", j, j);
        
        if (projY->GetEntries() == 0) {
            delete projY;
            continue;
        }

        double _center = h_RelRes_binned->GetXaxis()->GetBinCenter(j);
        
        // Check if we should skip fitting for this bin
        bool skipFit = m_disableFit;
        if (!skipFit && !m_fitRanges.empty() && j - 1 < static_cast<int>(m_fitRanges.size())) {
            if (TMath::AreEqualRel(m_fitRanges[j - 1].first, 0., 1e-6) &&
                TMath::AreEqualRel(m_fitRanges[j - 1].second, 0., 1e-6)) {
                skipFit = true;
            }
        }
        
        TF1* fitFunc = nullptr;
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
            
            fitFunc = BuildFitFunction(
                m_fitFunction,
                Form("fit_bin_%d", j),
                m_xMinFit,
                m_xMaxFit,
                projY
            );
            projY->Fit(fitFunc, "RQ");

            // Use fit moments (same definition as displayed in the plot header).
            mean = fitFunc->Mean(m_xMinFit, m_xMaxFit);
            const double fitVar = fitFunc->Variance(m_xMinFit, m_xMaxFit);
            sigma = (fitVar > 0.0) ? std::sqrt(fitVar) : -1.0;
            if (!std::isfinite(mean)) {
                mean = fitFunc->GetParameter(1);
            }
            if (!std::isfinite(sigma) || sigma < 0.0) {
                sigma = std::abs(fitFunc->GetParameter(2));
            }
        }
        
        // Always create and save the projection plot
        TCanvas* c_proj = new TCanvas(Form("c_proj_%s_%d", m_binSavePrefix, j),
                                      Form("Bin Projection %d", j), 800, 600);
        c_proj->SetBottomMargin(0.18);
        projY->SetTitle("");
        
        projY->SetStats(0);
        projY->Draw();
        
        double ymax = 1.1 * projY->GetMaximum();
        if (fitFunc) {
            ymax = 1.1 * std::fmax(projY->GetMaximum(), fitFunc->GetMaximum());
            fitFunc->SetLineColor(kRed);
            fitFunc->Draw("same");
        }
        projY->GetYaxis()->SetRangeUser(0, ymax);

        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.025);
        latex->SetTextColor(kBlack);
        latex->SetTextAlign(13);

        if (fitFunc) {
            latex->DrawLatex(
                0.12, 0.88,
                Form("Fit: %s | #chi^{2}/ndf=%.2f/%d | #mu_{fit}=%.5f | #mu_{hist}=%.5f | #sigma_{fit}=%.5f | RMS_{hist}=%.5f",
                     FitFunctionLabel(m_fitFunction),
                     fitFunc->GetChisquare(),
                     fitFunc->GetNDF(),
                     mean,
                     projY->GetMean(),
                     sigma,
                     projY->GetRMS())
            );
        } else {
            latex->DrawLatex(
                0.12, 0.88,
                Form("Fit: (skipped) | #mu_{hist}=%.5f | RMS_{hist}=%.5f",
                     projY->GetMean(),
                     projY->GetRMS())
            );
        }

        DrawSimLabels(inputFile);

        c_proj->Update();
        SaveCanvas(c_proj, Form("%s_bin_%d.png", m_binSavePrefix, j));
        
        // Add to fit graph only if fit was performed
        if (!skipFit) {
            g->SetPoint(j - 1, _center, mean);
            g->SetPointError(j - 1, 0.0, sigma);
        }
        
        // Always add RMS points
        const double xDraw = m_isLogX ? _center * 1.05 : _center + 0.005;
        g_RMS->SetPoint(j - 1, xDraw, projY->GetMean());
        g_RMS->SetPointError(j - 1, 0.0, projY->GetRMS());
        if (stitched) {
            if (isRPBin) {
                g_RMS_RP->SetPoint(nRP, xDraw, projY->GetMean());
                g_RMS_RP->SetPointError(nRP, 0.0, projY->GetRMS());
                ++nRP;
            } else {
                g_RMS_B0->SetPoint(nB0, xDraw, projY->GetMean());
                g_RMS_B0->SetPointError(nB0, 0.0, projY->GetRMS());
                ++nB0;
            }
        }

        delete c_proj;
        delete latex;
        if (fitFunc) delete fitFunc;
        delete projY;
    }
    
    // Rest of the plotting code remains the same...
    g_RMS->SetTitle(m_title ? m_title : "");
    g_RMS->SetMarkerStyle(20);
    g_RMS->SetMarkerColor(stitched ? kGray + 2 : kBlack);
    g_RMS->SetLineColor(stitched ? kGray + 2 : kBlack);
    g_RMS->SetLineWidth(2);
    if (stitched) {
        g_RMS->SetMarkerSize(0.0);
        g_RMS->SetLineStyle(0);
    }

    // Use base class range if set, otherwise fall back to constructor parameter
    auto xRange = m_rangeX ? *m_rangeX : m_xAxisRange;
    if (xRange.first != -999. && xRange.second != -999.) {
        g_RMS->GetXaxis()->SetLimits(xRange.first, xRange.second);
    }
    g_RMS->Draw("AP");

    const TString titleText = m_title ? m_title : "";
    const bool titleProvidesAxes = titleText.CountChar(';') >= 2;
    if (!titleProvidesAxes && m_xLabel && std::strlen(m_xLabel) > 0) {
        g_RMS->GetXaxis()->SetTitle(m_xLabel);
    }
    if (m_yLabel && std::strlen(m_yLabel) > 0) {
        g_RMS->GetYaxis()->SetTitle(m_yLabel);
    }
    if (m_rangeY) {
        g_RMS->GetYaxis()->SetRangeUser(m_rangeY->first, m_rangeY->second);
    }
    g_RMS->GetXaxis()->SetTitleSize(0.045);
    g_RMS->GetYaxis()->SetTitleSize(0.045);
    g_RMS->GetXaxis()->SetLabelSize(0.040);
    g_RMS->GetXaxis()->SetTitleOffset(1.2);
    g_RMS->GetYaxis()->SetLabelSize(0.040);
    g_RMS->GetYaxis()->SetTitleOffset(1.35);

    g_RMS->SetTitle("");

    if (!m_disableFit) {
        g->SetTitle("");
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlue + 2);
        g->SetLineColor(kBlue + 2);
        g->SetLineWidth(2);
        g->Draw("PSAME");
    }

    TLine* line = new TLine(0, 0, g_RMS->GetXaxis()->GetXmax(), 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    TLine* splitLine = nullptr;
    if (stitched) {
        if (g_RMS_RP && g_RMS_RP->GetN() > 0) {
            g_RMS_RP->SetMarkerStyle(20);
            g_RMS_RP->SetMarkerColor(kBlue + 1);
            g_RMS_RP->SetLineColor(kBlue + 1);
            g_RMS_RP->SetLineWidth(2);
            g_RMS_RP->Draw("PSAME");
        }
        if (g_RMS_B0 && g_RMS_B0->GetN() > 0) {
            g_RMS_B0->SetMarkerStyle(21);
            g_RMS_B0->SetMarkerColor(kRed + 1);
            g_RMS_B0->SetLineColor(kRed + 1);
            g_RMS_B0->SetLineWidth(2);
            g_RMS_B0->Draw("PSAME");
        }
        double yLo = g_RMS->GetYaxis()->GetXmin();
        double yHi = g_RMS->GetYaxis()->GetXmax();
        if (m_rangeY) { yLo = m_rangeY->first; yHi = m_rangeY->second; }
        splitLine = new TLine(m_boundaryValue, yLo, m_boundaryValue, yHi);
        splitLine->SetLineColor(kBlack);
        splitLine->SetLineStyle(7);
        splitLine->SetLineWidth(2);
        splitLine->Draw("SAME");
    }

    if (m_legendLB && m_legendRT) {
        TLegend* legend = new TLegend(m_legendLB->first, m_legendLB->second,
                                      m_legendRT->first, m_legendRT->second);
        if (stitched) {
            if (g_RMS_RP && g_RMS_RP->GetN() > 0) legend->AddEntry(g_RMS_RP, m_lowLabel ? m_lowLabel : "RP", "ep");
            if (g_RMS_B0 && g_RMS_B0->GetN() > 0) legend->AddEntry(g_RMS_B0, m_highLabel ? m_highLabel : "B0", "ep");
        } else {
            legend->AddEntry(g_RMS, "RMS (#mu #pm #sigma_{RMS})", "ep");
            if (!m_disableFit) {
                legend->AddEntry(g, Form("Fitted (#mu #pm #sigma_{fit}) [%s]", FitLegendLabel(m_fitFunction)), "ep");
            }
        }
        legend->Draw();
    }

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, m_saveName);
    delete c;
    delete g;
    delete g_RMS;
    if (g_RMS_RP) delete g_RMS_RP;
    if (g_RMS_B0) delete g_RMS_B0;
    delete line;
    if (splitLine) delete splitLine;
}
