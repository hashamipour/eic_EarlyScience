#include "Plotting.hpp"
#include <TF1.h>
#include <TLatex.h>
#include <iostream>
#include <cmath>
#include "Utility.hpp"

namespace {

TString NormalizeFitFunction(const char* fitFunction) {
    TString normalized = fitFunction ? fitFunction : "gaus";
    normalized.ToLower();
    return normalized;
}

bool UseCrystalBall(const TString& fitFunction) {
    return fitFunction == "crystalball" || fitFunction == "crystal_ball" || fitFunction == "cb";
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

TF1* BuildRelResFitFunction(const TString& fitFunction,
                            const TString& functionName,
                            const double xMinFit,
                            const double xMaxFit,
                            const double amplitude,
                            const double mean,
                            const double sigma) {
    if (UseCrystalBall(fitFunction)) {
        TF1* fit = new TF1(functionName, CrystalBallCore, xMinFit, xMaxFit, 5);
        fit->SetParNames("A", "Mean", "Sigma", "Alpha", "n");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 1.5, 3.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 0.01, 25.0);
        fit->SetParLimits(4, 1.01, 100.0);
        return fit;
    }

    TF1* fit = new TF1(functionName, "gaus", xMinFit, xMaxFit);
    fit->SetParameters(amplitude, mean, sigma);
    return fit;
}

}  // namespace

PlotOptionsRelRes::PlotOptionsRelRes(const TString& histName,
                                     const char* xLabel,
                                     const char* yLabel,
                                     double xMinFit,
                                     double xMaxFit,
                                     const char* saveName,
                                     const char* fitFunction)
    : m_histName(histName),
      m_xLabel(xLabel),
      m_yLabel(yLabel),
      m_fitFunction(NormalizeFitFunction(fitFunction)),
      m_xMinFit(xMinFit),
      m_xMaxFit(xMaxFit),
      m_saveName(saveName),
      m_bestMean(0),
      m_bestSigma(0),
      m_bestAmplitude(0) {}

void PlotOptionsRelRes::SetFitRangeByBins(TH1D* hist) {
    int peakBin = hist->GetMaximumBin();
    double peakCenter = hist->GetBinCenter(peakBin);

    struct FitResult {
        int nBinsLeft;
        int nBinsRight;
        int totalBins;
        double chi2_ndf;
        double yValueDiff;
        double yValueDiffPercent;
        double y_left;
        double y_right;
        double amplitude;
        double mean;
        double sigma;
        double xMinFit;
        double xMaxFit;
    };

    std::vector<FitResult> validFits;

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "   nBins_L | nBins_R | Total | Chi2/NDF | Y_L | Y_R | |f(x_L) - f(x_R)|% " << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    for (int nBinsLeft = 3; nBinsLeft <= 15; ++nBinsLeft) {
        for (int nBinsRight = 3; nBinsRight <= 15; ++nBinsRight) {
            
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
                result.yValueDiff = yValueDiff;
                result.yValueDiffPercent = yValueDiffPercent;
                result.y_left = y_left;
                result.y_right = y_right;
                result.amplitude = currentGaus->GetParameter(0);
                result.mean = currentGaus->GetParameter(1);
                result.sigma = currentGaus->GetParameter(2);
                result.xMinFit = x_min_fit;
                result.xMaxFit = x_max_fit;
                
                validFits.push_back(result);
                
                std::cout << TString::Format("   %7d | %7d | %5d | %8.4f | %6.1f | %6.1f | %8.2f%%", 
                    nBinsLeft, nBinsRight, result.totalBins, chi2_ndf, y_left, y_right, yValueDiffPercent) << std::endl;
            }
            
            delete currentGaus;
        }
    }
    
    if (validFits.empty()) {
        std::cerr << "Error: No valid fits found! (All fits had chi2/ndf > 10 or y-difference > 10%)" << std::endl;
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
    m_bestAmplitude = bestFit.amplitude;
    m_bestMean = bestFit.mean;
    m_bestSigma = bestFit.sigma;
    
    double finalXSymmetryValue = TMath::Abs(m_xMinFit + m_xMaxFit - 2.0 * peakCenter);
    
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
    std::cout << "Y-values: f(x_L)=" << bestFit.y_left << ", f(x_R)=" << bestFit.y_right << std::endl;
    std::cout << "Y-symmetry: " << bestFit.yValueDiffPercent << "%" << std::endl;
    std::cout << "X-symmetry value (original): " << finalXSymmetryValue << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}

void PlotOptionsRelRes::Plot(TFile* inputFile) {
    TCanvas* c = new TCanvas("c_relres", "", 1200, 800);
    c->SetBottomMargin(0.2);
    gStyle->SetOptTitle(0);

    TH1D* hist = (TH1D*)inputFile->Get(m_histName);
    if (!hist) {
        std::cerr << "Error: Histogram " << m_histName << " not found." << std::endl;
        delete c;
        return;
    }

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(m_xLabel);
    hist->GetYaxis()->SetTitle(m_yLabel);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(0.9);

    // Apply custom X range if set via base class
    if (m_rangeX) {
        hist->GetXaxis()->SetRangeUser(m_rangeX->first, m_rangeX->second);
    }

    hist->SetLineColor(kGreen + 3);
    hist->SetMarkerColor(kGreen + 3);
    hist->SetMarkerStyle(20);
    hist->Draw("pe");

    bool skipFit =
        (TMath::AreEqualRel(m_xMinFit, 0., 1e-6) && TMath::AreEqualRel(m_xMaxFit, 0., 1e-6));
    if (skipFit) {
        Logger::info("Skipping fitting as per user request: [" + std::string(m_saveName) + "]");

        DrawSimLabels(inputFile);

        SaveCanvas(c, m_saveName);
        delete c;
        return;
    }
    
    bool autoFitRange =
        (TMath::AreEqualRel(m_xMinFit, -999., 1e-6) && TMath::AreEqualRel(m_xMaxFit, -999., 1e-6));
    if (autoFitRange) {
        SetFitRangeByBins(hist);
        Logger::debug(" Auto. Chosen fit range: [" + std::to_string(m_xMinFit) + ", " + std::to_string(m_xMaxFit) + "]");

        if (autoFitRange) {
            Logger::warning("Skipping fitting as no fit range found. If you want a fit please provide fit range: [" + std::string(m_saveName) + "]");

            DrawSimLabels(inputFile);

            SaveCanvas(c, m_saveName);
            delete c;
            return;
        }
    } else {
        Logger::info("Using user-defined fit range: [" + std::to_string(m_xMinFit) + ", " + std::to_string(m_xMaxFit) + "] for " + std::string(m_saveName));
    }

    m_bestMean = hist->GetBinCenter(hist->GetMaximumBin());
    m_bestSigma = hist->GetRMS();
    m_bestAmplitude = hist->GetMaximum();
    
    TF1* fitFunction = BuildRelResFitFunction(
        m_fitFunction,
        "relres_fit",
        m_xMinFit,
        m_xMaxFit,
        m_bestAmplitude,
        m_bestMean,
        m_bestSigma
    );
    
    hist->Fit(fitFunction, "RQ+");

    fitFunction->SetLineColor(kRed);
    fitFunction->SetLineWidth(2);
    hist->GetYaxis()->SetRangeUser(0, fitFunction->GetMaximum() * 1.1);
    fitFunction->Draw("same");

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.04);
    latex->SetTextColor(kRed);
    latex->DrawLatex(0.2, 0.8, Form("RMS = %.2f", hist->GetRMS()));

    DrawSimLabels(inputFile);

    SaveCanvas(c, m_saveName);
    delete c;
    delete latex;
    delete fitFunction;
}
