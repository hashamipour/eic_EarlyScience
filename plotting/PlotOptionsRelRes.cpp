#include "Plotting.hpp"
#include <TF1.h>
#include <TLatex.h>
#include <TMath.h>
#include <RooArgSet.h>
#include <RooConstVar.h>
#include <RooDecay.h>
#include <RooNovosibirsk.h>
#include <RooRealVar.h>
#include <RooTruthModel.h>
#include <iostream>
#include <cmath>
#include <memory>
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

bool UseDSCB(const TString& fitFunction) {
    return fitFunction == "dscb" ||
           fitFunction == "doublecrystalball" ||
           fitFunction == "double_crystal_ball" ||
           fitFunction == "double_sided_crystalball" ||
           fitFunction == "double_sided_crystal_ball";
}

bool UseVoigt(const TString& fitFunction) {
    return fitFunction == "voigt" || fitFunction == "voigtian" || fitFunction == "voigt_profile";
}

bool UseBifurcatedGaussian(const TString& fitFunction) {
    return fitFunction == "bifurcatedgaussian" ||
           fitFunction == "bifurcated_gaussian" ||
           fitFunction == "bifurcated gaussian" ||
           fitFunction == "bifgaus" ||
           fitFunction == "bifurgaus";
}

bool UseNovosibirsk(const TString& fitFunction) {
    return fitFunction == "novosibirsk" ||
           fitFunction == "roo_novosibirsk" ||
           fitFunction == "roonovosibirsk" ||
           fitFunction == "novo";
}

bool UseGeneralizedGaussian(const TString& fitFunction) {
    return fitFunction == "generalizedgaussian" ||
           fitFunction == "generalized_gaussian" ||
           fitFunction == "generalized gaussian" ||
           fitFunction == "gen_gaus" ||
           fitFunction == "gengaus";
}

bool UseLaplaceDist(const TString& fitFunction) {
    return fitFunction == "laplace" ||
           fitFunction == "laplace_dist" ||
           fitFunction == "laplacedist" ||
           fitFunction == "laplace distribution";
}

bool UseBukinDist(const TString& fitFunction) {
    return fitFunction == "bukin" ||
           fitFunction == "bukin_dist" ||
           fitFunction == "bukindist" ||
           fitFunction == "bukin distribution";
}

bool UseRooJohnsonDist(const TString& fitFunction) {
    return fitFunction == "roojohnson" ||
           fitFunction == "roo_johnson" ||
           fitFunction == "johnson" ||
           fitFunction == "johnson_su" ||
           fitFunction == "johnsonsu";
}

bool UseRooDecayDist(const TString& fitFunction) {
    return fitFunction == "roodecay" ||
           fitFunction == "roo_decay" ||
           fitFunction == "bdecay" ||
           fitFunction == "b_decay" ||
           fitFunction == "roobdecay" ||
           fitFunction == "roo_bdecay";
}

bool HasVariableBinWidth(const TH1D* hist) {
    if (!hist || hist->GetNbinsX() < 2) return false;
    const double w0 = hist->GetXaxis()->GetBinWidth(1);
    for (int i = 2; i <= hist->GetNbinsX(); ++i) {
        const double wi = hist->GetXaxis()->GetBinWidth(i);
        const double scale = std::max(1.0, std::max(std::abs(w0), std::abs(wi)));
        if (std::abs(wi - w0) > 1e-9 * scale) return true;
    }
    return false;
}

const char* FitFunctionLabel(const TString& fitFunction) {
    if (UseDSCB(fitFunction)) return "Double-sided CrystalBall";
    if (UseCrystalBall(fitFunction)) return "CrystalBall";
    if (UseVoigt(fitFunction)) return "Voigt";
    if (UseBifurcatedGaussian(fitFunction)) return "Bifurcated Gaussian";
    if (UseNovosibirsk(fitFunction)) return "Novosibirsk";
    if (UseGeneralizedGaussian(fitFunction)) return "Generalized Gaussian";
    if (UseLaplaceDist(fitFunction)) return "Laplace Dist.";
    if (UseBukinDist(fitFunction)) return "Bukin";
    if (UseRooJohnsonDist(fitFunction)) return "RooJohnson";
    if (UseRooDecayDist(fitFunction)) return "RooDecay (Truth)";
    return "Gaussian";
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

double VoigtCore(double* x, double* p) {
    const double amplitude = p[0];
    const double mean = p[1];
    const double sigma = std::max(std::abs(p[2]), 1e-9);
    const double gamma = std::max(std::abs(p[3]), 1e-9);
    return amplitude * TMath::Voigt(x[0] - mean, sigma, gamma, 4);
}

double BifurcatedGaussianCore(double* x, double* p) {
    const double amplitude = p[0];
    const double mean = p[1];
    const double sigmaL = std::max(std::abs(p[2]), 1e-9);
    const double sigmaR = std::max(std::abs(p[3]), 1e-9);
    const double dx = x[0] - mean;
    const double sigma = (dx < 0.0) ? sigmaL : sigmaR;
    return amplitude * std::exp(-0.5 * (dx / sigma) * (dx / sigma));
}

double RooNovosibirskCore(double* x, double* p) {
    // Parameters:
    // p0=A, p1=peak, p2=width, p3=tail
    static RooRealVar obs("novo_obs", "novo_obs", 0.0, -10.0, 10.0);
    static RooRealVar peak("novo_peak", "novo_peak", 0.0, -5.0, 5.0);
    static RooRealVar width("novo_width", "novo_width", 0.01, 1e-6, 10.0);
    static RooRealVar tail("novo_tail", "novo_tail", 0.0, -10.0, 10.0);
    static RooNovosibirsk pdf("novo_pdf", "novo_pdf", obs, peak, width, tail);
    static RooArgSet vars(obs);

    const double amplitude = p[0];
    peak.setVal(p[1]);
    width.setVal(std::max(std::abs(p[2]), 1e-6));
    tail.setVal(std::clamp(p[3], -10.0, 10.0));
    obs.setVal(x[0]);

    return amplitude * pdf.getVal(&vars);
}

double BukinCore(double* x, double* p) {
    const double amplitude = p[0];
    const double Xp = p[1];
    const double sigp = std::max(std::abs(p[2]), 1e-9);
    const double xi = std::max(-0.999, std::min(0.999, p[3]));
    const double rho1 = std::max(-1.0, std::min(0.0, p[4]));
    const double rho2 = std::max(0.0, std::min(1.0, p[5]));

    const double consts = 2.0 * std::sqrt(2.0 * std::log(2.0));
    const double hp = sigp * consts;
    const double r3 = std::log(2.0);
    const double r4 = std::sqrt(xi * xi + 1.0);
    const double r1 = xi / r4;

    double r5 = 1.0;
    if (std::abs(xi) > std::exp(-6.0)) {
        const double denom = std::log(r4 + xi);
        if (std::abs(denom) > 1e-12) {
            r5 = xi / denom;
        }
    }

    const double x1 = Xp + (hp / 2.0) * (r1 - 1.0);
    const double x2 = Xp + (hp / 2.0) * (r1 + 1.0);

    double r2 = 0.0;
    if (x[0] < x1) {
        const double d = Xp - x1;
        if (std::abs(d) > 1e-12) {
            r2 = rho1 * (x[0] - x1) * (x[0] - x1) / (d * d)
               - r3
               + 4.0 * r3 * (x[0] - x1) / hp * r5 * r4 / ((r4 - xi) * (r4 - xi));
        }
    } else if (x[0] < x2) {
        if (std::abs(xi) > std::exp(-6.0)) {
            const double arg = 1.0 + 4.0 * xi * r4 * (x[0] - Xp) / hp;
            const double denom = std::log(1.0 + 2.0 * xi * (xi - r4));
            if (arg > 1e-12 && std::abs(denom) > 1e-12) {
                const double ratio = std::log(arg) / denom;
                r2 = -r3 * ratio * ratio;
            } else {
                r2 = -1e6;
            }
        } else {
            r2 = -4.0 * r3 * (x[0] - Xp) * (x[0] - Xp) / (hp * hp);
        }
    } else {
        const double d = Xp - x2;
        if (std::abs(d) > 1e-12) {
            r2 = rho2 * (x[0] - x2) * (x[0] - x2) / (d * d)
               - r3
               - 4.0 * r3 * (x[0] - x2) / hp * r5 * r4 / ((r4 + xi) * (r4 + xi));
        }
    }

    if (std::abs(r2) > 100.0) {
        return 0.0;
    }
    return amplitude * std::exp(r2);
}

double RooJohnsonCore(double* x, double* p) {
    // Johnson SU with an explicit amplitude scale:
    // f(x) = A * [delta/(lambda*sqrt(2pi)*sqrt(1+z^2))] * exp(-0.5*(gamma + delta*asinh(z))^2)
    // where z = (x-mu)/lambda.
    const double amplitude = p[0];
    const double mu = p[1];
    const double lambda = std::max(std::abs(p[2]), 1e-9);
    const double gamma = p[3];
    const double delta = std::max(std::abs(p[4]), 1e-9);
    const double z = (x[0] - mu) / lambda;
    const double norm = delta / (lambda * std::sqrt(2.0 * TMath::Pi()) * std::sqrt(1.0 + z * z));
    const double arg = gamma + delta * std::asinh(z);
    return amplitude * norm * std::exp(-0.5 * arg * arg);
}

double RooDecayTruthCore(double* x, double* p) {
    // Parameters:
    // p0=A, p1=mu (x-shift), p2=tau
    static RooRealVar t("rdec_truth_t", "rdec_truth_t", 0.0, -10.0, 10.0);
    static RooRealVar tau("rdec_truth_tau", "rdec_truth_tau", 0.03, 1e-6, 10.0);
    static RooTruthModel truthModel("rdec_truth_model", "rdec_truth_model", t);
    static RooDecay truthPdf("rdec_truth_pdf", "rdec_truth_pdf", t, tau, truthModel, RooDecay::DoubleSided);
    static RooArgSet obs(t);

    const double amplitude = p[0];
    const double mu = p[1];
    const double tauVal = std::max(std::abs(p[2]), 1e-6);

    t.setVal(x[0] - mu);
    tau.setVal(tauVal);

    return amplitude * truthPdf.getVal(&obs);
}

TF1* BuildRelResFitFunction(const TString& fitFunction,
                            const TString& functionName,
                            const double xMinFit,
                            const double xMaxFit,
                            const double amplitude,
                            const double mean,
                            const double sigma) {
    if (UseDSCB(fitFunction)) {
        TF1* fit = new TF1(functionName, DoubleSidedCrystalBallCore, xMinFit, xMaxFit, 7);
        fit->SetParNames("A", "Mean", "Sigma", "AlphaL", "nL", "AlphaR", "nR");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 1.5, 3.0, 1.5, 3.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 0.01, 25.0);
        fit->SetParLimits(4, 1.01, 100.0);
        fit->SetParLimits(5, 0.01, 25.0);
        fit->SetParLimits(6, 1.01, 100.0);
        return fit;
    }

    if (UseCrystalBall(fitFunction)) {
        TF1* fit = new TF1(functionName, CrystalBallCore, xMinFit, xMaxFit, 5);
        fit->SetParNames("A", "Mean", "Sigma", "Alpha", "n");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 1.5, 3.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 0.01, 25.0);
        fit->SetParLimits(4, 1.01, 100.0);
        return fit;
    }

    if (UseVoigt(fitFunction)) {
        TF1* fit = new TF1(functionName, VoigtCore, xMinFit, xMaxFit, 4);
        fit->SetParNames("A", "Mean", "Sigma", "Gamma");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), std::max(0.5 * sigma, 1e-6));
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 1e-6, 10.0);
        return fit;
    }

    if (UseBifurcatedGaussian(fitFunction)) {
        TF1* fit = new TF1(functionName, BifurcatedGaussianCore, xMinFit, xMaxFit, 4);
        fit->SetParNames("A", "Mean", "SigmaL", "SigmaR");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), std::max(sigma, 1e-6));
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 1e-6, 10.0);
        return fit;
    }

    if (UseNovosibirsk(fitFunction)) {
        TF1* fit = new TF1(functionName, RooNovosibirskCore, xMinFit, xMaxFit, 4);
        fit->SetParNames("A", "Peak", "Width", "Tail");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 0.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, -10.0, 10.0);
        return fit;
    }

    if (UseGeneralizedGaussian(fitFunction)) {
        TF1* fit = new TF1(
            functionName,
            "[0]*exp(-pow(abs((x-[1])/[2]),[3]))",
            xMinFit,
            xMaxFit
        );
        fit->SetParNames("A", "Mean", "Scale", "Beta");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 2.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, 0.1, 20.0);
        return fit;
    }

    if (UseLaplaceDist(fitFunction)) {
        TF1* fit = new TF1(
            functionName,
            "[0]*exp(-abs((x-[1])/[2]))",
            xMinFit,
            xMaxFit
        );
        fit->SetParNames("A", "Mean", "b");
        fit->SetParameters(amplitude, mean, std::max(sigma / std::sqrt(2.0), 1e-6));
        fit->SetParLimits(2, 1e-6, 10.0);
        return fit;
    }

    if (UseBukinDist(fitFunction)) {
        TF1* fit = new TF1(functionName, BukinCore, xMinFit, xMaxFit, 6);
        fit->SetParNames("A", "Xp", "Sigp", "Xi", "Rho1", "Rho2");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 0.0, -0.1, 0.1);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, -0.999, 0.999);
        fit->SetParLimits(4, -1.0, 0.0);
        fit->SetParLimits(5, 0.0, 1.0);
        return fit;
    }

    if (UseRooJohnsonDist(fitFunction)) {
        TF1* fit = new TF1(functionName, RooJohnsonCore, xMinFit, xMaxFit, 5);
        fit->SetParNames("A", "mu", "lambda", "gamma", "delta");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6), 0.0, 1.0);
        fit->SetParLimits(2, 1e-6, 10.0);
        fit->SetParLimits(3, -20.0, 20.0);
        fit->SetParLimits(4, 1e-3, 50.0);
        return fit;
    }

    if (UseRooDecayDist(fitFunction)) {
        TF1* fit = new TF1(functionName, RooDecayTruthCore, xMinFit, xMaxFit, 3);
        fit->SetParNames("A", "mu", "tau");
        fit->SetParameters(amplitude, mean, std::max(sigma, 1e-6));
        fit->SetParLimits(2, 1e-6, 10.0);
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
            
            hist->Fit(currentGaus, "RIQN", "", x_min_fit, x_max_fit);
            
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
        Logger::error("No valid fits found! (All fits had chi2/ndf > 10 or y-difference > 10%)");
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
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TH1D* sourceHist = (TH1D*)inputFile->Get(m_histName);
    if (!sourceHist) {
        std::cerr << "Error: Histogram " << m_histName << " not found." << std::endl;
        delete c;
        return;
    }
    std::unique_ptr<TH1D> histHolder(static_cast<TH1D*>(sourceHist->Clone(Form("%s_relres_work", m_histName.Data()))));
    TH1D* hist = histHolder.get();
    hist->SetDirectory(nullptr);

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(m_xLabel);
    hist->GetYaxis()->SetTitle(m_yLabel);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->SetStats(0);

    // Apply custom X range if set via base class
    if (m_rangeX) {
        hist->GetXaxis()->SetRangeUser(m_rangeX->first, m_rangeX->second);
    }

    hist->SetLineColor(kGreen + 3);
    hist->SetMarkerColor(kGreen + 3);
    hist->SetMarkerStyle(20);

    const double histMeanRaw = hist->GetMean();
    const double histRMSRaw = hist->GetRMS();
    const bool hasVariableBinWidth = HasVariableBinWidth(hist);
    if (hasVariableBinWidth) {
        // Width-normalize variable-bin histograms, equivalent to h->Scale(1.0, "width").
        hist->Scale(1.0, "width");
        if (std::string(m_yLabel) == "Counts") {
            hist->GetYaxis()->SetTitle("Counts / bin width");
        }
    }

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
    
    const char* fitOptions = hasVariableBinWidth ? "RQ+" : "RILQ+";
    hist->Fit(fitFunction, fitOptions);

    const bool useRooDecayTruth = UseRooDecayDist(m_fitFunction);
    fitFunction->SetLineColor(useRooDecayTruth ? kBlack : kRed);
    fitFunction->SetLineStyle(useRooDecayTruth ? 2 : 1);
    fitFunction->SetLineWidth(2);
    hist->GetYaxis()->SetRangeUser(0, fitFunction->GetMaximum(m_xMinFit, m_xMaxFit) * 1.1);
    fitFunction->Draw("same");

    double fitMean = fitFunction->Mean(m_xMinFit, m_xMaxFit);
    const double fitVar = fitFunction->Variance(m_xMinFit, m_xMaxFit);
    double fitSigma = (fitVar > 0.0) ? std::sqrt(fitVar) : -1.0;
    if (!std::isfinite(fitMean)) {
        fitMean = fitFunction->GetParameter(1);
    }
    if (!std::isfinite(fitSigma) || fitSigma < 0.0) {
        fitSigma = std::abs(fitFunction->GetParameter(2));
    }
    const double histMean = histMeanRaw;
    const double histRMS = histRMSRaw;
    const double chi2 = fitFunction->GetChisquare();
    const int ndf = fitFunction->GetNDF();

    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.025);
    latex->SetTextColor(kBlack);
    latex->SetTextAlign(13);
    latex->DrawLatex(
        0.12, 0.88,
        Form("Fit: %s | #chi^{2}/ndf=%.2f/%d | #mu_{fit}=%.5f | #mu_{hist}=%.5f | #sigma_{fit}=%.5f | RMS_{hist}=%.5f",
             FitFunctionLabel(m_fitFunction), chi2, ndf, fitMean, histMean, fitSigma, histRMS)
    );

    DrawSimLabels(inputFile);

    SaveCanvas(c, m_saveName);
    delete c;
    delete latex;
    delete fitFunction;
}
