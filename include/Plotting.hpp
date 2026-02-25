#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <vector>
#include <TStyle.h>
#include <TColor.h>
#include <string>
#include <optional>
#include "Utility.hpp"  // Use your existing utility functions

// Base class for plot options
class PlotOptions {
    protected:
        // std::pair<double, double> m_legendLB;// = {0.75, 0.8}; // Legend lower-left  corner
        // std::pair<double, double> m_legendRT;// = {1.05, 1.0}; // Legend upper-right corner
        // Now they are optional pairs
        std::optional<std::pair<double, double>> m_legendLB;
        std::optional<std::pair<double, double>> m_legendRT;
        std::optional<std::pair<double, double>> m_rangeX;
        std::optional<std::pair<double, double>> m_rangeY;
    public:
        virtual ~PlotOptions();
        virtual void Plot(TFile* inputFile) = 0;
        void SetLegendPosition(double xlb, double ylb, double xrt, double yrt) {
            m_legendLB = {xlb, ylb};
            m_legendRT = {xrt, yrt};
            // Logger::debug("Legend position set to: (" + std::to_string(xlb) + ", " + std::to_string(ylb) + ") - (" + std::to_string(xrt) + ", " + std::to_string(yrt) + ")");
        }
        void SetRangeX(double xmin, double xmax) {
            m_rangeX = {xmin, xmax};
        }
        void SetRangeY(double ymin, double ymax) {
            m_rangeY = {ymin, ymax};
        }
};

// 1D histogram plot options
class PlotOptions1D : public PlotOptions {
private:
    std::vector<TString> m_histNames;
    std::vector<const char*> m_legendEntries;
    std::vector<const char*> m_drawOptions;
    const char* m_canvasTitle;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    bool m_isLogX;
    bool m_isLogY;
    bool m_normalizeToPDF;

public:
    PlotOptions1D(const std::vector<TString>& histNames,
                  const std::vector<const char*>& legendEntries,
                  const std::vector<const char*>& drawOptions,
                  const char* canvasTitle,
                  const char* xLabel,
                  const char* yLabel,
                  const char* saveName,
                  const bool  isLogX = false,
                  const bool  isLogY = false,
                  const bool  normalizeToPDF = false);
    
    void Plot(TFile* inputFile) override;
};

// Relative resolution plot options
class PlotOptionsRelRes : public PlotOptions {
private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    TString m_fitFunction;
    double m_xMinFit;
    double m_xMaxFit;
    const char* m_saveName;
    // optimization variables for fit
    double m_bestMean;
    double m_bestSigma;
    double m_bestAmplitude;

    void SetFitRangeByBins(TH1D* hist);

public:
    PlotOptionsRelRes(const TString& histName,
                      const char* xLabel,
                      const char* yLabel,
                      double xMinFit,
                      double xMaxFit,
                      const char* saveName,
                      const char* fitFunction = "gaus"
                    );
    
    void Plot(TFile* inputFile) override;
};
//////////////////////////////////////////////////////////////////////////////////
// Binned relative resolution plot options
class PlotOptionsBinnedRelRes : public PlotOptions {
private:
    TString m_histName;
    const char* m_title;
    const char* m_xLabel;
    const char* m_yLabel;
    TString m_fitFunction;
    std::vector<std::pair<double, double>> m_fitRanges;
    double m_xMinFit;
    double m_xMaxFit;
    const char* m_saveName;
    const char* m_binSavePrefix;
    std::pair<double, double> m_xAxisRange;
    const bool m_isLogX;

    void SetFitRangeByBins(TH1D* hist);

public:
    PlotOptionsBinnedRelRes(const TString& histName,
                            const char* title,
                            const char* xLabel,
                            const char* yLabel,
                            const std::vector<std::pair<double, double>>& fitRanges,
                            const char* saveName,
                            const char* binSavePrefix,
                            const std::pair<double, double>& x_axis_range = {-999., -999.},
                            const bool  isLogX = false,
                            const char* fitFunction = "gaus"
                        );
    
    void Plot(TFile* inputFile) override;
};

// Response matrix plot options
class PlotOptionsResponseMatrix : public PlotOptions {
private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    bool m_isLogX;
    bool m_isLogY;
    std::pair<double, double> m_xRange;
    std::pair<double, double> m_yRange;

public:
    PlotOptionsResponseMatrix(const TString& histName,
                              const char* xLabel,
                              const char* yLabel,
                              const char* saveName,
                              const bool isLogX = false,
                              const bool isLogY = false,
                              const std::pair<double, double>& xRange = {-999., -999.},
                              const std::pair<double, double>& yRange = {-999., -999.});
    
    void Plot(TFile* inputFile) override;
};

// Simple 2D histogram plot (density)
class PlotOptions2D : public PlotOptions {
private:
    TString m_histName;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    bool m_isLogX;
    bool m_isLogY;
    std::pair<double, double> m_xRange;
    std::pair<double, double> m_yRange;

public:
    PlotOptions2D(const TString& histName,
                  const char* xLabel,
                  const char* yLabel,
                  const char* saveName,
                  const bool isLogX = false,
                  const bool isLogY = false,
                  const std::pair<double, double>& xRange = {-999., -999.},
                  const std::pair<double, double>& yRange = {-999., -999.});

    void Plot(TFile* inputFile) override;
};

// NEW CLASS: Combined correlation plot options
class PlotOptionsCombinedCorrelation : public PlotOptions {
private:
    std::vector<TString> m_histNames;
    std::vector<const char*> m_legendEntries;
    std::vector<Color_t> m_colors;
    std::vector<Style_t> m_markerStyles;
    const char* m_canvasTitle;
    const char* m_xLabel;
    const char* m_yLabel;
    const char* m_saveName;
    std::pair<double, double> m_xRange;
    std::pair<double, double> m_yRange;
    bool m_isLogX;
    bool m_isLogY;

public:
    PlotOptionsCombinedCorrelation(const std::vector<TString>& histNames,
                                   const std::vector<const char*>& legendEntries,
                                   const std::vector<Color_t>& colors,
                                   const std::vector<Style_t>& markerStyles,
                                   const char* canvasTitle,
                                   const char* xLabel,
                                   const char* yLabel,
                                   const char* saveName,
                                   const std::pair<double, double>& xRange = {-999., -999.},
                                   const std::pair<double, double>& yRange = {-999., -999.},
                                   const bool isLogX = false,
                                   const bool isLogY = false);
    
    void Plot(TFile* inputFile) override;
};

void SaveCanvas(TCanvas* canvas, const char* filename);
std::string BuildSimLabel(TFile* inputFile);
void DrawSimLabels(TFile* inputFile);

#endif // PLOTTING_HPP
