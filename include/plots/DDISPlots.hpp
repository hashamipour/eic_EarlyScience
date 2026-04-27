#pragma once

// Forward declarations for DDIS plotter entry points that have been split
// out of analysis/DDIS_Plot_Final.cpp into plotting/plots/*.cpp.

#include <string>
#include <vector>

class TFile;
class TH1;
class TH1D;
class TH2D;
class TH3D;

enum class TBinMetricKind {
    Acceptance,
    Purity,
    Efficiency
};

void PlotXQ2Density(TFile* inputFile);

void PlotMX2Comparison(TFile* inputFile,
                       const std::string& outpath,
                       bool logY);

void PlotEPzWithCuts(TFile* inputFile);

void PlotDensityFromHist(TFile* inputFile,
                         const char* histName,
                         const char* xLabel,
                         const char* yLabel,
                         const char* saveName,
                         const bool logX,
                         const bool logY);

void PlotDensityFromHistWithOverlay(TFile* inputFile,
                                    const char* histName,
                                    const char* xLabel,
                                    const char* yLabel,
                                    const char* saveName,
                                    const bool logX,
                                    const bool logY,
                                    const std::vector<double>& xbins,
                                    const std::vector<double>& ybins,
                                    const bool overlayLogX,
                                    const bool overlayLogY);

void PlotGraphDensity(TFile* inputFile,
                      const char* graphName,
                      const char* title,
                      const char* xLabel,
                      const char* yLabel,
                      const char* saveName,
                      double xmin,
                      double xmax,
                      int nBins,
                      bool logX,
                      bool logY);

void PlotPhaseSpaceSlices(TFile* inputFile, const std::string& yamlPath);

void PlotRecoSetUncorrected(TFile* inputFile,
                            const char* recoHistNameSetA,
                            const char* recoHistNameSetB,
                            const char* xLabel,
                            const char* title,
                            const char* saveName,
                            bool logX,
                            bool logY,
                            const char* setALabel,
                            const char* setBLabel,
                            bool normalizeToPDF = false);

void PlotRecoSetComparison(TFile* inputFile,
                           const char* effTruthHistNameSetA,
                           const char* effRecoHistNameSetA,
                           const char* recoHistNameSetA,
                           const char* recoHistNameSetB,
                           const char* xLabel,
                           const char* title,
                           const char* saveName,
                           const bool logX,
                           const bool logY,
                           const char* setALabel,
                           const char* setBLabel);

void PlotRecoSetComparisonBR(TFile* inputFile,
                             const char* effTruthHistNameSetA_B0,
                             const char* effRecoHistNameSetA_B0,
                             const char* recoHistNameSetA_B0,
                             const char* recoHistNameSetB_B0,
                             const char* effTruthHistNameSetA_RP,
                             const char* effRecoHistNameSetA_RP,
                             const char* recoHistNameSetA_RP,
                             const char* recoHistNameSetB_RP,
                             const char* xLabel,
                             const char* title,
                             const char* saveName,
                             const bool logX,
                             const bool logY,
                             const char* setALabel,
                             const char* setBLabel,
                             const bool includeSum,
                             const bool drawSetAAsTruth,
                             const bool drawSingleSetALine);

void PlotRecoSetComparisonBRWithSetBUncorrected(TFile* inputFile,
                                                const char* effTruthHistNameSetA_B0,
                                                const char* effRecoHistNameSetA_B0,
                                                const char* recoHistNameSetA_B0,
                                                const char* recoHistNameSetB_B0,
                                                const char* effTruthHistNameSetA_RP,
                                                const char* effRecoHistNameSetA_RP,
                                                const char* recoHistNameSetA_RP,
                                                const char* recoHistNameSetB_RP,
                                                const char* xLabel,
                                                const char* title,
                                                const char* saveName,
                                                const bool logX,
                                                const bool logY);

void PlotRecoSetComparisonBRSumOnly(TFile* inputFile,
                                    const char* effTruthHistNameSetA_B0,
                                    const char* effRecoHistNameSetA_B0,
                                    const char* recoHistNameSetA_B0,
                                    const char* recoHistNameSetB_B0,
                                    const char* effTruthHistNameSetA_RP,
                                    const char* effRecoHistNameSetA_RP,
                                    const char* recoHistNameSetA_RP,
                                    const char* recoHistNameSetB_RP,
                                    const char* xLabel,
                                    const char* title,
                                    const char* saveName,
                                    const bool logX,
                                    const bool logY,
                                    const char* setALabel,
                                    const char* setBLabel);

void PlotRelResVsK(TFile* inputFile,
                   const std::vector<std::string>& histNames,
                   const std::vector<std::string>& labels,
                   const std::string& title,
                   const std::string& outpath);

void PlotTBinMetric(TFile* inputFile,
                    const TBinMetricKind kind,
                    const char* title,
                    const char* yLabel,
                    const char* saveName);

void PlotResolutionComparison(TFile* inputFile,
                              const std::vector<const char*>& histNames,
                              const std::vector<const char*>& labels,
                              const char* xTitle,
                              const char* yTitle,
                              double xLo, double xHi,
                              double yLo, double yHi,
                              bool logX,
                              const char* saveName,
                              bool disableFit = true);

void Plot3DResponseMatrix(TFile* inputFile);

void PlotXQ2PhaseSpaceWithLines(TFile* inputFile, double sGeV2);

void PlotElectronIDComparison(TFile* inputFile, const std::string& outDir);
