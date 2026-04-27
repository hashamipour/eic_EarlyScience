#pragma once

#include <vector>

#include "YAMLBinning.hpp"  // BinDef + GetGlobalBinFromBinDef

class TH2;
class TH2D;
class TH3D;

// Grid and slice overlays on 2D/3D histograms used by the DDIS plotter.
// Rectangles are drawn per (xpom, Q2) or (iX, iY) cell; counts / metric
// values are printed inside each rectangle.

void DrawBinningGridWithCounts(TH2* hist,
                               const std::vector<double>& xbins,
                               const std::vector<double>& ybins,
                               bool logX,
                               bool logY,
                               bool showBinIds = true,
                               int  minEventsForBorder = 20);

void DrawBinsForSlice(TH2D* hist,
                      const std::vector<BinDef>& bins,
                      double beta_lo,
                      double beta_hi,
                      bool logX,
                      bool logY,
                      bool showBinIds = true,
                      int  minEventsForBorder = 20);

void DrawBinsForSliceMetric(TH2D* hist,
                            const std::vector<BinDef>& bins,
                            double beta_lo,
                            double beta_hi,
                            bool logX,
                            bool logY,
                            const std::vector<double>& metric_values,
                            const std::vector<double>& q2_edges,
                            const std::vector<double>& xpom_edges,
                            const std::vector<double>& beta_edges,
                            const std::vector<double>* metric_uncertainties = nullptr,
                            double valueTextScale = 1.0,
                            bool showBinIds = true,
                            int  minEventsForBorder = 20);

void DrawSliceGrid(TH3D* h3,
                   int sliceAxis,
                   const std::vector<double>& edges,
                   const char* projOpt,
                   const char* xTitle,
                   const char* yTitle,
                   const char* sliceLabel,
                   const char* saveName,
                   bool logX,
                   bool logY,
                   const std::vector<double>* overlayXBins = nullptr,
                   const std::vector<double>* overlayYBins = nullptr,
                   bool overlayLogX = false,
                   bool overlayLogY = false,
                   const std::vector<BinDef>* overlayBins = nullptr,
                   bool showBinIds = true);

void DrawSliceGridWithMetric(TH3D* h3,
                             int sliceAxis,
                             const std::vector<double>& edges,
                             const char* projOpt,
                             const char* xTitle,
                             const char* yTitle,
                             const char* sliceLabel,
                             const char* saveName,
                             bool logX,
                             bool logY,
                             const std::vector<BinDef>& overlayBins,
                             bool overlayLogX,
                             bool overlayLogY,
                             const std::vector<double>& metric_values,
                             const std::vector<double>& q2_edges,
                             const std::vector<double>& xpom_edges,
                             const std::vector<double>& beta_edges,
                             double zMin,
                             double zMax,
                             const std::vector<double>* metric_uncertainties = nullptr,
                             double valueTextScale = 1.0,
                             bool showBinIds = true);
