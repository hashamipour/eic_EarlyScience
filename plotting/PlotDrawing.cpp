#include "PlotDrawing.hpp"

#include <TColor.h>
#include <TMath.h>
#include <TStyle.h>

#include <iomanip>
#include <sstream>

std::vector<double> BuildLogEdges(double minVal, double maxVal, int nSlices) {
    std::vector<double> edges;
    if (nSlices <= 0 || minVal <= 0.0 || maxVal <= 0.0 || maxVal <= minVal) {
        return edges;
    }
    edges.reserve(nSlices + 1);
    const double logMin = TMath::Log10(minVal);
    const double logMax = TMath::Log10(maxVal);
    const double step = (logMax - logMin) / nSlices;
    for (int i = 0; i <= nSlices; ++i) {
        edges.push_back(TMath::Power(10.0, logMin + step * i));
    }
    return edges;
}

std::vector<double> BuildLinEdges(double minVal, double maxVal, int nSlices) {
    std::vector<double> edges;
    if (nSlices <= 0 || maxVal <= minVal) {
        return edges;
    }
    edges.reserve(nSlices + 1);
    const double step = (maxVal - minVal) / nSlices;
    for (int i = 0; i <= nSlices; ++i) {
        edges.push_back(minVal + step * i);
    }
    return edges;
}

std::string FormatRange(double lo, double hi) {
    auto fmt = [](double v) {
        std::ostringstream oss;
        if (v >= 0.01 && v < 1000.0) {
            oss << std::fixed << std::setprecision(2) << v;
        } else {
            oss << std::scientific << std::setprecision(2) << v;
        }
        return oss.str();
    };
    std::ostringstream out;
    out << "[" << fmt(lo) << ", " << fmt(hi) << "]";
    return out.str();
}

void SetGreenYellowRedPalette() {
    const Int_t nRGBs = 4;
    const Int_t nCont = 256;
    Double_t stops[nRGBs] = {0.0, 0.5, 0.75, 1.0};
    Double_t red[nRGBs]   = {0.18, 0.98, 0.98, 0.90};
    Double_t green[nRGBs] = {0.68, 0.90, 0.60, 0.20};
    Double_t blue[nRGBs]  = {0.38, 0.20, 0.20, 0.20};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
    gStyle->SetNumberContours(nCont);
}
