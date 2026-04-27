#include "plots/DDISPlots.hpp"

#include "Plotting.hpp"
#include "Utility.hpp"
#include "YAMLBinning.hpp"
#include "PlotDrawing.hpp"
#include "GridDrawing.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

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
                             bool logY) {
    if (!inputFile) return;
    TGraph* g = (TGraph*)inputFile->Get(graphName);
    if (!g) {
        Logger::warning(std::string("Graph ") + graphName + " not found; skipping density plot.");
        return;
    }
    std::vector<Double_t> xbins = logX ? GetLogBins(xmin, xmax, nBins)
                                       : BuildLinEdges(xmin, xmax, nBins);
    std::vector<Double_t> ybins = logY ? GetLogBins(xmin, xmax, nBins)
                                       : BuildLinEdges(xmin, xmax, nBins);

    TH2D* h = new TH2D(Form("h_%s_density", graphName),
                       title,
                       xbins.size() - 1, xbins.data(),
                       ybins.size() - 1, ybins.data());

    const int n = g->GetN();
    double x = 0.0, y = 0.0;
    for (int i = 0; i < n; ++i) {
        g->GetPoint(i, x, y);
        if (!std::isfinite(x) || !std::isfinite(y)) continue;
        if (x <= 0.0 || y <= 0.0) continue;
        h->Fill(x, y);
    }

    TCanvas* c = new TCanvas(Form("c_%s_density", graphName), title, 1200, 1000);
    gStyle->SetOptStat(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.1);
    if (logX) c->SetLogx();
    if (logY) c->SetLogy();

    h->SetTitle("");
    h->GetXaxis()->SetTitle(xLabel);
    h->GetYaxis()->SetTitle(yLabel);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->Draw("COLZ");

    TLine* diag = new TLine(xmin, xmin, xmax, xmax);
    diag->SetLineColor(kBlack);
    diag->SetLineStyle(2);
    diag->SetLineWidth(2);
    diag->Draw("same");

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete diag;
    delete c;
    delete h;
}
