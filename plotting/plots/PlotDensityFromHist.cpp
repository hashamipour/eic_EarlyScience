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

void PlotDensityFromHist(TFile* inputFile,
                                const char* histName,
                                const char* xLabel,
                                const char* yLabel,
                                const char* saveName,
                                const bool logX,
                                const bool logY) {
    if (!inputFile) return;
    TH2* h = (TH2*)inputFile->Get(histName);
    if (!h) {
        if (strstr(histName, "_reco") != nullptr) {
            std::string fallback = histName;
            fallback.replace(fallback.find("_reco"), 5, "_truth");
            h = (TH2*)inputFile->Get(fallback.c_str());
            if (h) {
                Logger::warning(std::string("Histogram ") + histName + " not found; using " +
                                fallback + " instead.");
            }
        }
        if (!h) {
            Logger::warning(std::string("Histogram ") + histName + " not found; skipping density plot.");
            return;
        }
    }

    TCanvas* c = new TCanvas(Form("c_%s", histName), histName, 1200, 1000);
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

    DrawSimLabels(inputFile);

    c->Update();
    SaveCanvas(c, saveName);
    delete c;
}
