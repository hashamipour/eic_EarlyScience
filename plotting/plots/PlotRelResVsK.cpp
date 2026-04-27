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

void PlotRelResVsK(TFile* inputFile,
                          const std::vector<std::string>& histNames,
                          const std::vector<std::string>& labels,
                          const std::string& title,
                          const std::string& outpath) {
    if (!inputFile) return;
    std::vector<TH1*> hists;
    hists.reserve(histNames.size());
    for (const auto& name : histNames) {
        TH1* h = dynamic_cast<TH1*>(inputFile->Get(name.c_str()));
        if (h) hists.push_back(h);
    }
    if (hists.empty()) return;

    TCanvas c("c_relres_vs_k", "c_relres_vs_k", 3000, 600);
    gStyle->SetOptTitle(0);
    c.SetGridx();
    c.SetGridy();

    double ymax = 0.0;
    for (const auto* h : hists) {
        if (!h) continue;
        ymax = std::max(ymax, h->GetMaximum());
    }
    if (ymax <= 0.0) ymax = 0.1;

    const int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
    for (size_t i = 0; i < hists.size(); i++) {
        TH1* h = hists[i];
        if (!h) continue;
        h->SetTitle("");
        h->SetMarkerStyle(20 + static_cast<int>(i));
        h->SetMarkerSize(0.9);
        h->SetLineWidth(2);
        h->SetMarkerColor(colors[i % 4]);
        h->SetLineColor(colors[i % 4]);
        h->GetYaxis()->SetRangeUser(0.0, 1.2 * ymax);
        if (i == 0) {
            h->Draw("E1");
        } else {
            h->Draw("E1 SAME");
        }
    }

    if (labels.size() == hists.size()) {
        TLegend leg(0.7, 0.78, 0.92, 0.92);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        for (size_t i = 0; i < hists.size(); i++) {
            leg.AddEntry(hists[i], labels[i].c_str(), "lep");
        }
        leg.Draw();
    }

    DrawSimLabels(inputFile);
    SaveCanvas(&c, outpath.c_str());
}
