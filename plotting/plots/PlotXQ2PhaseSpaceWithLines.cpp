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

void PlotXQ2PhaseSpaceWithLines(TFile* inputFile, double sGeV2) {
    if (!inputFile || sGeV2 <= 0.0) return;

    TH2D* h_reco = (TH2D*)inputFile->Get("xQ2_reco");
    TH2D* h_truth = (TH2D*)inputFile->Get("xQ2_truth");

    // Build from tree if not stored
    TH2D* h_tmp = nullptr;
    if (!h_reco && !h_truth) {
        const int nx = 120, nq = 120;
        std::vector<Double_t> xb = GetLogBins(1.0e-4, 1.0, nx);
        std::vector<Double_t> qb = GetLogBins(1.0, 1000.0, nq);
        TTree* tree = (TTree*)inputFile->Get("Q2_tree");
        if (!tree) return;
        h_tmp = new TH2D("xQ2_phase_tmp", "", xb.size()-1, xb.data(), qb.size()-1, qb.data());
        float xv, qv;
        if (tree->GetBranch("x_truth") && tree->GetBranch("Q2_truth")) {
            tree->SetBranchAddress("x_truth", &xv);
            tree->SetBranchAddress("Q2_truth", &qv);
            for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
                tree->GetEntry(i);
                if (std::isfinite(xv) && xv > 0 && std::isfinite(qv) && qv > 0)
                    h_tmp->Fill(xv, qv);
            }
        }
    }

    TH2* h = h_reco ? (TH2*)h_reco : (h_truth ? (TH2*)h_truth : (TH2*)h_tmp);
    if (!h) { delete h_tmp; return; }

    TCanvas* c = new TCanvas("c_phase_lines", "", 1200, 1000);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);
    c->SetLogx();
    c->SetLogy();

    TH2* hd = (TH2*)h->Clone("xQ2_phase_draw");
    hd->SetDirectory(nullptr);
    hd->SetTitle("");
    hd->GetXaxis()->SetTitle("x_{Bj}");
    hd->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
    hd->GetXaxis()->SetTitleOffset(1.3);
    hd->GetYaxis()->SetTitleOffset(1.3);
    hd->Draw("COLZ");

    const double xlo = 1.0e-4, xhi = 1.0, q2lo = 1.0, q2hi = 1000.0;

    // iso-y lines: Q2 = s * x * y
    auto makeIsoY = [&](const char* name, double yval) {
        TF1* f = new TF1(name, "[0]*x*[1]", xlo, xhi);
        f->SetParameters(sGeV2, yval);
        f->SetLineColor(kBlack);
        f->SetLineStyle(2);
        f->SetLineWidth(2);
        f->SetRange(xlo, xhi);
        f->SetMinimum(q2lo);
        f->SetMaximum(q2hi);
        return f;
    };

    TF1* f_ylo = makeIsoY("iso_y_lo", 0.01);
    TF1* f_yhi = makeIsoY("iso_y_hi", 0.95);
    f_ylo->Draw("same");
    f_yhi->Draw("same");

    // iso-W2 line: Q2 = W2 * x / (1 - x)  (neglecting m_p^2)
    TF1* f_w2 = new TF1("iso_w2", "[0]*x/(1.0-x)", xlo, 0.98);
    f_w2->SetParameter(0, 20.0);
    f_w2->SetLineColor(kBlue + 2);
    f_w2->SetLineStyle(3);
    f_w2->SetLineWidth(2);
    f_w2->SetMinimum(q2lo);
    f_w2->SetMaximum(q2hi);
    f_w2->Draw("same");

    // Labels
    TLatex* lat = new TLatex();
    lat->SetNDC(false);
    lat->SetTextSize(0.030);
    lat->SetTextColor(kBlack);
    lat->DrawLatex(2e-2, sGeV2 * 2e-2 * 0.01 * 1.5, "y = 0.01");
    lat->DrawLatex(3e-4, std::min(sGeV2 * 3e-4 * 0.95 * 1.8, q2hi * 0.9), "y = 0.95");
    lat->SetTextColor(kBlue + 2);
    lat->DrawLatex(2e-2, 20.0 * 2e-2 / (1.0 - 2e-2) * 1.5, "W^{2} = 20 GeV^{2}");

    DrawSimLabels(inputFile);
    c->Update();
    SaveCanvas(c, "figs/Q2_x_phase_space.pdf");

    delete lat; delete f_ylo; delete f_yhi; delete f_w2; delete hd; delete c; delete h_tmp;
}
