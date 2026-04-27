#include "plots/DDISPlots.hpp"

#include "Plotting.hpp"
#include "Utility.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <vector>

void PlotXQ2Density(TFile* inputFile) {
    if (!inputFile) return;
    const int n_x_bins = 120;
    const int n_q2_bins = 120;
    const double q2_min = 1.0;
    const double q2_max = 1000.0;
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, n_x_bins);
    std::vector<Double_t> q2_bins = GetLogBins(q2_min, q2_max, n_q2_bins);

    auto makeHist = [&](const char* name, const char* title) {
        return new TH2D(name,
                        title,
                        x_bins.size() - 1, x_bins.data(),
                        q2_bins.size() - 1, q2_bins.data());
    };

    TH2D* h_reco_file = (TH2D*)inputFile->Get("xQ2_reco");
    TH2D* h_truth_file = (TH2D*)inputFile->Get("xQ2_truth");
    TH2D* h_reco_tmp = nullptr;
    TH2D* h_truth_tmp = nullptr;

    // Build missing reco/truth histograms from Q2_tree only when needed.
    if (!h_reco_file || !h_truth_file) {
        TTree* tree = (TTree*)inputFile->Get("Q2_tree");
        if (tree) {
            const bool hasRecoBranches = tree->GetBranch("x_EM") && tree->GetBranch("Q2_EM");
            const bool hasTruthBranches = tree->GetBranch("x_truth") && tree->GetBranch("Q2_truth");

            if (!h_reco_file && hasRecoBranches) {
                h_reco_tmp = makeHist("xQ2_reco_tmp", "Reco Event Density;x_{Bj};Q^{2} [GeV^{2}]");
            }
            if (!h_truth_file && hasTruthBranches) {
                h_truth_tmp = makeHist("xQ2_truth_tmp", "Truth Event Density;x_{Bj};Q^{2} [GeV^{2}]");
            }

            if (h_reco_tmp || h_truth_tmp) {
                float x_em = -999.0f, q2_em = -999.0f;
                float x_truth = -999.0f, q2_truth = -999.0f;
                if (h_reco_tmp) {
                    tree->SetBranchAddress("x_EM", &x_em);
                    tree->SetBranchAddress("Q2_EM", &q2_em);
                }
                if (h_truth_tmp) {
                    tree->SetBranchAddress("x_truth", &x_truth);
                    tree->SetBranchAddress("Q2_truth", &q2_truth);
                }

                const Long64_t nEntries = tree->GetEntries();
                for (Long64_t i = 0; i < nEntries; ++i) {
                    tree->GetEntry(i);
                    if (h_reco_tmp &&
                        std::isfinite(x_em) && x_em > 0.0f && x_em < 1.0f &&
                        std::isfinite(q2_em) && q2_em > 0.0f) {
                        h_reco_tmp->Fill(x_em, q2_em);
                    }
                    if (h_truth_tmp &&
                        std::isfinite(x_truth) && x_truth > 0.0f && x_truth < 1.0f &&
                        std::isfinite(q2_truth) && q2_truth > 0.0f) {
                        h_truth_tmp->Fill(x_truth, q2_truth);
                    }
                }
            }
        }
    }

    TH2* h_reco = h_reco_file ? static_cast<TH2*>(h_reco_file) : static_cast<TH2*>(h_reco_tmp);
    TH2* h_truth = h_truth_file ? static_cast<TH2*>(h_truth_file) : static_cast<TH2*>(h_truth_tmp);

    if (!h_reco && !h_truth) {
        Logger::warning("Could not find or build xQ2_reco/xQ2_truth; skipping x-Q2 density plots.");
        delete h_reco_tmp;
        delete h_truth_tmp;
        return;
    }

    auto drawDensity = [&](TH2* h_src, const char* canvasName, const char* saveName) {
        if (!h_src) return;
        TH2* h = (TH2*)h_src->Clone(Form("%s_draw_%s", h_src->GetName(), canvasName));
        h->SetDirectory(nullptr);
        TCanvas* c = new TCanvas(canvasName, "x-Q2 Density", 1200, 1000);
        gStyle->SetOptStat(0);
        c->SetRightMargin(0.15);
        c->SetLeftMargin(0.15);
        c->SetTopMargin(0.1);
        c->SetBottomMargin(0.1);
        c->SetLogx();
        c->SetLogy();

        h->SetTitle("");
        h->GetXaxis()->SetTitle("x_{Bj}");
        h->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h->GetXaxis()->SetTitleOffset(1.3);
        h->GetYaxis()->SetTitleOffset(1.3);
        h->Draw("COLZ");

        DrawSimLabels(inputFile);

        c->Update();
        SaveCanvas(c, saveName);

        delete h;
        delete c;
    };

    drawDensity(h_reco, "c_xq2_density_reco", "figs/inclusive/distributions/xbj_q2_density_reco.png");
    drawDensity(h_truth, "c_xq2_density_truth", "figs/inclusive/distributions/xbj_q2_density_truth.png");

    // Backward-compatible alias used by existing slides/docs.
    if (h_reco) {
        drawDensity(h_reco, "c_xq2_density_legacy", "figs/inclusive/distributions/xbj_q2_density.png");
    } else if (h_truth) {
        drawDensity(h_truth, "c_xq2_density_legacy", "figs/inclusive/distributions/xbj_q2_density.png");
    }

    delete h_reco_tmp;
    delete h_truth_tmp;
}
