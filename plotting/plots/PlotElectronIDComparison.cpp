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

void PlotElectronIDComparison(TFile* inputFile, const std::string& outDir) {
    if (!inputFile) return;
    gSystem->mkdir(outDir.c_str(), kTRUE);

    auto Get1D = [&](const char* n) { return dynamic_cast<TH1*>(inputFile->Get(n)); };
    auto Get2D = [&](const char* n) { return dynamic_cast<TH2*>(inputFile->Get(n)); };

    auto Overlay1D = [&](const char* oldName, const char* eidName,
                         const char* xTitle, const std::string& fileName,
                         bool logY) {
        TH1* hOld = Get1D(oldName);
        TH1* hEid = Get1D(eidName);
        if (!hOld || !hEid) {
            Logger::warning(std::string("PlotElectronIDComparison: missing '") + oldName +
                            "' or '" + eidName + "'");
            return;
        }
        hOld = (TH1*)hOld->Clone(TString(oldName) + "_cmp");
        hEid = (TH1*)hEid->Clone(TString(eidName) + "_cmp");
        hOld->SetDirectory(nullptr);
        hEid->SetDirectory(nullptr);

        TCanvas c("c_eid_cmp_1d", "", 800, 600);
        if (logY) c.SetLogy();
        hOld->SetLineColor(kBlack); hOld->SetLineWidth(2);
        hEid->SetLineColor(kRed+1); hEid->SetLineWidth(2); hEid->SetLineStyle(2);
        const double yMax = 1.25 * std::max(hOld->GetMaximum(), hEid->GetMaximum());
        hOld->GetYaxis()->SetRangeUser(logY ? 0.5 : 0.0, yMax > 0 ? yMax : 1.0);
        hOld->GetXaxis()->SetTitle(xTitle);
        hOld->GetYaxis()->SetTitle("Events");
        hOld->SetTitle("");
        hOld->Draw("hist");
        hEid->Draw("hist same");

        TLegend leg(0.58, 0.72, 0.88, 0.88);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(hOld, Form("old (objIdx)  N=%.0f", hOld->GetEntries()), "l");
        leg.AddEntry(hEid, Form("ElectronID    N=%.0f", hEid->GetEntries()), "l");
        leg.Draw();

        c.SaveAs((outDir + "/" + fileName).c_str());
        delete hOld; delete hEid;
    };

    Overlay1D("Ep_e",        "Ep_e_eid",        "E'_{e} [GeV]",       "e_energy_overlay.png", false);
    Overlay1D("pT_e",        "pT_e_eid",        "p_{T}^{e} [GeV]",    "e_pT_overlay.png",     false);
    Overlay1D("phi_e",       "phi_e_eid",       "#phi_{e} [rad]",     "e_phi_overlay.png",    false);
    Overlay1D("EPz_reco_mc", "EPz_eid",         "#Sigma(E-p_{z}) [GeV]", "event_EmPz_overlay.png", false);

    auto Draw2DCorrelation = [&](const char* name, const char* xTitle,
                                 const char* yTitle, const std::string& fileName,
                                 bool logZ) {
        TH2* h = Get2D(name);
        if (!h) {
            Logger::warning(std::string("PlotElectronIDComparison: missing 2D '") + name + "'");
            return;
        }
        h = (TH2*)h->Clone(TString(name) + "_cmp");
        h->SetDirectory(nullptr);

        TCanvas c("c_eid_cmp_2d", "", 700, 650);
        c.SetRightMargin(0.15);
        if (logZ) c.SetLogz();
        h->GetXaxis()->SetTitle(xTitle);
        h->GetYaxis()->SetTitle(yTitle);
        h->SetTitle("");
        h->Draw("colz");

        const double lo = std::min(h->GetXaxis()->GetXmin(), h->GetYaxis()->GetXmin());
        const double hi = std::max(h->GetXaxis()->GetXmax(), h->GetYaxis()->GetXmax());
        TLine diag(lo, lo, hi, hi);
        diag.SetLineColor(kBlack); diag.SetLineStyle(2); diag.SetLineWidth(2);
        diag.Draw("same");

        c.SaveAs((outDir + "/" + fileName).c_str());
        delete h;
    };

    Draw2DCorrelation("pT_e_old_vs_eid",
                      "p_{T}^{e} old (objIdx) [GeV]",
                      "p_{T}^{e} ElectronID [GeV]",
                      "e_pT_old_vs_eid.png", true);
    Draw2DCorrelation("Ep_e_old_vs_eid",
                      "E'_{e} old (objIdx) [GeV]",
                      "E'_{e} ElectronID [GeV]",
                      "e_energy_old_vs_eid.png", true);

    auto DrawDiff = [&](const char* name, const char* xTitle,
                        const std::string& fileName) {
        TH1* h = Get1D(name);
        if (!h) return;
        h = (TH1*)h->Clone(TString(name) + "_cmp");
        h->SetDirectory(nullptr);
        TCanvas c("c_eid_cmp_diff", "", 800, 600);
        h->SetLineColor(kAzure+2); h->SetLineWidth(2);
        h->GetXaxis()->SetTitle(xTitle);
        h->GetYaxis()->SetTitle("Events");
        h->SetTitle("");
        h->Draw("hist");
        TLatex lat;
        lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.036);
        lat.DrawLatex(0.16, 0.86, Form("mean = %.3g", h->GetMean()));
        lat.DrawLatex(0.16, 0.81, Form("rms = %.3g",  h->GetRMS()));
        c.SaveAs((outDir + "/" + fileName).c_str());
        delete h;
    };

    DrawDiff("dpT_e_old_eid",  "p_{T}^{eid} - p_{T}^{old} [GeV]", "e_dpT_old_eid.png");
    DrawDiff("dphi_e_old_eid", "#Delta#phi(old, eid) [rad]",      "e_dphi_old_eid.png");

    if (TH1* hCat = Get1D("e_finder_category")) {
        hCat = (TH1*)hCat->Clone("e_finder_category_cmp");
        hCat->SetDirectory(nullptr);
        TCanvas c("c_eid_cmp_cat", "", 800, 600);
        hCat->SetFillColor(kOrange-3); hCat->SetLineColor(kBlack);
        hCat->SetBarWidth(0.8); hCat->SetBarOffset(0.1);
        hCat->GetYaxis()->SetTitle("Events");
        hCat->SetTitle("");
        hCat->Draw("bar");
        const double tot = hCat->Integral();
        TLatex lat; lat.SetTextFont(42); lat.SetTextSize(0.035); lat.SetTextAlign(22);
        if (tot > 0) {
            for (int b = 1; b <= hCat->GetNbinsX(); ++b) {
                const double n = hCat->GetBinContent(b);
                const double x = hCat->GetXaxis()->GetBinCenter(b);
                lat.DrawLatex(x, n + 0.02 * hCat->GetMaximum(),
                              Form("%.0f (%.1f%%)", n, 100.0 * n / tot));
            }
        }
        c.SaveAs((outDir + "/e_finder_category.png").c_str());
        delete hCat;
    }

    Logger::info("PlotElectronIDComparison: wrote comparison plots to " + outDir);
}
