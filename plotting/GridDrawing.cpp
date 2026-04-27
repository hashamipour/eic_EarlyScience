#include "GridDrawing.hpp"
#include "PlotDrawing.hpp"
#include "Plotting.hpp"
#include "Utility.hpp"

#include <TBox.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TPad.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>

void DrawBinningGridWithCounts(TH2* hist,
                               const std::vector<double>& xbins,
                               const std::vector<double>& ybins,
                               bool logX,
                               bool logY,
                               bool showBinIds,
                               int  minEventsForBorder) {
    if (!hist || xbins.size() < 2 || ybins.size() < 2) {
        if (hist && (xbins.size() < 2 || ybins.size() < 2)) {
            Logger::warning("bin overlay skipped (empty bin edges).");
        }
        return;
    }

    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double ymin = hist->GetYaxis()->GetXmin();
    const double ymax = hist->GetYaxis()->GetXmax();

    const int lineColor = kRed;
    const int lineWidth = 3;

    TLatex countText;
    countText.SetTextColor(kBlack);
    countText.SetTextSize(0.022);
    countText.SetTextAlign(22);
    countText.SetTextFont(42);

    TLatex idText;
    idText.SetTextColor(kBlue + 2);
    idText.SetTextSize(0.018);
    idText.SetTextAlign(22);
    idText.SetTextFont(42);

    const double tiny = 1e-12;
    const int nx = static_cast<int>(xbins.size()) - 1;
    for (size_t ix = 0; ix + 1 < xbins.size(); ++ix) {
        const double xlow = xbins[ix];
        const double xhigh = xbins[ix + 1];
        if (xhigh <= xmin || xlow >= xmax) continue;
        const double xlo = std::max(xlow, xmin);
        const double xhi = std::min(xhigh, xmax);
        const double xcenter = logX ? std::sqrt(xlo * xhi) : 0.5 * (xlo + xhi);
        const int binxlow = hist->GetXaxis()->FindBin(xlo + tiny);
        const int binxhigh = hist->GetXaxis()->FindBin(xhi - tiny);

        for (size_t iy = 0; iy + 1 < ybins.size(); ++iy) {
            const double ylow = ybins[iy];
            const double yhigh = ybins[iy + 1];
            if (yhigh <= ymin || ylow >= ymax) continue;
            const double ylo = std::max(ylow, ymin);
            const double yhi = std::min(yhigh, ymax);
            const double ycenter = logY ? std::sqrt(ylo * yhi) : 0.5 * (ylo + yhi);
            const int binylow = hist->GetYaxis()->FindBin(ylo + tiny);
            const int binyhigh = hist->GetYaxis()->FindBin(yhi - tiny);
            const int bin_id = static_cast<int>(ix) + static_cast<int>(iy) * nx;
            const double y_id = logY ? (ycenter * std::pow(yhi / ylo, 0.18))
                                     : (ycenter + 0.16 * (yhi - ylo));
            const double y_count = logY ? (ycenter / std::pow(yhi / ylo, 0.10))
                                        : (ycenter - 0.08 * (yhi - ylo));

            const double count = hist->Integral(binxlow, binxhigh, binylow, binyhigh);

            if (count >= minEventsForBorder) {
                TBox* box = new TBox(xlo, ylo, xhi, yhi);
                box->SetFillStyle(0);
                box->SetLineColor(lineColor);
                box->SetLineWidth(lineWidth);
                box->Draw("same");
            }
            if (count > 0.0) {
                countText.DrawLatex(xcenter, y_count, Form("%d", (int)std::lround(count)));
            }
            if (showBinIds) {
                idText.DrawLatex(xcenter, y_id, Form("%d", bin_id));
            }
        }
    }
}

void DrawBinsForSlice(TH2D* hist,
                      const std::vector<BinDef>& bins,
                      double beta_lo,
                      double beta_hi,
                      bool logX,
                      bool logY,
                      bool showBinIds,
                      int  minEventsForBorder) {
    if (!hist || bins.empty()) return;
    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double ymin = hist->GetYaxis()->GetXmin();
    const double ymax = hist->GetYaxis()->GetXmax();
    const double eps = 1e-12;

    TLatex countText;
    countText.SetTextColor(kBlack);
    countText.SetTextSize(0.050);
    countText.SetTextAlign(22);
    countText.SetTextFont(62);

    TLatex idText;
    idText.SetTextColor(kBlue + 2);
    idText.SetTextSize(0.030);
    idText.SetTextAlign(22);
    idText.SetTextFont(42);

    for (const auto& b : bins) {
        if (b.beta_max <= beta_lo + eps || b.beta_min >= beta_hi - eps) {
            continue;
        }
        double xlow = b.xpom_min;
        double xhigh = b.xpom_max;
        double ylow = b.Q2_min;
        double yhigh = b.Q2_max;
        if (xhigh <= xmin || xlow >= xmax || yhigh <= ymin || ylow >= ymax) {
            continue;
        }
        xlow = std::max(xlow, xmin);
        xhigh = std::min(xhigh, xmax);
        ylow = std::max(ylow, ymin);
        yhigh = std::min(yhigh, ymax);
        if (xlow <= 0.0 || ylow <= 0.0) continue;

        const double xcenter = logX ? std::sqrt(xlow * xhigh) : 0.5 * (xlow + xhigh);
        const double ycenter = logY ? std::sqrt(ylow * yhigh) : 0.5 * (ylow + yhigh);
        const int binxlow = hist->GetXaxis()->FindBin(xlow + eps);
        const int binxhigh = hist->GetXaxis()->FindBin(xhigh - eps);
        const int binylow = hist->GetYaxis()->FindBin(ylow + eps);
        const int binyhigh = hist->GetYaxis()->FindBin(yhigh - eps);
        const double count = hist->Integral(binxlow, binxhigh, binylow, binyhigh);

        if (count >= minEventsForBorder) {
            TBox* box = new TBox(xlow, ylow, xhigh, yhigh);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->Draw("same");
        }
        if (count > 0.0) {
            countText.DrawLatex(xcenter, ycenter, Form("%d", (int)std::lround(count)));
        }
        if (showBinIds && b.bin_id >= 0) {
            const double y_id = logY ? (ycenter * 1.25) : (ycenter + 0.06 * (yhigh - ylow));
            idText.DrawLatex(xcenter, y_id, Form("%d", b.bin_id));
        }
    }
}

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
                            const std::vector<double>* metric_uncertainties,
                            double valueTextScale,
                            bool showBinIds,
                            int  minEventsForBorder) {
    if (!hist || bins.empty()) return;
    const double xmin = hist->GetXaxis()->GetXmin();
    const double xmax = hist->GetXaxis()->GetXmax();
    const double ymin = hist->GetYaxis()->GetXmin();
    const double ymax = hist->GetYaxis()->GetXmax();
    const double eps = 1e-12;

    TLatex metricText;
    metricText.SetTextColor(kBlack);
    metricText.SetTextSize(0.050 * valueTextScale);
    metricText.SetTextAlign(22);
    metricText.SetTextFont(62);

    TLatex idText;
    idText.SetTextColor(kBlue + 2);
    idText.SetTextSize(0.030);
    idText.SetTextAlign(22);
    idText.SetTextFont(42);

    for (const auto& b : bins) {
        if (b.beta_max <= beta_lo + eps || b.beta_min >= beta_hi - eps) {
            continue;
        }
        double xlow = b.xpom_min;
        double xhigh = b.xpom_max;
        double ylow = b.Q2_min;
        double yhigh = b.Q2_max;
        if (xhigh <= xmin || xlow >= xmax || yhigh <= ymin || ylow >= ymax) {
            continue;
        }
        xlow = std::max(xlow, xmin);
        xhigh = std::min(xhigh, xmax);
        ylow = std::max(ylow, ymin);
        yhigh = std::min(yhigh, ymax);
        if (xlow <= 0.0 || ylow <= 0.0) continue;

        const double xcenter = logX ? std::sqrt(xlow * xhigh) : 0.5 * (xlow + xhigh);
        const double ycenter = logY ? std::sqrt(ylow * yhigh) : 0.5 * (ylow + yhigh);

        const int binxlow  = hist->GetXaxis()->FindBin(xlow + eps);
        const int binxhigh = hist->GetXaxis()->FindBin(xhigh - eps);
        const int binylow  = hist->GetYaxis()->FindBin(ylow + eps);
        const int binyhigh = hist->GetYaxis()->FindBin(yhigh - eps);
        const double count = hist->Integral(binxlow, binxhigh, binylow, binyhigh);

        if (count >= minEventsForBorder) {
            TBox* box = new TBox(xlow, ylow, xhigh, yhigh);
            box->SetFillStyle(0);
            box->SetLineColor(kRed);
            box->SetLineWidth(2);
            box->Draw("same");
        }

        const int k = GetGlobalBinFromBinDef(b, q2_edges, xpom_edges, beta_edges);
        if (k >= 0 && k < static_cast<int>(metric_values.size()) &&
            std::isfinite(metric_values[k]) && metric_values[k] >= 0.0) {
            const bool hasUnc = (metric_uncertainties &&
                                 k < static_cast<int>(metric_uncertainties->size()) &&
                                 std::isfinite((*metric_uncertainties)[k]) &&
                                 (*metric_uncertainties)[k] >= 0.0);
            if (hasUnc) {
                metricText.DrawLatex(xcenter, ycenter,
                                     Form("%.2f#pm%.2f", metric_values[k], (*metric_uncertainties)[k]));
            } else {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << metric_values[k];
                metricText.DrawLatex(xcenter, ycenter, oss.str().c_str());
            }
        }

        if (showBinIds && b.bin_id >= 0) {
            const double y_id = logY ? (ycenter * 1.25) : (ycenter + 0.06 * (yhigh - ylow));
            idText.DrawLatex(xcenter, y_id, Form("%d", b.bin_id));
        }
    }
}

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
                   const std::vector<double>* overlayXBins,
                   const std::vector<double>* overlayYBins,
                   bool overlayLogX,
                   bool overlayLogY,
                   const std::vector<BinDef>* overlayBins,
                   bool showBinIds) {
    if (!h3 || edges.size() < 2) return;
    const int nSlices = static_cast<int>(edges.size()) - 1;
    const int nCols = 2;
    const int nRows = (nSlices + nCols - 1) / nCols;

    TCanvas* c = new TCanvas(Form("c_phase_slices_%d", sliceAxis), saveName, 2400, 2000);
    c->Divide(nCols, nRows, 0.002, 0.002);

    TAxis* axis = nullptr;
    if (sliceAxis == 1) axis = h3->GetXaxis();
    if (sliceAxis == 2) axis = h3->GetYaxis();
    if (sliceAxis == 3) axis = h3->GetZaxis();
    if (!axis) {
        delete c;
        return;
    }

    std::vector<TH2D*> projections;
    projections.reserve(nSlices);

    for (int i = 0; i < nSlices; ++i) {
        const double lo = edges[i];
        const double hi = edges[i + 1];
        const double tiny = 1e-12;
        const int binLo = axis->FindBin(lo + tiny);
        const int binHi = axis->FindBin(hi - tiny);
        axis->SetRange(binLo, binHi);

        TH2D* h2raw = (TH2D*)h3->Project3D(projOpt);
        if (!h2raw) continue;
        TH2D* h2 = (TH2D*)h2raw->Clone(Form("slice_%s_%d", projOpt, i));
        h2->SetDirectory(nullptr);
        delete h2raw;
        h2->SetTitle("");
        h2->GetXaxis()->SetTitle(xTitle);
        h2->GetYaxis()->SetTitle(yTitle);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.2);
        if (strcmp(xTitle, "#beta") == 0) {
            h2->GetXaxis()->SetRangeUser(0.0, 1.0);
        }
        if (strcmp(yTitle, "#beta") == 0) {
            h2->GetYaxis()->SetRangeUser(0.0, 1.0);
        }

        TPad* pad = (TPad*)c->cd(i + 1);
        pad->SetRightMargin(0.14);
        pad->SetLeftMargin(0.14);
        pad->SetTopMargin(0.12);
        pad->SetBottomMargin(0.12);
        if (logX) pad->SetLogx();
        if (logY) pad->SetLogy();
        h2->Draw("COLZ");

        if (overlayBins) {
            DrawBinsForSlice(h2, *overlayBins, lo, hi, overlayLogX, overlayLogY, showBinIds);
        } else if (overlayXBins && overlayYBins) {
            DrawBinningGridWithCounts(h2, *overlayXBins, *overlayYBins, overlayLogX, overlayLogY, showBinIds);
        }

        TLatex latex;
        latex.SetTextSize(0.06);
        latex.SetNDC();
        latex.SetTextColor(kBlack);
        const std::string label = std::string(sliceLabel) + " " + FormatRange(lo, hi);
        latex.DrawLatex(0.15, 0.88, label.c_str());
        projections.push_back(h2);
    }

    axis->SetRange(0, 0);
    c->Update();
    SaveCanvas(c, saveName);
    for (TH2D* proj : projections) {
        delete proj;
    }
    delete c;
}

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
                             const std::vector<double>* metric_uncertainties,
                             double valueTextScale,
                             bool showBinIds) {
    if (!h3 || edges.size() < 2) return;
    const int nSlices = static_cast<int>(edges.size()) - 1;
    const int nCols = 2;
    const int nRows = (nSlices + nCols - 1) / nCols;

    TCanvas* c = new TCanvas(Form("c_phase_slices_metric_%d", sliceAxis), saveName, 2400, 2000);
    c->Divide(nCols, nRows, 0.002, 0.002);

    TAxis* axis = nullptr;
    if (sliceAxis == 1) axis = h3->GetXaxis();
    if (sliceAxis == 2) axis = h3->GetYaxis();
    if (sliceAxis == 3) axis = h3->GetZaxis();
    if (!axis) {
        delete c;
        return;
    }

    std::vector<TH2D*> projections;
    projections.reserve(nSlices);

    for (int i = 0; i < nSlices; ++i) {
        const double lo = edges[i];
        const double hi = edges[i + 1];
        const double tiny = 1e-12;
        const int binLo = axis->FindBin(lo + tiny);
        const int binHi = axis->FindBin(hi - tiny);
        axis->SetRange(binLo, binHi);

        TH2D* h2raw = (TH2D*)h3->Project3D(projOpt);
        if (!h2raw) continue;
        TH2D* h2 = (TH2D*)h2raw->Clone(Form("slice_metric_%s_%d", projOpt, i));
        h2->SetDirectory(nullptr);
        delete h2raw;
        h2->Reset("ICES");
        h2->SetTitle("");
        h2->GetXaxis()->SetTitle(xTitle);
        h2->GetYaxis()->SetTitle(yTitle);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.2);
        if (strcmp(xTitle, "#beta") == 0) {
            h2->GetXaxis()->SetRangeUser(0.0, 1.0);
        }
        if (strcmp(yTitle, "#beta") == 0) {
            h2->GetYaxis()->SetRangeUser(0.0, 1.0);
        }

        const double eps = 1e-12;
        const double xmin = h2->GetXaxis()->GetXmin();
        const double xmax = h2->GetXaxis()->GetXmax();
        const double ymin = h2->GetYaxis()->GetXmin();
        const double ymax = h2->GetYaxis()->GetXmax();
        for (const auto& b : overlayBins) {
            if (b.beta_max <= lo + eps || b.beta_min >= hi - eps) continue;
            double xlow = std::max(b.xpom_min, xmin);
            double xhigh = std::min(b.xpom_max, xmax);
            double ylow = std::max(b.Q2_min, ymin);
            double yhigh = std::min(b.Q2_max, ymax);
            if (xhigh <= xlow || yhigh <= ylow) continue;

            const int k = GetGlobalBinFromBinDef(b, q2_edges, xpom_edges, beta_edges);
            if (k < 0 || k >= static_cast<int>(metric_values.size())) continue;
            const double val = metric_values[k];
            if (!std::isfinite(val) || val < 0.0) continue;

            const int binxlow = h2->GetXaxis()->FindBin(xlow + eps);
            const int binxhigh = h2->GetXaxis()->FindBin(xhigh - eps);
            const int binylow = h2->GetYaxis()->FindBin(ylow + eps);
            const int binyhigh = h2->GetYaxis()->FindBin(yhigh - eps);
            for (int bx = binxlow; bx <= binxhigh; ++bx) {
                for (int by = binylow; by <= binyhigh; ++by) {
                    h2->SetBinContent(bx, by, val);
                }
            }
        }

        h2->SetMinimum(zMin);
        h2->SetMaximum(zMax);

        TPad* pad = (TPad*)c->cd(i + 1);
        pad->SetRightMargin(0.14);
        pad->SetLeftMargin(0.14);
        pad->SetTopMargin(0.12);
        pad->SetBottomMargin(0.12);
        if (logX) pad->SetLogx();
        if (logY) pad->SetLogy();
        h2->Draw("COLZ");

        DrawBinsForSliceMetric(h2,
                               overlayBins,
                               lo, hi,
                               overlayLogX, overlayLogY,
                               metric_values,
                               q2_edges, xpom_edges, beta_edges,
                               metric_uncertainties, valueTextScale,
                               showBinIds);

        TLatex latex;
        latex.SetTextSize(0.06);
        latex.SetNDC();
        latex.SetTextColor(kBlack);
        const std::string label = std::string(sliceLabel) + " " + FormatRange(lo, hi);
        latex.DrawLatex(0.15, 0.88, label.c_str());
        projections.push_back(h2);
    }

    axis->SetRange(0, 0);
    c->Update();
    SaveCanvas(c, saveName);
    for (TH2D* proj : projections) {
        delete proj;
    }
    delete c;
}
