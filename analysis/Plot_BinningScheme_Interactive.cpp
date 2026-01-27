//g++ -std=c++17 -o Plot_BinningScheme_Interactive Plot_BinningScheme_Interactive.cpp `root-config --cflags --libs`

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TPad.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <TMath.h>
#include <ctime>

//=========================================================================
// Structure to track bin modification history
//=========================================================================
struct BinModificationRecord {
    std::string timestamp;
    std::string variable;
    std::string operation;  // "merge" or "delete"
    double edgeValue;
    int binsBefore;
    int binsAfter;
    std::string reason;
};

// Global history vector
std::vector<BinModificationRecord> g_modificationHistory;

//=========================================================================
// Structure to represent a 3D bin
//=========================================================================
struct Bin3D {
    size_t i, j, k;  // Indices in the original grid
    double Q2_low, Q2_high;
    double beta_low, beta_high;
    double xpom_low, xpom_high;
    bool deleted;  // Mark if this bin has been deleted

    Bin3D(size_t ii, size_t jj, size_t kk,
          double q2l, double q2h,
          double bl, double bh,
          double xl, double xh)
        : i(ii), j(jj), k(kk),
          Q2_low(q2l), Q2_high(q2h),
          beta_low(bl), beta_high(bh),
          xpom_low(xl), xpom_high(xh),
          deleted(false) {}
};

// Global list of 3D bins (for irregular binning support)
std::vector<Bin3D> g_bins3D;

//=========================================================================
// Function to export modification history to CSV
//=========================================================================
void ExportModificationHistory(const std::string& filename = "bin_modification_history.csv") {
    if(g_modificationHistory.empty()) {
        std::cout << "No modifications to export." << std::endl;
        return;
    }

    std::ofstream outfile(filename);
    if(!outfile.is_open()) {
        std::cerr << "ERROR: Cannot open " << filename << " for writing." << std::endl;
        return;
    }

    // Write header
    outfile << "Timestamp,Variable,Operation,EdgeValue,BinsBefore,BinsAfter,Reason" << std::endl;

    // Write records
    for(const auto& record : g_modificationHistory) {
        outfile << record.timestamp << ","
                << record.variable << ","
                << record.operation << ","
                << record.edgeValue << ","
                << record.binsBefore << ","
                << record.binsAfter << ","
                << record.reason << std::endl;
    }

    outfile.close();
    std::cout << "\n✓ Modification history exported to: " << filename << std::endl;
    std::cout << "  Total operations: " << g_modificationHistory.size() << std::endl;
}

//=========================================================================
// Function to get current timestamp
//=========================================================================
std::string GetCurrentTimestamp() {
    time_t now = time(0);
    char buf[80];
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime(&now));
    return std::string(buf);
}

//=========================================================================
// Function to draw bin edges as text on plot
//=========================================================================
void DrawBinEdgesText(const std::vector<double>& bins, const std::string& varname,
                      double x, double y, double textsize = 0.02) {
    TPaveText* pt = new TPaveText(x, y, x + 0.25, y + 0.15, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetFillStyle(1001);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextSize(textsize);
    pt->SetTextFont(42);

    TString header = Form("%s bins (%d):", varname.c_str(), (int)(bins.size()-1));
    pt->AddText(header);

    // Show first few and last few edges if too many
    if(bins.size() <= 8) {
        TString edgesStr = "";
        for(size_t i = 0; i < bins.size(); i++) {
            if(TMath::Abs(bins[i]) < 0.01 || TMath::Abs(bins[i]) > 100) {
                edgesStr += Form("%.2e", bins[i]);
            } else {
                edgesStr += Form("%.3g", bins[i]);
            }
            if(i < bins.size()-1) edgesStr += ", ";
        }
        pt->AddText(edgesStr);
    } else {
        // Show first 3, ..., last 3
        TString edgesStr = "";
        for(int i = 0; i < 3; i++) {
            edgesStr += Form("%.3g, ", bins[i]);
        }
        edgesStr += "...";
        for(size_t i = bins.size()-3; i < bins.size(); i++) {
            edgesStr += Form(", %.3g", bins[i]);
        }
        pt->AddText(edgesStr);
    }

    pt->Draw();
}

//=========================================================================
// Function to draw binning grid with event counts
//=========================================================================
void DrawBinningGridWithCounts(TH2D* hist, std::vector<double> xbins, std::vector<double> ybins,
                                bool isLogX = false, double xmin = 0, double xmax = 1,
                                double ymin = 0, double ymax = 1, bool showBinEdges = true,
                                std::string xvarname = "x", std::string yvarname = "y",
                                int lowStatsThreshold = 10) {

    int lineColor = kRed;
    int lineWidth = 3;
    int lineStyle = 1;

    // Find actual data range (where bins are defined)
    double xmin_data = xbins.front();
    double xmax_data = xbins.back();
    double ymin_data = ybins.front();
    double ymax_data = ybins.back();

    if(isLogX) {
        xmin_data = TMath::Log10(xmin_data);
        xmax_data = TMath::Log10(xmax_data);
    }

    // Clip to histogram range
    xmin_data = TMath::Max(xmin_data, xmin);
    xmax_data = TMath::Min(xmax_data, xmax);
    ymin_data = TMath::Max(ymin_data, ymin);
    ymax_data = TMath::Min(ymax_data, ymax);

    // Draw vertical lines for x bins (clipped to data range)
    for(size_t i = 0; i < xbins.size(); i++) {
        double xval = xbins[i];
        if(isLogX) {
            xval = TMath::Log10(xval);
        }

        if(xval >= xmin && xval <= xmax) {
            TLine* line = new TLine(xval, ymin_data, xval, ymax_data);
            line->SetLineColor(lineColor);
            line->SetLineWidth(lineWidth);
            line->SetLineStyle(lineStyle);
            line->Draw("same");
        }
    }

    // Draw horizontal lines for y bins (clipped to data range)
    for(size_t i = 0; i < ybins.size(); i++) {
        double yval = ybins[i];

        if(yval >= ymin && yval <= ymax) {
            TLine* line = new TLine(xmin_data, yval, xmax_data, yval);
            line->SetLineColor(lineColor);
            line->SetLineWidth(lineWidth);
            line->SetLineStyle(lineStyle);
            line->Draw("same");
        }
    }

    // Add event counts in each bin
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.020);      // Decreased from 0.025 to 0.020
    latex->SetTextAlign(22);
    latex->SetTextFont(42);

    for(size_t i = 0; i < xbins.size() - 1; i++) {
        for(size_t j = 0; j < ybins.size() - 1; j++) {
            double xlow = xbins[i];
            double xhigh = xbins[i+1];
            double ylow = ybins[j];
            double yhigh = ybins[j+1];

            if(isLogX) {
                xlow = TMath::Log10(xbins[i]);
                xhigh = TMath::Log10(xbins[i+1]);
            }

            double xcenter = (xlow + xhigh) / 2.0;
            double ycenter = (ylow + yhigh) / 2.0;

            int count = 0;
            int binxlow = hist->GetXaxis()->FindBin(xlow);
            int binxhigh = hist->GetXaxis()->FindBin(xhigh);
            int binylow = hist->GetYaxis()->FindBin(ylow);
            int binyhigh = hist->GetYaxis()->FindBin(yhigh);

            for(int ix = binxlow; ix <= binxhigh; ix++) {
                for(int iy = binylow; iy <= binyhigh; iy++) {
                    count += hist->GetBinContent(ix, iy);
                }
            }

            if(count > 0) {
                // All numbers in black
                latex->SetTextColor(kBlack);
                TString countStr = Form("%d", count);
                latex->DrawLatex(xcenter, ycenter, countStr);
            }
        }
    }

    // Optionally draw bin edges as text
    if(showBinEdges) {
        DrawBinEdgesText(xbins, xvarname, 0.15, 0.75, 0.018);
        DrawBinEdgesText(ybins, yvarname, 0.15, 0.55, 0.018);
    }
}

//=========================================================================
// Function to draw binning grid with 3D event counts (minimum across 3rd dim)
//=========================================================================
void DrawBinningGridWithCounts3D(TH3D* hist3d, std::vector<double> xbins, std::vector<double> ybins,
                                  std::vector<double> zbins, int projectionType,
                                  bool isLogX = false, double xmin = 0, double xmax = 1,
                                  double ymin = 0, double ymax = 1, bool showBinEdges = true,
                                  std::string xvarname = "x", std::string yvarname = "y",
                                  int lowStatsThreshold = 10) {
    // projectionType: 0 = Q2 vs beta (iterate over xpom/z)
    //                 1 = Q2 vs xpom (iterate over beta/y)
    //                 2 = beta vs xpom (iterate over Q2/x)

    int lineColor = kRed;
    int lineWidth = 3;
    int lineStyle = 1;

    // Find actual data range
    double xmin_data = xbins.front();
    double xmax_data = xbins.back();
    double ymin_data = ybins.front();
    double ymax_data = ybins.back();

    if(isLogX) {
        xmin_data = TMath::Log10(xmin_data);
        xmax_data = TMath::Log10(xmax_data);
    }

    // Clip to histogram range
    xmin_data = TMath::Max(xmin_data, xmin);
    xmax_data = TMath::Min(xmax_data, xmax);
    ymin_data = TMath::Max(ymin_data, ymin);
    ymax_data = TMath::Min(ymax_data, ymax);

    // Draw grid lines
    for(size_t i = 0; i < xbins.size(); i++) {
        double xval = xbins[i];
        if(isLogX) xval = TMath::Log10(xval);
        if(xval >= xmin && xval <= xmax) {
            TLine* line = new TLine(xval, ymin_data, xval, ymax_data);
            line->SetLineColor(lineColor);
            line->SetLineWidth(lineWidth);
            line->SetLineStyle(lineStyle);
            line->Draw("same");
        }
    }

    for(size_t i = 0; i < ybins.size(); i++) {
        double yval = ybins[i];
        if(yval >= ymin && yval <= ymax) {
            TLine* line = new TLine(xmin_data, yval, xmax_data, yval);
            line->SetLineColor(lineColor);
            line->SetLineWidth(lineWidth);
            line->SetLineStyle(lineStyle);
            line->Draw("same");
        }
    }

    // Add event counts showing MINIMUM across 3rd dimension
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.020);
    latex->SetTextAlign(22);
    latex->SetTextFont(42);

    for(size_t i = 0; i < xbins.size() - 1; i++) {
        for(size_t j = 0; j < ybins.size() - 1; j++) {
            double xlow_orig = xbins[i];
            double xhigh_orig = xbins[i+1];
            double ylow = ybins[j];
            double yhigh = ybins[j+1];

            double xlow = xlow_orig;
            double xhigh = xhigh_orig;
            if(isLogX) {
                xlow = TMath::Log10(xlow_orig);
                xhigh = TMath::Log10(xhigh_orig);
            }

            double xcenter = (xlow + xhigh) / 2.0;
            double ycenter = (ylow + yhigh) / 2.0;

            // Find minimum event count across 3rd dimension
            int minCount = 999999;
            int totalCount = 0;

            for(size_t k = 0; k < zbins.size() - 1; k++) {
                double zlow = zbins[k];
                double zhigh = zbins[k+1];

                int binx, biny, binz;
                int count = 0;

                // Get bin indices based on projection type
                if(projectionType == 0) {  // Q2 vs beta (x=beta, y=Q2, z=xpom)
                    binx = hist3d->GetXaxis()->FindBin((ylow + yhigh) / 2.0);  // Q2
                    biny = hist3d->GetYaxis()->FindBin((xlow_orig + xhigh_orig) / 2.0);  // beta
                    binz = hist3d->GetZaxis()->FindBin((zlow + zhigh) / 2.0);  // xpom
                } else if(projectionType == 1) {  // Q2 vs xpom (x=xpom, y=Q2, z=beta)
                    binx = hist3d->GetXaxis()->FindBin((ylow + yhigh) / 2.0);  // Q2
                    biny = hist3d->GetYaxis()->FindBin((zlow + zhigh) / 2.0);  // beta
                    binz = hist3d->GetZaxis()->FindBin((xlow_orig + xhigh_orig) / 2.0);  // xpom
                } else {  // beta vs xpom (x=xpom, y=beta, z=Q2)
                    binx = hist3d->GetXaxis()->FindBin((zlow + zhigh) / 2.0);  // Q2
                    biny = hist3d->GetYaxis()->FindBin((ylow + yhigh) / 2.0);  // beta
                    binz = hist3d->GetZaxis()->FindBin((xlow_orig + xhigh_orig) / 2.0);  // xpom
                }

                count = hist3d->GetBinContent(binx, biny, binz);
                totalCount += count;
                if(count > 0 && count < minCount) minCount = count;
            }

            if(minCount < 999999 && minCount > 0) {
                // All numbers in black
                latex->SetTextColor(kBlack);
                TString countStr = Form("%d", minCount);
                latex->DrawLatex(xcenter, ycenter, countStr);
            }
        }
    }

    if(showBinEdges) {
        DrawBinEdgesText(xbins, xvarname, 0.15, 0.75, 0.018);
        DrawBinEdgesText(ybins, yvarname, 0.15, 0.55, 0.018);
    }
}

//=========================================================================
// Display bin edges nicely
//=========================================================================
void DisplayBinEdges(const std::vector<double>& bins, const std::string& varname) {
    std::cout << "\n  " << varname << " (" << bins.size()-1 << " bins):" << std::endl;
    std::cout << "  [";
    for(size_t i = 0; i < bins.size(); i++) {
        if(TMath::Abs(bins[i]) < 0.001 || TMath::Abs(bins[i]) > 1000) {
            std::cout << std::scientific << std::setprecision(3) << bins[i];
        } else {
            std::cout << std::fixed << std::setprecision(4) << bins[i];
        }
        if(i < bins.size()-1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

//=========================================================================
// Show 3D bin statistics
//=========================================================================
void Show3DBinStatistics(TH3D* hist3d,
                         const std::vector<double>& Q2_bins,
                         const std::vector<double>& beta_bins,
                         const std::vector<double>& xpom_bins) {

    std::cout << "\n========================================" << std::endl;
    std::cout << "3D Bin Statistics (Q², β, x_pom)" << std::endl;
    std::cout << "========================================" << std::endl;

    // Count non-deleted bins
    int total_bins_count = 0;
    for(const auto& bin : g_bins3D) {
        if(!bin.deleted) total_bins_count++;
    }

    std::cout << "\nBinning:" << std::endl;
    std::cout << "  Q²: " << Q2_bins.size()-1 << " bins" << std::endl;
    std::cout << "  β: " << beta_bins.size()-1 << " bins" << std::endl;
    std::cout << "  x_pom: " << xpom_bins.size()-1 << " bins" << std::endl;
    std::cout << "  Total 3D bins (active): " << total_bins_count << std::endl;

    std::cout << "\n" << std::string(100, '=') << std::endl;
    std::cout << std::setw(8) << "Bin #"
              << std::setw(20) << "Q² range"
              << std::setw(20) << "β range"
              << std::setw(25) << "x_pom range"
              << std::setw(15) << "Events" << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    int binNum = 0;
    int total_events = 0;
    int non_empty_bins = 0;

    // Iterate over g_bins3D instead of the regular grid to respect deleted bins
    for(const auto& bin3d : g_bins3D) {
        // Skip deleted bins
        if(bin3d.deleted) {
            continue;
        }

        binNum++;

        double Q2_low = bin3d.Q2_low;
        double Q2_high = bin3d.Q2_high;
        double beta_low = bin3d.beta_low;
        double beta_high = bin3d.beta_high;
        double xpom_low = bin3d.xpom_low;
        double xpom_high = bin3d.xpom_high;

        // Find all histogram bins that overlap with this custom bin
        int binx_low = hist3d->GetXaxis()->FindBin(Q2_low);
        int binx_high = hist3d->GetXaxis()->FindBin(Q2_high);
        int biny_low = hist3d->GetYaxis()->FindBin(beta_low);
        int biny_high = hist3d->GetYaxis()->FindBin(beta_high);
        int binz_low = hist3d->GetZaxis()->FindBin(xpom_low);
        int binz_high = hist3d->GetZaxis()->FindBin(xpom_high);

        // Sum events in all histogram bins within this custom bin
        int count = 0;
        for(int ix = binx_low; ix <= binx_high; ix++) {
            for(int iy = biny_low; iy <= biny_high; iy++) {
                for(int iz = binz_low; iz <= binz_high; iz++) {
                    count += hist3d->GetBinContent(ix, iy, iz);
                }
            }
        }

        total_events += count;

        if(count > 0) {
            non_empty_bins++;
            std::cout << std::setw(8) << binNum
                      << std::setw(9) << std::fixed << std::setprecision(1) << Q2_low
                      << "-" << std::setw(8) << Q2_high
                      << std::setw(9) << std::fixed << std::setprecision(2) << beta_low
                      << "-" << std::setw(8) << beta_high
                      << std::setw(12) << std::scientific << std::setprecision(2) << xpom_low
                      << "-" << std::setw(10) << xpom_high
                      << std::setw(15) << count << std::endl;
        }
    }
    std::cout << std::string(100, '=') << std::endl;
    std::cout << "Total events: " << total_events << std::endl;
    std::cout << "Non-empty bins: " << non_empty_bins << " / " << binNum << std::endl;
}

//=========================================================================
// Forward declaration for CreatePlots
//=========================================================================
void CreatePlots(TFile* inputFile, TH3D* h_d3sigma,
                const std::vector<double>& Q2_bins,
                const std::vector<double>& beta_bins,
                const std::vector<double>& xpom_bins,
                const std::vector<double>& y_bins,
                bool saveFinal);

//=========================================================================
// Show compact 3D bin statistics summary
//=========================================================================
void Show3DBinSummary(TH3D* hist3d,
                      const std::vector<double>& Q2_bins,
                      const std::vector<double>& beta_bins,
                      const std::vector<double>& xpom_bins) {

    int total_3d_bins = 0;  // Count only non-deleted bins
    int non_empty_bins = 0;
    int total_events = 0;
    int min_event_count = 999999;
    int max_event_count = 0;

    // Calculate statistics from g_bins3D (respecting deleted bins)
    for(const auto& bin3d : g_bins3D) {
        // Skip deleted bins
        if(bin3d.deleted) {
            continue;
        }

        total_3d_bins++;

        double Q2_low = bin3d.Q2_low;
        double Q2_high = bin3d.Q2_high;
        double beta_low = bin3d.beta_low;
        double beta_high = bin3d.beta_high;
        double xpom_low = bin3d.xpom_low;
        double xpom_high = bin3d.xpom_high;

        // Find all histogram bins that overlap with this custom bin
        int binx_low = hist3d->GetXaxis()->FindBin(Q2_low);
        int binx_high = hist3d->GetXaxis()->FindBin(Q2_high);
        int biny_low = hist3d->GetYaxis()->FindBin(beta_low);
        int biny_high = hist3d->GetYaxis()->FindBin(beta_high);
        int binz_low = hist3d->GetZaxis()->FindBin(xpom_low);
        int binz_high = hist3d->GetZaxis()->FindBin(xpom_high);

        // Sum events in all histogram bins within this custom bin
        int count = 0;
        for(int ix = binx_low; ix <= binx_high; ix++) {
            for(int iy = biny_low; iy <= biny_high; iy++) {
                for(int iz = binz_low; iz <= binz_high; iz++) {
                    count += hist3d->GetBinContent(ix, iy, iz);
                }
            }
        }

        total_events += count;
        if(count > 0) {
            non_empty_bins++;
            if(count < min_event_count) min_event_count = count;
            if(count > max_event_count) max_event_count = count;
        }
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "3D Binning Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Q²: " << Q2_bins.size()-1 << " bins  ["
              << std::fixed << std::setprecision(1) << Q2_bins.front() << " - "
              << Q2_bins.back() << " GeV²]" << std::endl;
    std::cout << "  β: " << beta_bins.size()-1 << " bins  ["
              << std::fixed << std::setprecision(3) << beta_bins.front() << " - "
              << beta_bins.back() << "]" << std::endl;
    std::cout << "  x_pom: " << xpom_bins.size()-1 << " bins  ["
              << std::scientific << std::setprecision(2) << xpom_bins.front() << " - "
              << xpom_bins.back() << "]" << std::endl;
    std::cout << "\n  Total 3D bins: " << total_3d_bins << std::endl;
    std::cout << "  Non-empty bins: " << non_empty_bins << " ("
              << std::fixed << std::setprecision(1)
              << 100.0 * non_empty_bins / total_3d_bins << "%)" << std::endl;
    std::cout << "  Empty bins: " << total_3d_bins - non_empty_bins << std::endl;
    std::cout << "\n  Total events: " << total_events << std::endl;
    if(non_empty_bins > 0) {
        std::cout << "  Min events/bin: " << min_event_count;
        if(min_event_count < 100) std::cout << " ⚠ LOW";
        std::cout << std::endl;
        std::cout << "  Max events/bin: " << max_event_count << std::endl;
        std::cout << "  Avg events/bin: " << std::fixed << std::setprecision(0)
                  << (double)total_events / non_empty_bins << std::endl;
    }
    std::cout << "========================================" << std::endl;
}

//=========================================================================
// Calculate and display bin statistics
//=========================================================================
void ShowBinStatistics(TH2D* hist, const std::vector<double>& xbins,
                       const std::vector<double>& ybins, bool isLogX,
                       const std::string& xname, const std::string& yname) {
    std::cout << "\n  Bin Statistics for " << xname << " vs " << yname << ":" << std::endl;
    std::cout << "  " << std::string(60, '-') << std::endl;
    std::cout << "  Bin #  | " << std::setw(15) << xname << " range | "
              << std::setw(15) << yname << " range | Events" << std::endl;
    std::cout << "  " << std::string(60, '-') << std::endl;

    int total_events = 0;
    int bin_num = 0;

    for(size_t j = 0; j < ybins.size() - 1; j++) {
        for(size_t i = 0; i < xbins.size() - 1; i++) {
            double xlow_orig = xbins[i];
            double xhigh_orig = xbins[i+1];
            double ylow = ybins[j];
            double yhigh = ybins[j+1];

            double xlow = isLogX ? TMath::Log10(xbins[i]) : xbins[i];
            double xhigh = isLogX ? TMath::Log10(xbins[i+1]) : xbins[i+1];

            int count = 0;
            int binxlow = hist->GetXaxis()->FindBin(xlow);
            int binxhigh = hist->GetXaxis()->FindBin(xhigh);
            int binylow = hist->GetYaxis()->FindBin(ylow);
            int binyhigh = hist->GetYaxis()->FindBin(yhigh);

            for(int ix = binxlow; ix <= binxhigh; ix++) {
                for(int iy = binylow; iy <= binyhigh; iy++) {
                    count += hist->GetBinContent(ix, iy);
                }
            }

            if(count > 0) {
                bin_num++;
                std::cout << "  " << std::setw(5) << bin_num << "  | ";
                std::cout << std::setw(6) << std::setprecision(3) << xlow_orig << "-"
                          << std::setw(6) << xhigh_orig << " | ";
                std::cout << std::setw(6) << ylow << "-" << std::setw(6) << yhigh << " | ";
                std::cout << std::setw(6) << count << std::endl;
                total_events += count;
            }
        }
    }
    std::cout << "  " << std::string(60, '-') << std::endl;
    std::cout << "  Total events: " << total_events << std::endl;
}

//=========================================================================
// Safe integer input with error handling
//=========================================================================
int GetIntInput(const std::string& prompt = "") {
    int value;
    while(true) {
        if(!prompt.empty()) {
            std::cout << prompt;
        }

        if(std::cin >> value) {
            // Clear any remaining input
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return value;
        } else {
            // Clear error state
            std::cin.clear();
            // Discard bad input
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter a number: ";
        }
    }
}

//=========================================================================
// Generate logarithmic bins
//=========================================================================
std::vector<double> GenerateLogBins(int nbins, double min, double max) {
    std::vector<double> bins;
    double logmin = TMath::Log10(min);
    double logmax = TMath::Log10(max);
    double step = (logmax - logmin) / nbins;

    for(int i = 0; i <= nbins; i++) {
        bins.push_back(TMath::Power(10, logmin + i * step));
    }
    return bins;
}

//=========================================================================
// Generate linear bins
//=========================================================================
std::vector<double> GenerateLinearBins(int nbins, double min, double max) {
    std::vector<double> bins;
    double step = (max - min) / nbins;

    for(int i = 0; i <= nbins; i++) {
        bins.push_back(min + i * step);
    }
    return bins;
}

//=========================================================================
// Read manual bin edges from user
//=========================================================================
std::vector<double> ReadManualBins(const std::string& varname) {
    std::vector<double> bins;
    std::string input;

    std::cout << "\nEnter bin edges for " << varname << " (space-separated, in ascending order):" << std::endl;
    std::cout << "Example: 1.0 2.0 5.0 10.0 20.0" << std::endl;
    std::cout << "> ";

    std::getline(std::cin, input);

    std::istringstream iss(input);
    double value;
    while(iss >> value) {
        bins.push_back(value);
    }

    if(bins.size() < 2) {
        std::cerr << "Error: Need at least 2 bin edges. Using default." << std::endl;
        bins.clear();
    }

    return bins;
}

//=========================================================================
// Merge/Delete bin edges
//=========================================================================
std::vector<double> MergeDeleteBins(const std::vector<double>& current_bins, const std::string& varname) {
    std::vector<double> bins = current_bins;

    while(true) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Merge/Delete Bin Edges: " << varname << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Current bin edges (" << bins.size() << " edges, " << bins.size()-1 << " bins):" << std::endl;

        for(size_t i = 0; i < bins.size(); i++) {
            std::cout << "  [" << i << "] " << std::setw(12) << std::setprecision(6) << bins[i];
            if(i > 0 && i < bins.size()-1) {
                std::cout << "  (remove to merge bins " << i-1 << " and " << i << ")";
                std::cout << "  Range: [" << bins[i-1] << " - " << bins[i] << "] + ["
                          << bins[i] << " - " << bins[i+1] << "]";
            } else if(i == 0) {
                std::cout << "  (lower boundary - remove to delete lowest bin)";
            } else {
                std::cout << "  (upper boundary - remove to delete highest bin)";
            }
            std::cout << std::endl;
        }

        std::cout << "\nOptions:" << std::endl;
        std::cout << "  Enter edge index to remove (0-" << bins.size()-1 << ")" << std::endl;
        std::cout << "  Or -1 to finish and return" << std::endl;

        int choice = GetIntInput("Choice: ");

        if(choice == -1) {
            break;
        } else if(choice < 0 || choice >= (int)bins.size()) {
            std::cout << "Invalid choice." << std::endl;
        } else if(bins.size() <= 2) {
            std::cout << "Cannot remove more edges. Minimum 1 bin required." << std::endl;
        } else {
            // Show what will happen
            std::cout << "\n⚠ Confirmation Required" << std::endl;
            std::cout << "  Edge to remove: " << bins[choice] << std::endl;

            if(choice == 0) {
                // Deleting lower boundary - removes lowest bin
                std::cout << "  This will DELETE the lowest bin:" << std::endl;
                std::cout << "    Bin 0: [" << bins[0] << " - " << bins[1] << "] will be DELETED" << std::endl;
                std::cout << "  New range will be: [" << bins[1] << " - " << bins[bins.size()-1] << "]" << std::endl;
            } else if(choice == (int)bins.size()-1) {
                // Deleting upper boundary - removes highest bin
                std::cout << "  This will DELETE the highest bin:" << std::endl;
                std::cout << "    Bin " << bins.size()-2 << ": [" << bins[bins.size()-2] << " - " << bins[bins.size()-1] << "] will be DELETED" << std::endl;
                std::cout << "  New range will be: [" << bins[0] << " - " << bins[bins.size()-2] << "]" << std::endl;
            } else {
                // Merging two adjacent bins
                std::cout << "  This will merge:" << std::endl;
                std::cout << "    Bin " << choice-1 << ": [" << bins[choice-1] << " - " << bins[choice] << "]" << std::endl;
                std::cout << "    Bin " << choice << ": [" << bins[choice] << " - " << bins[choice+1] << "]" << std::endl;
                std::cout << "  Into new bin: [" << bins[choice-1] << " - " << bins[choice+1] << "]" << std::endl;
            }
            std::cout << "  Total bins will change: " << bins.size()-1 << " → " << bins.size()-2 << std::endl;

            std::cout << "\nProceed with removal? (1=Yes, 0=No): ";
            int confirm = GetIntInput("");

            if(confirm == 1) {
                // Ask for optional reason
                std::cout << "Optional: Enter reason for this change (or press Enter to skip): ";
                std::string reason;
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::getline(std::cin, reason);

                // Determine operation type
                std::string operationType;
                if(choice == 0) {
                    operationType = "delete_lower";
                    if(reason.empty()) reason = "Deleted lowest bin";
                } else if(choice == (int)bins.size()-1) {
                    operationType = "delete_upper";
                    if(reason.empty()) reason = "Deleted highest bin";
                } else {
                    operationType = "merge";
                    if(reason.empty()) reason = "User-initiated merge";
                }

                // Record the operation
                BinModificationRecord record;
                record.timestamp = GetCurrentTimestamp();
                record.variable = varname;
                record.operation = operationType;
                record.edgeValue = bins[choice];
                record.binsBefore = bins.size() - 1;
                record.binsAfter = bins.size() - 2;
                record.reason = reason;
                g_modificationHistory.push_back(record);

                // Perform the removal
                if(choice == 0) {
                    std::cout << "✓ Deleting lowest bin by removing edge at " << bins[choice] << std::endl;
                } else if(choice == (int)bins.size()-1) {
                    std::cout << "✓ Deleting highest bin by removing edge at " << bins[choice] << std::endl;
                } else {
                    std::cout << "✓ Merging bins by removing edge at " << bins[choice] << std::endl;
                }
                bins.erase(bins.begin() + choice);
                std::cout << "✓ New configuration: " << bins.size()-1 << " bins" << std::endl;
                std::cout << "✓ Operation logged to history" << std::endl;
            } else {
                std::cout << "✗ Operation cancelled" << std::endl;
            }
        }
    }

    return bins;
}

//=========================================================================
// Get binning choice from user
//=========================================================================
std::vector<double> GetBinningChoice(const std::string& varname,
                                     const std::vector<double>& current_bins,
                                     TH2D* histogram, bool isXaxis) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Configure binning for: " << varname << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Current binning: " << current_bins.size() - 1 << " bins" << std::endl;
    std::cout << "  Range: [" << current_bins.front() << ", " << current_bins.back() << "]" << std::endl;

    DisplayBinEdges(current_bins, varname);

    std::cout << "\nOptions:" << std::endl;
    std::cout << "  1. Keep current binning" << std::endl;
    std::cout << "  2. Linear binning (equal width)" << std::endl;
    std::cout << "  3. Logarithmic binning" << std::endl;
    std::cout << "  4. Manual bin edges" << std::endl;
    std::cout << "  5. Merge/Delete bin edges" << std::endl;

    int choice = GetIntInput("Enter choice (1-5): ");

    if(choice == 1) {
        return current_bins;
    }

    // Get range from histogram
    TAxis* axis = isXaxis ? histogram->GetXaxis() : histogram->GetYaxis();
    double hmin = axis->GetXmin();
    double hmax = axis->GetXmax();

    if(choice == 2) {
        int nbins = GetIntInput("Enter number of bins: ");
        return GenerateLinearBins(nbins, hmin, hmax);
    }
    else if(choice == 3) {
        int nbins = GetIntInput("Enter number of bins: ");

        double min, max;
        std::cout << "Enter min value (>0 for log scale): ";
        std::cin >> min;
        std::cout << "Enter max value: ";
        std::cin >> max;

        if(min <= 0) {
            std::cerr << "Error: Min value must be > 0 for log scale. Using histogram range." << std::endl;
            min = TMath::Max(1e-6, hmin);
        }
        return GenerateLogBins(nbins, min, max);
    }
    else if(choice == 4) {
        std::vector<double> manual_bins = ReadManualBins(varname);
        if(manual_bins.empty()) {
            return current_bins;
        }
        return manual_bins;
    }
    else if(choice == 5) {
        return MergeDeleteBins(current_bins, varname);
    }
    else {
        std::cout << "Invalid choice. Keeping current binning." << std::endl;
        return current_bins;
    }
}

//=========================================================================
// Function to rebuild 3D bins from the regular grid
//=========================================================================
void RebuildBins3D(const std::vector<double>& Q2_bins,
                   const std::vector<double>& beta_bins,
                   const std::vector<double>& xpom_bins) {
    g_bins3D.clear();

    for(size_t i = 0; i < Q2_bins.size() - 1; i++) {
        for(size_t j = 0; j < beta_bins.size() - 1; j++) {
            for(size_t k = 0; k < xpom_bins.size() - 1; k++) {
                g_bins3D.emplace_back(i, j, k,
                                      Q2_bins[i], Q2_bins[i+1],
                                      beta_bins[j], beta_bins[j+1],
                                      xpom_bins[k], xpom_bins[k+1]);
            }
        }
    }

    std::cout << "Rebuilt 3D bins: " << g_bins3D.size() << " bins total" << std::endl;
}

//=========================================================================
// Export bin edges to file (now includes 3D bins)
//=========================================================================
void ExportBinEdges(const std::string& filename,
                    TH3D* hist3d,
                    const std::vector<double>& Q2_bins,
                    const std::vector<double>& beta_bins,
                    const std::vector<double>& xpom_bins,
                    const std::vector<double>& y_bins) {

    // Ask user if they want to include empty bins
    std::cout << "\nInclude empty bins in 3D bin export? (y/n): ";
    char includeEmpty;
    std::cin >> includeEmpty;
    bool include_empty = (includeEmpty == 'y' || includeEmpty == 'Y');

    std::ofstream outfile(filename);

    if(!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }

    outfile << "# Binning scheme export" << std::endl;
    outfile << "# Generated by Plot_BinningScheme_Interactive" << std::endl;
    outfile << std::endl;

    // Q² bins
    outfile << "Q2_bins," << Q2_bins.size() - 1 << std::endl;
    for(size_t i = 0; i < Q2_bins.size(); i++) {
        outfile << std::setprecision(10) << Q2_bins[i];
        if(i < Q2_bins.size() - 1) outfile << ",";
    }
    outfile << std::endl << std::endl;

    // β bins
    outfile << "beta_bins," << beta_bins.size() - 1 << std::endl;
    for(size_t i = 0; i < beta_bins.size(); i++) {
        outfile << std::setprecision(10) << beta_bins[i];
        if(i < beta_bins.size() - 1) outfile << ",";
    }
    outfile << std::endl << std::endl;

    // x_pom bins
    outfile << "xpom_bins," << xpom_bins.size() - 1 << std::endl;
    for(size_t i = 0; i < xpom_bins.size(); i++) {
        outfile << std::setprecision(10) << xpom_bins[i];
        if(i < xpom_bins.size() - 1) outfile << ",";
    }
    outfile << std::endl << std::endl;

    // y bins
    outfile << "y_bins," << y_bins.size() - 1 << std::endl;
    for(size_t i = 0; i < y_bins.size(); i++) {
        outfile << std::setprecision(10) << y_bins[i];
        if(i < y_bins.size() - 1) outfile << ",";
    }
    outfile << std::endl << std::endl;

    // 3D bins (complete list of all bins with their boundaries)
    if(include_empty) {
        outfile << "# 3D Bin Structure (non-deleted bins - all bins including empty)" << std::endl;
    } else {
        outfile << "# 3D Bin Structure (non-deleted bins - only bins with events)" << std::endl;
    }
    outfile << "bin_number,Q2_low,Q2_high,beta_low,beta_high,xpom_low,xpom_high,event_count" << std::endl;

    int bin_number = 0;
    int written_bins = 0;
    for(const auto& bin : g_bins3D) {
        if(bin.deleted) {
            continue;
        }

        // Calculate event count for this bin
        int binx_low = hist3d->GetXaxis()->FindBin(bin.Q2_low);
        int binx_high = hist3d->GetXaxis()->FindBin(bin.Q2_high);
        int biny_low = hist3d->GetYaxis()->FindBin(bin.beta_low);
        int biny_high = hist3d->GetYaxis()->FindBin(bin.beta_high);
        int binz_low = hist3d->GetZaxis()->FindBin(bin.xpom_low);
        int binz_high = hist3d->GetZaxis()->FindBin(bin.xpom_high);

        int count = 0;
        for(int ix = binx_low; ix <= binx_high; ix++) {
            for(int iy = biny_low; iy <= biny_high; iy++) {
                for(int iz = binz_low; iz <= binz_high; iz++) {
                    count += hist3d->GetBinContent(ix, iy, iz);
                }
            }
        }

        // Skip empty bins if user chose not to include them
        if(!include_empty && count == 0) {
            continue;
        }

        bin_number++;
        written_bins++;
        outfile << bin_number << ","
               << std::setprecision(10) << bin.Q2_low << ","
               << std::setprecision(10) << bin.Q2_high << ","
               << std::setprecision(10) << bin.beta_low << ","
               << std::setprecision(10) << bin.beta_high << ","
               << std::setprecision(10) << bin.xpom_low << ","
               << std::setprecision(10) << bin.xpom_high << ","
               << count << std::endl;
    }

    outfile.close();
    std::cout << "Bin edges exported to: " << filename << std::endl;
    std::cout << "  - 1D bin edges (Q2, beta, xpom, y)" << std::endl;
    std::cout << "  - 3D bin structure: " << written_bins << " bins written";
    if(!include_empty) {
        std::cout << " (non-empty only)";
    }
    std::cout << std::endl;
}

//=========================================================================
// Display all current bin edges
//=========================================================================
void DisplayAllBinEdges(const std::vector<double>& Q2_bins,
                        const std::vector<double>& beta_bins,
                        const std::vector<double>& xpom_bins,
                        const std::vector<double>& y_bins) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Current Binning Configuration" << std::endl;
    std::cout << "========================================" << std::endl;

    DisplayBinEdges(Q2_bins, "Q²");
    DisplayBinEdges(beta_bins, "β");
    DisplayBinEdges(xpom_bins, "x_pom");
    DisplayBinEdges(y_bins, "y");

    std::cout << "\nPress Enter to continue...";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cin.get();
}

//=========================================================================
// Merge individual 3D bin with adjacent bin (upward or downward)
//=========================================================================
void MergeIndividual3DBin(TFile* inputFile, TH3D* hist3d,
                          std::vector<double>& Q2_bins,
                          std::vector<double>& beta_bins,
                          std::vector<double>& xpom_bins,
                          const std::vector<double>& y_bins) {

    const int LOW_STATS_THRESHOLD = 100;

    struct BinInfo {
        int binNum;
        size_t i, j, k;
        double Q2_low, Q2_high, beta_low, beta_high, xpom_low, xpom_high;
        int count;
    };

    bool done = false;

    while(!done) {
        // Step 1: Find and display all low-statistics 3D bins
        std::cout << "\n========================================" << std::endl;
        std::cout << "Low-Statistics 3D Bins (< " << LOW_STATS_THRESHOLD << " events)" << std::endl;
        std::cout << "========================================" << std::endl;

        std::vector<BinInfo> lowStatsBins;
        int binNum = 0;

    // Iterate over g_bins3D instead of the regular grid to respect deleted bins
    for(const auto& bin3d : g_bins3D) {
        // Skip deleted bins
        if(bin3d.deleted) {
            continue;
        }

        binNum++;

        // Calculate event count for this 3D bin
        double Q2_low = bin3d.Q2_low;
        double Q2_high = bin3d.Q2_high;
        double beta_low = bin3d.beta_low;
        double beta_high = bin3d.beta_high;
        double xpom_low = bin3d.xpom_low;
        double xpom_high = bin3d.xpom_high;

        int binx_low = hist3d->GetXaxis()->FindBin(Q2_low);
        int binx_high = hist3d->GetXaxis()->FindBin(Q2_high);
        int biny_low = hist3d->GetYaxis()->FindBin(beta_low);
        int biny_high = hist3d->GetYaxis()->FindBin(beta_high);
        int binz_low = hist3d->GetZaxis()->FindBin(xpom_low);
        int binz_high = hist3d->GetZaxis()->FindBin(xpom_high);

        int count = 0;
        for(int ix = binx_low; ix <= binx_high; ix++) {
            for(int iy = biny_low; iy <= biny_high; iy++) {
                for(int iz = binz_low; iz <= binz_high; iz++) {
                    count += hist3d->GetBinContent(ix, iy, iz);
                }
            }
        }

        if(count < LOW_STATS_THRESHOLD && count > 0) {
            BinInfo info;
            info.binNum = binNum;
            info.i = bin3d.i;
            info.j = bin3d.j;
            info.k = bin3d.k;
            info.Q2_low = Q2_low;
            info.Q2_high = Q2_high;
            info.beta_low = beta_low;
            info.beta_high = beta_high;
            info.xpom_low = xpom_low;
            info.xpom_high = xpom_high;
            info.count = count;
            lowStatsBins.push_back(info);
        }
    }

    if(lowStatsBins.empty()) {
        std::cout << "\n  ✓ No low-statistics bins found!" << std::endl;
        std::cout << "  All 3D bins have >= " << LOW_STATS_THRESHOLD << " events." << std::endl;
        std::cout << "\nReturning to workshop menu...";
        break;  // Exit the loop and return to workshop
    }

    // Display low-stats bins
    std::cout << "\n  Bin #  | Q² range      | β range      | x_pom range   | Events" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;

    for(const auto& bin : lowStatsBins) {
        std::cout << "  " << std::setw(5) << bin.binNum << "  | "
                  << std::setw(5) << std::fixed << std::setprecision(1) << bin.Q2_low << "-"
                  << std::setw(6) << bin.Q2_high << " | "
                  << std::setw(4) << std::setprecision(2) << bin.beta_low << "-"
                  << std::setw(4) << bin.beta_high << " | "
                  << std::setw(5) << std::setprecision(3) << bin.xpom_low << "-"
                  << std::setw(5) << bin.xpom_high << " | "
                  << std::setw(6) << bin.count << " ⚠" << std::endl;
    }

    // Step 2: User selects a bin
    std::cout << "\nEnter bin number to merge/delete (or 0 to return to workshop): ";
    int selectedBinNum = GetIntInput("");

    if(selectedBinNum == 0) {
        std::cout << "Returning to workshop menu..." << std::endl;
        break;  // Exit to workshop menu
    }

    // Find the selected bin
    BinInfo* selectedBin = nullptr;
    for(auto& bin : lowStatsBins) {
        if(bin.binNum == selectedBinNum) {
            selectedBin = &bin;
            break;
        }
    }

    if(!selectedBin) {
        std::cout << "Invalid bin number. Try again." << std::endl;
        continue;  // Go back to start of loop
    }

    // Step 3: Display selected bin details
    std::cout << "\n========================================" << std::endl;
    std::cout << "Selected Bin #" << selectedBin->binNum << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Q²: [" << selectedBin->Q2_low << " - " << selectedBin->Q2_high << "]" << std::endl;
    std::cout << "  β: [" << selectedBin->beta_low << " - " << selectedBin->beta_high << "]" << std::endl;
    std::cout << "  x_pom: [" << selectedBin->xpom_low << " - " << selectedBin->xpom_high << "]" << std::endl;
    std::cout << "  Events: " << selectedBin->count << std::endl;

    // Step 4: Choose merge dimension
    std::cout << "\nChoose dimension to merge along:" << std::endl;
    std::cout << "  1. Merge along Q² (keep β and x_pom)" << std::endl;
    std::cout << "  2. Merge along β (keep Q² and x_pom)" << std::endl;
    std::cout << "  3. Merge along x_pom (keep Q² and β)" << std::endl;
    std::cout << "  0. Cancel" << std::endl;

    int dimChoice = GetIntInput("Enter choice (0-3): ");

    if(dimChoice == 0) {
        std::cout << "Cancelled." << std::endl;
        continue;  // Go back to bin selection
    }

    if(dimChoice < 1 || dimChoice > 3) {
        std::cout << "Invalid choice. Try again." << std::endl;
        continue;  // Go back to start
    }

    // Step 5: Choose merge/delete action
    std::cout << "\nChoose action:" << std::endl;
    std::cout << "  1. Merge UPWARD (with higher bin value)" << std::endl;
    std::cout << "  2. Merge DOWNWARD (with lower bin value)" << std::endl;
    std::cout << "  3. DELETE this bin (remove from binning entirely)" << std::endl;
    std::cout << "  0. Cancel" << std::endl;

    int dirChoice = GetIntInput("Enter choice (0-3): ");

    if(dirChoice == 0) {
        std::cout << "Cancelled." << std::endl;
        continue;  // Go back to bin selection
    }

    if(dirChoice < 1 || dirChoice > 3) {
        std::cout << "Invalid choice. Try again." << std::endl;
        continue;  // Go back to start
    }

    // Step 6: Perform the merge
    bool success = false;
    std::string mergeInfo;

    if(dimChoice == 1) {  // Merge along Q²
        if(dirChoice == 1) {  // Upward - remove upper edge
            if(selectedBin->i < Q2_bins.size() - 2) {
                double edgeToRemove = Q2_bins[selectedBin->i + 1];
                mergeInfo = "Merging Q² bins [" + std::to_string(selectedBin->Q2_low) + "-" +
                           std::to_string(selectedBin->Q2_high) + "] and [" +
                           std::to_string(selectedBin->Q2_high) + "-" +
                           std::to_string(Q2_bins[selectedBin->i + 2]) + "]";
                Q2_bins.erase(Q2_bins.begin() + selectedBin->i + 1);
                success = true;
            } else {
                std::cout << "Cannot merge upward - already at highest Q² bin." << std::endl;
            }
        } else if(dirChoice == 2) {  // Downward - remove lower edge
            if(selectedBin->i > 0) {
                double edgeToRemove = Q2_bins[selectedBin->i];
                mergeInfo = "Merging Q² bins [" + std::to_string(Q2_bins[selectedBin->i - 1]) + "-" +
                           std::to_string(selectedBin->Q2_low) + "] and [" +
                           std::to_string(selectedBin->Q2_low) + "-" +
                           std::to_string(selectedBin->Q2_high) + "]";
                Q2_bins.erase(Q2_bins.begin() + selectedBin->i);
                success = true;
            } else {
                std::cout << "Cannot merge downward - already at lowest Q² bin." << std::endl;
            }
        } else {  // Delete - mark this specific 3D bin as deleted
            // Find and mark the specific 3D bin as deleted
            for(auto& bin : g_bins3D) {
                if(bin.i == selectedBin->i && bin.j == selectedBin->j && bin.k == selectedBin->k) {
                    if(!bin.deleted) {
                        bin.deleted = true;
                        mergeInfo = "Deleted 3D bin: Q²[" + std::to_string(selectedBin->Q2_low) + "-" +
                                   std::to_string(selectedBin->Q2_high) + "], β[" +
                                   std::to_string(selectedBin->beta_low) + "-" +
                                   std::to_string(selectedBin->beta_high) + "], x_pom[" +
                                   std::to_string(selectedBin->xpom_low) + "-" +
                                   std::to_string(selectedBin->xpom_high) + "]";
                        success = true;
                    } else {
                        std::cout << "This bin is already deleted." << std::endl;
                    }
                    break;
                }
            }
        }
    } else if(dimChoice == 2) {  // Merge along β
        if(dirChoice == 1) {  // Upward
            if(selectedBin->j < beta_bins.size() - 2) {
                double edgeToRemove = beta_bins[selectedBin->j + 1];
                mergeInfo = "Merging β bins [" + std::to_string(selectedBin->beta_low) + "-" +
                           std::to_string(selectedBin->beta_high) + "] and [" +
                           std::to_string(selectedBin->beta_high) + "-" +
                           std::to_string(beta_bins[selectedBin->j + 2]) + "]";
                beta_bins.erase(beta_bins.begin() + selectedBin->j + 1);
                success = true;
            } else {
                std::cout << "Cannot merge upward - already at highest β bin." << std::endl;
            }
        } else if(dirChoice == 2) {  // Downward
            if(selectedBin->j > 0) {
                double edgeToRemove = beta_bins[selectedBin->j];
                mergeInfo = "Merging β bins [" + std::to_string(beta_bins[selectedBin->j - 1]) + "-" +
                           std::to_string(selectedBin->beta_low) + "] and [" +
                           std::to_string(selectedBin->beta_low) + "-" +
                           std::to_string(selectedBin->beta_high) + "]";
                beta_bins.erase(beta_bins.begin() + selectedBin->j);
                success = true;
            } else {
                std::cout << "Cannot merge downward - already at lowest β bin." << std::endl;
            }
        } else {  // Delete - mark this specific 3D bin as deleted
            // Find and mark the specific 3D bin as deleted
            for(auto& bin : g_bins3D) {
                if(bin.i == selectedBin->i && bin.j == selectedBin->j && bin.k == selectedBin->k) {
                    if(!bin.deleted) {
                        bin.deleted = true;
                        mergeInfo = "Deleted 3D bin: Q²[" + std::to_string(selectedBin->Q2_low) + "-" +
                                   std::to_string(selectedBin->Q2_high) + "], β[" +
                                   std::to_string(selectedBin->beta_low) + "-" +
                                   std::to_string(selectedBin->beta_high) + "], x_pom[" +
                                   std::to_string(selectedBin->xpom_low) + "-" +
                                   std::to_string(selectedBin->xpom_high) + "]";
                        success = true;
                    } else {
                        std::cout << "This bin is already deleted." << std::endl;
                    }
                    break;
                }
            }
        }
    } else {  // Merge along x_pom
        if(dirChoice == 1) {  // Upward
            if(selectedBin->k < xpom_bins.size() - 2) {
                double edgeToRemove = xpom_bins[selectedBin->k + 1];
                mergeInfo = "Merging x_pom bins [" + std::to_string(selectedBin->xpom_low) + "-" +
                           std::to_string(selectedBin->xpom_high) + "] and [" +
                           std::to_string(selectedBin->xpom_high) + "-" +
                           std::to_string(xpom_bins[selectedBin->k + 2]) + "]";
                xpom_bins.erase(xpom_bins.begin() + selectedBin->k + 1);
                success = true;
            } else {
                std::cout << "Cannot merge upward - already at highest x_pom bin." << std::endl;
            }
        } else if(dirChoice == 2) {  // Downward
            if(selectedBin->k > 0) {
                double edgeToRemove = xpom_bins[selectedBin->k];
                mergeInfo = "Merging x_pom bins [" + std::to_string(xpom_bins[selectedBin->k - 1]) + "-" +
                           std::to_string(selectedBin->xpom_low) + "] and [" +
                           std::to_string(selectedBin->xpom_low) + "-" +
                           std::to_string(selectedBin->xpom_high) + "]";
                xpom_bins.erase(xpom_bins.begin() + selectedBin->k);
                success = true;
            } else {
                std::cout << "Cannot merge downward - already at lowest x_pom bin." << std::endl;
            }
        } else {  // Delete - mark this specific 3D bin as deleted
            // Find and mark the specific 3D bin as deleted
            for(auto& bin : g_bins3D) {
                if(bin.i == selectedBin->i && bin.j == selectedBin->j && bin.k == selectedBin->k) {
                    if(!bin.deleted) {
                        bin.deleted = true;
                        mergeInfo = "Deleted 3D bin: Q²[" + std::to_string(selectedBin->Q2_low) + "-" +
                                   std::to_string(selectedBin->Q2_high) + "], β[" +
                                   std::to_string(selectedBin->beta_low) + "-" +
                                   std::to_string(selectedBin->beta_high) + "], x_pom[" +
                                   std::to_string(selectedBin->xpom_low) + "-" +
                                   std::to_string(selectedBin->xpom_high) + "]";
                        success = true;
                    } else {
                        std::cout << "This bin is already deleted." << std::endl;
                    }
                    break;
                }
            }
        }
    }

    if(success) {
        std::cout << "\n✓ " << mergeInfo << std::endl;

        // If this was a merge operation (not delete), rebuild the 3D bins
        if(dirChoice != 3) {
            std::cout << "\nRebuilding 3D bin structure after merge..." << std::endl;
            RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
        }

        std::cout << "\nRegenerating plots with updated binning..." << std::endl;
        CreatePlots(inputFile, hist3d, Q2_bins, beta_bins, xpom_bins, y_bins, false);
        std::cout << "✓ Plots regenerated. Check figs/*_preview.png" << std::endl;

        // Show updated 3D statistics
        std::cout << "\nUpdated 3D Bin Statistics:" << std::endl;
        Show3DBinSummary(hist3d, Q2_bins, beta_bins, xpom_bins);

        // Loop back to show updated list
        std::cout << "\nContinuing to next bin (updated list will be shown)..." << std::endl;
    }
    }  // End of while loop
}

//=========================================================================
// 3D Binning Workshop - Interactive 3D binning with live feedback
//=========================================================================
void BinningWorkshop3D(TFile* inputFile, TH3D* h_d3sigma,
                       std::vector<double>& Q2_bins,
                       std::vector<double>& beta_bins,
                       std::vector<double>& xpom_bins,
                       const std::vector<double>& y_bins) {

    bool done = false;

    while(!done) {
        // Show compact 3D statistics summary
        Show3DBinSummary(h_d3sigma, Q2_bins, beta_bins, xpom_bins);

        std::cout << "\n3D Binning Workshop" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  1. Merge/Delete Q² bins" << std::endl;
        std::cout << "  2. Merge/Delete β bins" << std::endl;
        std::cout << "  3. Merge/Delete x_pom bins" << std::endl;
        std::cout << "  4. Merge individual 3D bins (surgical control)" << std::endl;
        std::cout << "  5. Show detailed 3D statistics table" << std::endl;
        std::cout << "  6. Return to main menu" << std::endl;

        int choice = GetIntInput("Enter choice (1-6): ");

        switch(choice) {
            case 1: {
                std::cout << "\n=== Adjusting Q² binning ===" << std::endl;
                std::vector<double> new_Q2_bins = MergeDeleteBins(Q2_bins, "Q²");
                if(new_Q2_bins.size() != Q2_bins.size()) {
                    Q2_bins = new_Q2_bins;
                    std::cout << "\n→ Rebuilding 3D bins with new Q² binning..." << std::endl;
                    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
                    std::cout << "\n→ Regenerating plots with new Q² binning..." << std::endl;
                    CreatePlots(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins, false);
                    std::cout << "✓ Plots updated! Check figs/ directory." << std::endl;
                }
                break;
            }
            case 2: {
                std::cout << "\n=== Adjusting β binning ===" << std::endl;
                std::vector<double> new_beta_bins = MergeDeleteBins(beta_bins, "β");
                if(new_beta_bins.size() != beta_bins.size()) {
                    beta_bins = new_beta_bins;
                    std::cout << "\n→ Rebuilding 3D bins with new β binning..." << std::endl;
                    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
                    std::cout << "\n→ Regenerating plots with new β binning..." << std::endl;
                    CreatePlots(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins, false);
                    std::cout << "✓ Plots updated! Check figs/ directory." << std::endl;
                }
                break;
            }
            case 3: {
                std::cout << "\n=== Adjusting x_pom binning ===" << std::endl;
                std::vector<double> new_xpom_bins = MergeDeleteBins(xpom_bins, "x_pom");
                if(new_xpom_bins.size() != xpom_bins.size()) {
                    xpom_bins = new_xpom_bins;
                    std::cout << "\n→ Rebuilding 3D bins with new x_pom binning..." << std::endl;
                    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
                    std::cout << "\n→ Regenerating plots with new x_pom binning..." << std::endl;
                    CreatePlots(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins, false);
                    std::cout << "✓ Plots updated! Check figs/ directory." << std::endl;
                }
                break;
            }
            case 4: {
                std::cout << "\n=== Merge Individual 3D Bins ===" << std::endl;
                MergeIndividual3DBin(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins);
                break;
            }
            case 5:
                Show3DBinStatistics(h_d3sigma, Q2_bins, beta_bins, xpom_bins);
                std::cout << "\nPress Enter to continue...";
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cin.get();
                break;
            case 6:
                done = true;
                std::cout << "\nReturning to main menu..." << std::endl;
                break;
            default:
                std::cout << "Invalid choice. Please try again." << std::endl;
                break;
        }
    }
}

//=========================================================================
// Create all plots with current binning
//=========================================================================
void CreatePlots(TFile* inputFile,
                TH3D* h_d3sigma,
                const std::vector<double>& Q2_bins,
                const std::vector<double>& beta_bins,
                const std::vector<double>& xpom_bins,
                const std::vector<double>& y_bins,
                bool saveFinal = false) {

    std::cout << "\n========================================" << std::endl;
    std::cout << "Generating plots with current binning..." << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  Q² bins: " << Q2_bins.size()-1 << std::endl;
    std::cout << "  β bins: " << beta_bins.size()-1 << std::endl;
    std::cout << "  x_pom bins: " << xpom_bins.size()-1 << std::endl;
    std::cout << "  y bins: " << y_bins.size()-1 << std::endl;

    std::string suffix = saveFinal ? "_final" : "_preview";

    //=========================================================================
    // Plot 1: Q² vs x_IP
    //=========================================================================
    TH2D* h_Q2_vs_xpom = (TH2D*)inputFile->Get("Q2_vs_xpom_MC");
    if(h_Q2_vs_xpom) {
        TCanvas* c1 = new TCanvas("c1", "Q2 vs xIP", 800, 700);
        c1->SetLogy();
        c1->SetLogz();
        c1->SetRightMargin(0.15);

        h_Q2_vs_xpom->SetTitle("Q^{2} vs x_{IP} (with event counts)");
        h_Q2_vs_xpom->GetXaxis()->SetTitle("log_{10}(x_{IP})");
        h_Q2_vs_xpom->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_Q2_vs_xpom->GetXaxis()->SetTitleSize(0.05);
        h_Q2_vs_xpom->GetYaxis()->SetTitleSize(0.05);
        h_Q2_vs_xpom->GetZaxis()->SetTitleSize(0.05);
        h_Q2_vs_xpom->GetXaxis()->SetLabelSize(0.045);
        h_Q2_vs_xpom->GetYaxis()->SetLabelSize(0.045);
        h_Q2_vs_xpom->Draw("COLZ");

        DrawBinningGridWithCounts3D(h_d3sigma, xpom_bins, Q2_bins, beta_bins, 1,
                                     true,
                                     h_Q2_vs_xpom->GetXaxis()->GetXmin(),
                                     h_Q2_vs_xpom->GetXaxis()->GetXmax(),
                                     h_Q2_vs_xpom->GetYaxis()->GetXmin(),
                                     h_Q2_vs_xpom->GetYaxis()->GetXmax(),
                                     false, "x_IP", "Q²");

        if(saveFinal) {
            c1->SaveAs("figs/binning_Q2_vs_xIP_counts.png");
            c1->SaveAs("figs/binning_Q2_vs_xIP_counts.pdf");
        } else {
            c1->SaveAs(("figs/binning_Q2_vs_xIP" + suffix + ".png").c_str());
        }
        delete c1;
    }

    //=========================================================================
    // Plot 2: Q² vs β
    //=========================================================================
    TH2D* h_Q2_vs_beta = (TH2D*)inputFile->Get("Q2_vs_beta_MC");
    if(h_Q2_vs_beta) {
        TCanvas* c2 = new TCanvas("c2", "Q2 vs beta", 800, 700);
        c2->SetLogy();
        c2->SetLogz();
        c2->SetRightMargin(0.15);

        h_Q2_vs_beta->SetTitle("Q^{2} vs #beta (with event counts)");
        h_Q2_vs_beta->GetXaxis()->SetTitle("#beta");
        h_Q2_vs_beta->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_Q2_vs_beta->GetXaxis()->SetTitleSize(0.05);
        h_Q2_vs_beta->GetYaxis()->SetTitleSize(0.05);
        h_Q2_vs_beta->GetZaxis()->SetTitleSize(0.05);
        h_Q2_vs_beta->GetXaxis()->SetLabelSize(0.045);
        h_Q2_vs_beta->GetYaxis()->SetLabelSize(0.045);
        h_Q2_vs_beta->Draw("COLZ");

        DrawBinningGridWithCounts3D(h_d3sigma, beta_bins, Q2_bins, xpom_bins, 0,
                                     false,
                                     h_Q2_vs_beta->GetXaxis()->GetXmin(),
                                     h_Q2_vs_beta->GetXaxis()->GetXmax(),
                                     h_Q2_vs_beta->GetYaxis()->GetXmin(),
                                     h_Q2_vs_beta->GetYaxis()->GetXmax(),
                                     false, "β", "Q²");

        if(saveFinal) {
            c2->SaveAs("figs/binning_Q2_vs_beta_counts.png");
            c2->SaveAs("figs/binning_Q2_vs_beta_counts.pdf");
        } else {
            c2->SaveAs(("figs/binning_Q2_vs_beta" + suffix + ".png").c_str());
        }
        delete c2;
    }

    //=========================================================================
    // Plot 3: β vs x_IP
    //=========================================================================
    TH2D* h_beta_vs_xpom = (TH2D*)inputFile->Get("beta_vs_xpom_MC");
    if(h_beta_vs_xpom) {
        TCanvas* c3 = new TCanvas("c3", "beta vs xIP", 800, 700);
        c3->SetLogz();
        c3->SetRightMargin(0.15);

        h_beta_vs_xpom->SetTitle("#beta vs x_{IP} (with event counts)");
        h_beta_vs_xpom->GetXaxis()->SetTitle("log_{10}(x_{IP})");
        h_beta_vs_xpom->GetYaxis()->SetTitle("#beta");
        h_beta_vs_xpom->GetXaxis()->SetTitleSize(0.05);
        h_beta_vs_xpom->GetYaxis()->SetTitleSize(0.05);
        h_beta_vs_xpom->GetZaxis()->SetTitleSize(0.05);
        h_beta_vs_xpom->GetXaxis()->SetLabelSize(0.045);
        h_beta_vs_xpom->GetYaxis()->SetLabelSize(0.045);
        h_beta_vs_xpom->Draw("COLZ");

        DrawBinningGridWithCounts3D(h_d3sigma, xpom_bins, beta_bins, Q2_bins, 2,
                                     true,
                                     h_beta_vs_xpom->GetXaxis()->GetXmin(),
                                     h_beta_vs_xpom->GetXaxis()->GetXmax(),
                                     h_beta_vs_xpom->GetYaxis()->GetXmin(),
                                     h_beta_vs_xpom->GetYaxis()->GetXmax(),
                                     false, "x_IP", "β");

        if(saveFinal) {
            c3->SaveAs("figs/binning_beta_vs_xIP_counts.png");
            c3->SaveAs("figs/binning_beta_vs_xIP_counts.pdf");
        } else {
            c3->SaveAs(("figs/binning_beta_vs_xIP" + suffix + ".png").c_str());
        }
        delete c3;
    }

    //=========================================================================
    // Plot 4: Q² vs y
    //=========================================================================
    TH2D* h_Q2_vs_y = (TH2D*)inputFile->Get("Q2_vs_y_MC");
    if(h_Q2_vs_y) {
        TCanvas* c4 = new TCanvas("c4", "Q2 vs y", 800, 700);
        c4->SetLogy();
        c4->SetLogz();
        c4->SetRightMargin(0.15);

        h_Q2_vs_y->SetTitle("Q^{2} vs y (with event counts)");
        h_Q2_vs_y->GetXaxis()->SetTitle("y");
        h_Q2_vs_y->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_Q2_vs_y->GetXaxis()->SetTitleSize(0.05);
        h_Q2_vs_y->GetYaxis()->SetTitleSize(0.05);
        h_Q2_vs_y->GetZaxis()->SetTitleSize(0.05);
        h_Q2_vs_y->GetXaxis()->SetLabelSize(0.045);
        h_Q2_vs_y->GetYaxis()->SetLabelSize(0.045);
        h_Q2_vs_y->Draw("COLZ");

        DrawBinningGridWithCounts(h_Q2_vs_y, y_bins, Q2_bins, false,
                                   h_Q2_vs_y->GetXaxis()->GetXmin(),
                                   h_Q2_vs_y->GetXaxis()->GetXmax(),
                                   h_Q2_vs_y->GetYaxis()->GetXmin(),
                                   h_Q2_vs_y->GetYaxis()->GetXmax(),
                                   false, "y", "Q²");

        if(saveFinal) {
            c4->SaveAs("figs/binning_Q2_vs_y_counts.png");
            c4->SaveAs("figs/binning_Q2_vs_y_counts.pdf");
        } else {
            c4->SaveAs(("figs/binning_Q2_vs_y" + suffix + ".png").c_str());
        }
        delete c4;
    }

    //=========================================================================
    // Combined 2x2 plot
    //=========================================================================
    TCanvas* c_combined = new TCanvas("c_combined", "Binning Scheme", 1600, 1400);
    c_combined->Divide(2, 2);

    c_combined->cd(1);
    gPad->SetLogy();
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_Q2_vs_xpom) {
        h_Q2_vs_xpom->Draw("COLZ");
        DrawBinningGridWithCounts3D(h_d3sigma, xpom_bins, Q2_bins, beta_bins, 1,
                                     true,
                                     h_Q2_vs_xpom->GetXaxis()->GetXmin(),
                                     h_Q2_vs_xpom->GetXaxis()->GetXmax(),
                                     h_Q2_vs_xpom->GetYaxis()->GetXmin(),
                                     h_Q2_vs_xpom->GetYaxis()->GetXmax(),
                                     false, "x_IP", "Q²");
    }

    c_combined->cd(2);
    gPad->SetLogy();
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_Q2_vs_beta) {
        h_Q2_vs_beta->Draw("COLZ");
        DrawBinningGridWithCounts3D(h_d3sigma, beta_bins, Q2_bins, xpom_bins, 0,
                                     false,
                                     h_Q2_vs_beta->GetXaxis()->GetXmin(),
                                     h_Q2_vs_beta->GetXaxis()->GetXmax(),
                                     h_Q2_vs_beta->GetYaxis()->GetXmin(),
                                     h_Q2_vs_beta->GetYaxis()->GetXmax(),
                                     false, "β", "Q²");
    }

    c_combined->cd(3);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_beta_vs_xpom) {
        h_beta_vs_xpom->Draw("COLZ");
        DrawBinningGridWithCounts3D(h_d3sigma, xpom_bins, beta_bins, Q2_bins, 2,
                                     true,
                                     h_beta_vs_xpom->GetXaxis()->GetXmin(),
                                     h_beta_vs_xpom->GetXaxis()->GetXmax(),
                                     h_beta_vs_xpom->GetYaxis()->GetXmin(),
                                     h_beta_vs_xpom->GetYaxis()->GetXmax(),
                                     false, "x_IP", "β");
    }

    c_combined->cd(4);
    gPad->SetLogy();
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_Q2_vs_y) {
        h_Q2_vs_y->Draw("COLZ");
        DrawBinningGridWithCounts(h_Q2_vs_y, y_bins, Q2_bins, false,
                                   h_Q2_vs_y->GetXaxis()->GetXmin(),
                                   h_Q2_vs_y->GetXaxis()->GetXmax(),
                                   h_Q2_vs_y->GetYaxis()->GetXmin(),
                                   h_Q2_vs_y->GetYaxis()->GetXmax(),
                                   false, "y", "Q²");
    }

    if(saveFinal) {
        c_combined->SaveAs("figs/binning_scheme_combined_counts.png");
        c_combined->SaveAs("figs/binning_scheme_combined_counts.pdf");
    } else {
        c_combined->SaveAs(("figs/binning_scheme_combined" + suffix + ".png").c_str());
    }
    delete c_combined;

    std::cout << "Plots created in figs/ directory" << std::endl;
}

//=========================================================================
// Main function
//=========================================================================
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    const char* inputFileName = argv[1];

    // Open input file
    TFile* inputFile = TFile::Open(inputFileName);
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFileName << std::endl;
        return 1;
    }

    // Set global style
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(256);
    gStyle->SetPalette(kBird);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.3, "Y");

    // Create output directory
    gSystem->Exec("mkdir -p figs");

    std::cout << "\n========================================" << std::endl;
    std::cout << "Interactive Binning Scheme Plotter" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Input file: " << inputFileName << std::endl;

    //=========================================================================
    // Extract initial binning from 3D histogram
    //=========================================================================
    TH3D* h_d3sigma = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_MC");
    if(!h_d3sigma) {
        std::cerr << "Error: Could not find d3sigma_dQ2dbeta_dxpom_MC histogram" << std::endl;
        inputFile->Close();
        return 1;
    }

    // Extract initial binning
    std::vector<double> Q2_bins, beta_bins, xpom_bins, y_bins;

    TAxis* xaxis = h_d3sigma->GetXaxis();
    for(int i = 1; i <= xaxis->GetNbins() + 1; i++) {
        Q2_bins.push_back(xaxis->GetBinLowEdge(i));
    }

    TAxis* yaxis = h_d3sigma->GetYaxis();
    for(int i = 1; i <= yaxis->GetNbins() + 1; i++) {
        beta_bins.push_back(yaxis->GetBinLowEdge(i));
    }

    TAxis* zaxis = h_d3sigma->GetZaxis();
    for(int i = 1; i <= zaxis->GetNbins() + 1; i++) {
        xpom_bins.push_back(zaxis->GetBinLowEdge(i));
    }

    TH1D* h_y = (TH1D*)inputFile->Get("y_truth");
    if(h_y) {
        TAxis* y_ax = h_y->GetXaxis();
        for(int i = 1; i <= y_ax->GetNbins() + 1; i++) {
            y_bins.push_back(y_ax->GetBinLowEdge(i));
        }
    } else {
        for(int i = 0; i <= 10; i++) {
            y_bins.push_back(i * 0.1);
        }
    }

    std::cout << "\nInitial binning extracted from histograms:" << std::endl;
    DisplayBinEdges(Q2_bins, "Q²");
    DisplayBinEdges(beta_bins, "β");
    DisplayBinEdges(xpom_bins, "x_pom");
    DisplayBinEdges(y_bins, "y");

    // Initialize 3D bin structure
    std::cout << "\nInitializing 3D bin structure..." << std::endl;
    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);

    // Get required histograms for binning choices
    TH2D* h_Q2_vs_xpom = (TH2D*)inputFile->Get("Q2_vs_xpom_MC");
    TH2D* h_Q2_vs_beta = (TH2D*)inputFile->Get("Q2_vs_beta_MC");
    TH2D* h_beta_vs_xpom = (TH2D*)inputFile->Get("beta_vs_xpom_MC");
    TH2D* h_Q2_vs_y = (TH2D*)inputFile->Get("Q2_vs_y_MC");

    //=========================================================================
    // Interactive loop
    //=========================================================================
    bool satisfied = false;

    while(!satisfied) {
        // Generate preview plots
        CreatePlots(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins, false);

        std::cout << "\n========================================" << std::endl;
        std::cout << "Main Menu" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Please check the preview plots in figs/ directory" << std::endl;
        std::cout << "\nOptions:" << std::endl;
        std::cout << "  1. Adjust Q² binning" << std::endl;
        std::cout << "  2. Adjust β binning" << std::endl;
        std::cout << "  3. Adjust x_pom binning" << std::endl;
        std::cout << "  4. Adjust y binning" << std::endl;
        std::cout << "  5. Display all current bin edges" << std::endl;
        std::cout << "  6. Show bin statistics" << std::endl;
        std::cout << "  7. Accept current binning and save final plots" << std::endl;
        std::cout << "  8. Exit without saving" << std::endl;
        std::cout << "  9. 3D Binning Workshop (adjust Q², β, x_pom with live feedback)" << std::endl;

        int mainChoice = GetIntInput("Enter choice (1-9): ");

        switch(mainChoice) {
            case 1:
                if(h_Q2_vs_beta) {
                    Q2_bins = GetBinningChoice("Q²", Q2_bins, h_Q2_vs_beta, false);
                    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
                }
                break;
            case 2:
                if(h_Q2_vs_beta) {
                    beta_bins = GetBinningChoice("β", beta_bins, h_Q2_vs_beta, true);
                    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
                }
                break;
            case 3:
                if(h_Q2_vs_xpom) {
                    xpom_bins = GetBinningChoice("x_pom", xpom_bins, h_Q2_vs_xpom, true);
                    RebuildBins3D(Q2_bins, beta_bins, xpom_bins);
                }
                break;
            case 4:
                if(h_Q2_vs_y) {
                    y_bins = GetBinningChoice("y", y_bins, h_Q2_vs_y, true);
                }
                break;
            case 5:
                DisplayAllBinEdges(Q2_bins, beta_bins, xpom_bins, y_bins);
                break;
            case 6: {
                std::cout << "\nWhich statistics to show?" << std::endl;
                std::cout << "  1. 3D bins (Q², β, x_pom) - FULL 3D STATISTICS" << std::endl;
                std::cout << "  2. Q² vs x_IP (2D projection)" << std::endl;
                std::cout << "  3. Q² vs β (2D projection)" << std::endl;
                std::cout << "  4. β vs x_IP (2D projection)" << std::endl;
                std::cout << "  5. Q² vs y (2D projection)" << std::endl;
                int statsChoice = GetIntInput("Choice: ");
                switch(statsChoice) {
                    case 1:
                        Show3DBinStatistics(h_d3sigma, Q2_bins, beta_bins, xpom_bins);
                        break;
                    case 2:
                        if(h_Q2_vs_xpom) ShowBinStatistics(h_Q2_vs_xpom, xpom_bins, Q2_bins, true, "x_IP", "Q²");
                        break;
                    case 3:
                        if(h_Q2_vs_beta) ShowBinStatistics(h_Q2_vs_beta, beta_bins, Q2_bins, false, "β", "Q²");
                        break;
                    case 4:
                        if(h_beta_vs_xpom) ShowBinStatistics(h_beta_vs_xpom, xpom_bins, beta_bins, true, "x_IP", "β");
                        break;
                    case 5:
                        if(h_Q2_vs_y) ShowBinStatistics(h_Q2_vs_y, y_bins, Q2_bins, false, "y", "Q²");
                        break;
                }
                std::cout << "\nPress Enter to continue...";
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cin.get();
                break;
            }
            case 7:
                satisfied = true;
                std::cout << "\nSaving final plots..." << std::endl;
                CreatePlots(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins, true);

                // Ask if user wants to export bin edges
                std::cout << "\nWould you like to export bin edges to a file? (y/n): ";
                char exportChoice;
                std::cin >> exportChoice;
                if(exportChoice == 'y' || exportChoice == 'Y') {
                    ExportBinEdges("figs/bin_edges.csv", h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins);
                }

                // Ask if user wants to export modification history
                if(!g_modificationHistory.empty()) {
                    std::cout << "\nWould you like to export modification history? (y/n): ";
                    char historyChoice;
                    std::cin >> historyChoice;
                    if(historyChoice == 'y' || historyChoice == 'Y') {
                        ExportModificationHistory("figs/bin_modification_history.csv");
                    }
                } else {
                    std::cout << "\nNo bin modifications were made during this session." << std::endl;
                }

                std::cout << "\n========================================" << std::endl;
                std::cout << "All final plots saved successfully!" << std::endl;
                std::cout << "========================================" << std::endl;
                break;
            case 8:
                std::cout << "Exiting without saving final plots." << std::endl;
                inputFile->Close();
                return 0;
            case 9:
                BinningWorkshop3D(inputFile, h_d3sigma, Q2_bins, beta_bins, xpom_bins, y_bins);
                break;
            default:
                std::cout << "Invalid choice. Please try again." << std::endl;
        }
    }

    inputFile->Close();
    delete inputFile;

    return 0;
}
