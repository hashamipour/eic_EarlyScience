//g++ -std=c++17 -o Plot_BinningScheme_WithCounts Plot_BinningScheme_WithCounts.cpp `root-config --cflags --libs`

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
#include <iostream>
#include <vector>
#include <TMath.h>

void DrawBinningGridWithCounts(TH2D* hist, std::vector<double> xbins, std::vector<double> ybins,
                                bool isLogX = false, double xmin = 0, double xmax = 1,
                                double ymin = 0, double ymax = 1) {

    // Set line properties
    int lineColor = kRed;       // Changed from kBlack to kRed
    int lineWidth = 3;          // Increased from 2 to 3
    int lineStyle = 1; // solid

    // Draw vertical lines for x bins
    for(size_t i = 0; i < xbins.size(); i++) {
        double xval = xbins[i];
        if(isLogX) {
            xval = TMath::Log10(xval);
        }

        // Only draw if within range
        if(xval >= xmin && xval <= xmax) {
            TLine* line = new TLine(xval, ymin, xval, ymax);
            line->SetLineColor(lineColor);
            line->SetLineWidth(lineWidth);
            line->SetLineStyle(lineStyle);
            line->Draw("same");
        }
    }

    // Draw horizontal lines for y bins
    for(size_t i = 0; i < ybins.size(); i++) {
        double yval = ybins[i];

        // Only draw if within range
        if(yval >= ymin && yval <= ymax) {
            TLine* line = new TLine(xmin, yval, xmax, yval);
            line->SetLineColor(lineColor);
            line->SetLineWidth(lineWidth);
            line->SetLineStyle(lineStyle);
            line->Draw("same");
        }
    }

    // Now add event counts in each bin
    TLatex* latex = new TLatex();
    latex->SetTextColor(kBlack);    // Changed from kRed to kBlack
    latex->SetTextSize(0.020);      // Decreased from 0.025 to 0.020
    latex->SetTextAlign(22); // Center alignment
    latex->SetTextFont(42);

    // Loop over bins and calculate counts
    for(size_t i = 0; i < xbins.size() - 1; i++) {
        for(size_t j = 0; j < ybins.size() - 1; j++) {
            double xlow = xbins[i];
            double xhigh = xbins[i+1];
            double ylow = ybins[j];
            double yhigh = ybins[j+1];

            // Convert to log if needed for x
            if(isLogX) {
                xlow = TMath::Log10(xbins[i]);
                xhigh = TMath::Log10(xbins[i+1]);
            }

            // Calculate center of bin
            double xcenter = (xlow + xhigh) / 2.0;
            double ycenter = (ylow + yhigh) / 2.0;

            // Count events in this bin
            int count = 0;

            // Find histogram bins that correspond to this range
            int binxlow = hist->GetXaxis()->FindBin(xlow);
            int binxhigh = hist->GetXaxis()->FindBin(xhigh);
            int binylow = hist->GetYaxis()->FindBin(ylow);
            int binyhigh = hist->GetYaxis()->FindBin(yhigh);

            for(int ix = binxlow; ix <= binxhigh; ix++) {
                for(int iy = binylow; iy <= binyhigh; iy++) {
                    count += hist->GetBinContent(ix, iy);
                }
            }

            // Draw count if non-zero
            if(count > 0) {
                TString countStr = Form("%d", count);
                latex->DrawLatex(xcenter, ycenter, countStr);
            }
        }
    }
}

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

    //=========================================================================
    // Extract binning schemes from the 3D histogram (the actual analysis binning)
    //=========================================================================
    TH3D* h_d3sigma = (TH3D*)inputFile->Get("d3sigma_dQ2dbeta_dxpom_MC");
    if(!h_d3sigma) {
        std::cerr << "Error: Could not find d3sigma_dQ2dbeta_dxpom_MC histogram" << std::endl;
        std::cerr << "This histogram is needed to extract the binning scheme" << std::endl;
        return 1;
    }

    // Extract Q2 bins from X axis
    std::vector<double> Q2_bins;
    TAxis* xaxis = h_d3sigma->GetXaxis();
    for(int i = 1; i <= xaxis->GetNbins() + 1; i++) {
        Q2_bins.push_back(xaxis->GetBinLowEdge(i));
    }

    // Extract beta bins from Y axis
    std::vector<double> beta_bins;
    TAxis* yaxis = h_d3sigma->GetYaxis();
    for(int i = 1; i <= yaxis->GetNbins() + 1; i++) {
        beta_bins.push_back(yaxis->GetBinLowEdge(i));
    }

    // Extract x_pom bins from Z axis
    std::vector<double> xpom_bins;
    TAxis* zaxis = h_d3sigma->GetZaxis();
    for(int i = 1; i <= zaxis->GetNbins() + 1; i++) {
        xpom_bins.push_back(zaxis->GetBinLowEdge(i));
    }

    // Extract y bins from y histogram
    std::vector<double> y_bins;
    TH1D* h_y = (TH1D*)inputFile->Get("y_truth");
    if(h_y) {
        TAxis* y_ax = h_y->GetXaxis();
        for(int i = 1; i <= y_ax->GetNbins() + 1; i++) {
            y_bins.push_back(y_ax->GetBinLowEdge(i));
        }
    } else {
        // Fallback if histogram not found
        std::cerr << "Warning: y_truth histogram not found, using default y binning" << std::endl;
        for(int i = 0; i <= 10; i++) {
            y_bins.push_back(i * 0.1);
        }
    }

    std::cout << "Extracted binning from histograms:" << std::endl;
    std::cout << "  Q² bins: " << Q2_bins.size()-1 << " bins" << std::endl;
    std::cout << "  β bins: " << beta_bins.size()-1 << " bins" << std::endl;
    std::cout << "  x_pom bins: " << xpom_bins.size()-1 << " bins" << std::endl;
    std::cout << "  y bins: " << y_bins.size()-1 << " bins" << std::endl;

    //=========================================================================
    // Create output directory if it doesn't exist
    //=========================================================================
    gSystem->Exec("mkdir -p figs");

    //=========================================================================
    // Plot 1: Q² vs x_IP with counts
    //=========================================================================
    TH2D* h_Q2_vs_xpom = (TH2D*)inputFile->Get("Q2_vs_xpom_MC");
    if(h_Q2_vs_xpom) {
        TCanvas* c1 = new TCanvas("c1", "Q2 vs xIP with counts", 800, 700);
        c1->SetLogy();  // Logarithmic Y axis (Q²)
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

        // Overlay binning grid with counts
        DrawBinningGridWithCounts(h_Q2_vs_xpom, xpom_bins, Q2_bins, true,
                                   h_Q2_vs_xpom->GetXaxis()->GetXmin(),
                                   h_Q2_vs_xpom->GetXaxis()->GetXmax(),
                                   h_Q2_vs_xpom->GetYaxis()->GetXmin(),
                                   h_Q2_vs_xpom->GetYaxis()->GetXmax());

        c1->SaveAs("figs/binning_Q2_vs_xIP_counts.png");
        c1->SaveAs("figs/binning_Q2_vs_xIP_counts.pdf");
        delete c1;
        std::cout << "Created: figs/binning_Q2_vs_xIP_counts.png" << std::endl;
    } else {
        std::cerr << "Warning: Q2_vs_xpom_MC histogram not found" << std::endl;
    }

    //=========================================================================
    // Plot 2: Q² vs β with counts
    //=========================================================================
    TH2D* h_Q2_vs_beta = (TH2D*)inputFile->Get("Q2_vs_beta_MC");
    if(h_Q2_vs_beta) {
        TCanvas* c2 = new TCanvas("c2", "Q2 vs beta with counts", 800, 700);
        c2->SetLogy();  // Logarithmic Y axis (Q²)
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

        // Overlay binning grid with counts
        DrawBinningGridWithCounts(h_Q2_vs_beta, beta_bins, Q2_bins, false,
                                   h_Q2_vs_beta->GetXaxis()->GetXmin(),
                                   h_Q2_vs_beta->GetXaxis()->GetXmax(),
                                   h_Q2_vs_beta->GetYaxis()->GetXmin(),
                                   h_Q2_vs_beta->GetYaxis()->GetXmax());

        c2->SaveAs("figs/binning_Q2_vs_beta_counts.png");
        c2->SaveAs("figs/binning_Q2_vs_beta_counts.pdf");
        delete c2;
        std::cout << "Created: figs/binning_Q2_vs_beta_counts.png" << std::endl;
    } else {
        std::cerr << "Warning: Q2_vs_beta_MC histogram not found" << std::endl;
    }

    //=========================================================================
    // Plot 3: β vs x_IP with counts
    //=========================================================================
    TH2D* h_beta_vs_xpom = (TH2D*)inputFile->Get("beta_vs_xpom_MC");
    if(h_beta_vs_xpom) {
        TCanvas* c3 = new TCanvas("c3", "beta vs xIP with counts", 800, 700);
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

        // Overlay binning grid with counts
        DrawBinningGridWithCounts(h_beta_vs_xpom, xpom_bins, beta_bins, true,
                                   h_beta_vs_xpom->GetXaxis()->GetXmin(),
                                   h_beta_vs_xpom->GetXaxis()->GetXmax(),
                                   h_beta_vs_xpom->GetYaxis()->GetXmin(),
                                   h_beta_vs_xpom->GetYaxis()->GetXmax());

        c3->SaveAs("figs/binning_beta_vs_xIP_counts.png");
        c3->SaveAs("figs/binning_beta_vs_xIP_counts.pdf");
        delete c3;
        std::cout << "Created: figs/binning_beta_vs_xIP_counts.png" << std::endl;
    } else {
        std::cerr << "Warning: beta_vs_xpom_MC histogram not found" << std::endl;
    }

    //=========================================================================
    // Plot 4: Q² vs y with counts
    //=========================================================================
    TH2D* h_Q2_vs_y = (TH2D*)inputFile->Get("Q2_vs_y_MC");
    if(h_Q2_vs_y) {
        TCanvas* c4 = new TCanvas("c4", "Q2 vs y with counts", 800, 700);
        c4->SetLogy();  // Logarithmic Y axis (Q²)
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

        // Overlay binning grid with counts
        DrawBinningGridWithCounts(h_Q2_vs_y, y_bins, Q2_bins, false,
                                   h_Q2_vs_y->GetXaxis()->GetXmin(),
                                   h_Q2_vs_y->GetXaxis()->GetXmax(),
                                   h_Q2_vs_y->GetYaxis()->GetXmin(),
                                   h_Q2_vs_y->GetYaxis()->GetXmax());

        c4->SaveAs("figs/binning_Q2_vs_y_counts.png");
        c4->SaveAs("figs/binning_Q2_vs_y_counts.pdf");
        delete c4;
        std::cout << "Created: figs/binning_Q2_vs_y_counts.png" << std::endl;
    } else {
        std::cerr << "Warning: Q2_vs_y_MC histogram not found" << std::endl;
    }

    //=========================================================================
    // Create combined 2x2 plot with counts
    //=========================================================================
    TCanvas* c_combined = new TCanvas("c_combined", "Binning Scheme with Counts", 1600, 1400);
    c_combined->Divide(2, 2);

    // Top left: Q² vs x_IP (Q² on Y axis - logarithmic)
    c_combined->cd(1);
    gPad->SetLogy();  // Logarithmic Y axis (Q²)
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_Q2_vs_xpom) {
        h_Q2_vs_xpom->Draw("COLZ");
        DrawBinningGridWithCounts(h_Q2_vs_xpom, xpom_bins, Q2_bins, true,
                                   h_Q2_vs_xpom->GetXaxis()->GetXmin(),
                                   h_Q2_vs_xpom->GetXaxis()->GetXmax(),
                                   h_Q2_vs_xpom->GetYaxis()->GetXmin(),
                                   h_Q2_vs_xpom->GetYaxis()->GetXmax());
    }

    // Top right: Q² vs β (Q² on Y axis - logarithmic)
    c_combined->cd(2);
    gPad->SetLogy();  // Logarithmic Y axis (Q²)
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_Q2_vs_beta) {
        h_Q2_vs_beta->Draw("COLZ");
        DrawBinningGridWithCounts(h_Q2_vs_beta, beta_bins, Q2_bins, false,
                                   h_Q2_vs_beta->GetXaxis()->GetXmin(),
                                   h_Q2_vs_beta->GetXaxis()->GetXmax(),
                                   h_Q2_vs_beta->GetYaxis()->GetXmin(),
                                   h_Q2_vs_beta->GetYaxis()->GetXmax());
    }

    // Bottom left: β vs x_IP (β on Y axis - linear)
    c_combined->cd(3);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_beta_vs_xpom) {
        h_beta_vs_xpom->Draw("COLZ");
        DrawBinningGridWithCounts(h_beta_vs_xpom, xpom_bins, beta_bins, true,
                                   h_beta_vs_xpom->GetXaxis()->GetXmin(),
                                   h_beta_vs_xpom->GetXaxis()->GetXmax(),
                                   h_beta_vs_xpom->GetYaxis()->GetXmin(),
                                   h_beta_vs_xpom->GetYaxis()->GetXmax());
    }

    // Bottom right: Q² vs y (Q² on Y axis - logarithmic)
    c_combined->cd(4);
    gPad->SetLogy();  // Logarithmic Y axis (Q²)
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(h_Q2_vs_y) {
        h_Q2_vs_y->Draw("COLZ");
        DrawBinningGridWithCounts(h_Q2_vs_y, y_bins, Q2_bins, false,
                                   h_Q2_vs_y->GetXaxis()->GetXmin(),
                                   h_Q2_vs_y->GetXaxis()->GetXmax(),
                                   h_Q2_vs_y->GetYaxis()->GetXmin(),
                                   h_Q2_vs_y->GetYaxis()->GetXmax());
    }

    c_combined->SaveAs("figs/binning_scheme_combined_counts.png");
    c_combined->SaveAs("figs/binning_scheme_combined_counts.pdf");
    delete c_combined;
    std::cout << "Created: figs/binning_scheme_combined_counts.png" << std::endl;

    inputFile->Close();
    delete inputFile;

    std::cout << "\nAll binning scheme plots with event counts created successfully!" << std::endl;

    return 0;
}
