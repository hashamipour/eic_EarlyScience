#include "Utility.hpp"
#include <vector>
#include <map>
#include <TStyle.h>
#include <TColor.h>
#include <TMath.h>
#include <iostream>

void SetCustomPalette(const std::string& paletteName) {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    UInt_t colors[NRGBs];

    // Define palette map
    std::map<std::string, std::vector<UInt_t>> paletteMap = {
        { "PetalFlare", {
            0xE85A78, 0xF57C3A, 0xFFE94F, 0xA4D65E, 0x2FAE5E
        }},
        { "CitrusFade", {
            0xFFFF00, 0xC6F94D, 0x7ED957, 0x4CCB5B, 0x249E57
        }},
        { "LemonLush", {
            0xFEFE69, 0xDDF969, 0xA9F36A, 0x78EC6C, 0x57E86B
        }},
        { "CandyEclipse", {
            0xFEE27A, 0xFEA959, 0xFF605D, 0xF13484, 0x8D379E
        }},
        { "SolarSplash", {
            0x35BBCA, 0x0191B4, 0xF8D90F, 0xD3DD18, 0xFE7A15
        }},
        { "SolarBloom", {
            0xE85A78, // Rosy Crimson
            0xF57C3A, // Early Orange
            0xFFE94F, // Bright Yellow
            0xA4D65E, // Zesty Green
            0x2FAE5E  // Slightly Darker Green
        }}
    };

    // Fallback to PetalFlare if name not found
    auto it = paletteMap.find(paletteName);
    if (it == paletteMap.end()) {
        it = paletteMap.find("PetalFlare");
    }

    for (int i = 0; i < NRGBs; ++i) {
        colors[i] = it->second[i];
    }

    // Convert to RGB
    Double_t stops[NRGBs], red[NRGBs], green[NRGBs], blue[NRGBs];
    for (int i = 0; i < NRGBs; ++i) {
        stops[i] = (Double_t)i / (NRGBs - 1.0);
        UChar_t r = (colors[i] >> 16) & 0xFF;
        UChar_t g = (colors[i] >> 8) & 0xFF;
        UChar_t b = colors[i] & 0xFF;
        red[i] = (Double_t)r / 255.0;
        green[i] = (Double_t)g / 255.0;
        blue[i] = (Double_t)b / 255.0;
    }

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void SetCustomPalette(const int& paletteID = 0) {
    std::map<int, std::string> idToName = {
        {0, "PetalFlare"},
        {1, "CitrusFade"},
        {2, "LemonLush"},
        {3, "CandyEclipse"},
        {4, "SolarSplash"},
        {5, "SolarBloom"}
    };
    SetCustomPalette(idToName.count(paletteID) ? idToName[paletteID] : "PetalFlare");
    // Default to PetalFlare if ID not found
}



/**
 * @brief Generates logarithmically spaced bin edges with custom rounding using TMath::Nint.
 *
 * This function calculates bin edges on a logarithmic scale, applies custom rounding,
 * and then processes the resulting vector to ensure all bin edges are unique and monotonic.
 * for now it's only intended for Q2 binning.
 *
 * @param min_val The overall minimum value for the histogram range. Must be > 0.
 * @param max_val The overall maximum value for the histogram range.
 * @param n_bins The desired number of bins.
 * @return A vector of doubles containing the correctly ordered and rounded bin edges.
 */
std::vector<Double_t> GetRoundedLogBins(double min_val, double max_val, int n_bins) {

    if (min_val <= 0) {
        Logger::error("The minimum value for logarithmic binning must be greater than zero.");
        return std::vector<Double_t>();
    }

    std::vector<Double_t> bin_edges;
    bin_edges.reserve(n_bins + 1);

    double log_min = TMath::Log10(min_val);
    double log_max = TMath::Log10(max_val);
    double log_step = (log_max - log_min) / n_bins;

    // Step 1: Generate all rounded bin edges without checking for uniqueness
    for (int i = 0; i <= n_bins; ++i) {
        double current_log_val = log_min + i * log_step;
        double unrounded_edge = TMath::Power(10, current_log_val);
        double rounded_edge;

        if (unrounded_edge <= 10.0) {
            rounded_edge = TMath::Nint(unrounded_edge * 10.0) / 10.0;
        } else if (unrounded_edge > 10.0 && unrounded_edge <= 20.0) {
            rounded_edge = TMath::Nint(unrounded_edge * 2.0) / 2.0;
        } else if (unrounded_edge > 20.0 && unrounded_edge <= 40.0) {
            rounded_edge = TMath::Nint(unrounded_edge);
        } else if (unrounded_edge > 40.0 && unrounded_edge <= 100.0) {
            rounded_edge = TMath::Nint(unrounded_edge / 5) * 5.0;
        } else if (unrounded_edge > 100.0 && unrounded_edge <= 150.0) {
            rounded_edge = TMath::Nint(unrounded_edge / 10.0) * 10.0;
        } else {
            rounded_edge = TMath::Nint(unrounded_edge / 50.0) * 50.0;
        }
        
        bin_edges.push_back(rounded_edge);
    }
    
    // Step 2: Remove duplicates to ensure monotonicity
    // TMath::Sort might be an option, but std::sort is more standard.
    std::sort(bin_edges.begin(), bin_edges.end());
    
    // std::unique moves all unique elements to the front and returns an iterator to the new end.
    auto last = std::unique(bin_edges.begin(), bin_edges.end());
    
    // Erase the elements from the new end to the original end.
    bin_edges.erase(last, bin_edges.end());
    
    return bin_edges;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Manually defines and returns a vector of bin edges for Q2 binning.
 *
 * This function provides a predefined vector of Q2 bin edges, which can be
 * easily modified or replaced as needed. The returned vector is ready to be
 * used for histogram creation or re-binning.
 *
 * @return A std::vector<Double_t> containing the custom bin edges.
 */
std::vector<Double_t> GetManualQ2Bins() {
    
    // Define bin edges here.
    std::vector<Double_t> q2_bins = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
         21, 22.5, 23, 24, 25, 26, 27.5, 29, 31, 33, 40, 50, 60, 70, 85, 100, 120, 140, 170, 200};
    
    // Note: ensure the vector is sorted and unique
    
    return q2_bins;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Generates logarithmically spaced bin edges.
 *
 * @param min_val The overall minimum value for the histogram range. Must be > 0.
 * @param max_val The overall maximum value for the histogram range.
 * @param n_bins The desired number of bins.
 * @return A vector of doubles containing the correctly ordered bin edges.
 */
std::vector<Double_t> GetLogBins(double min_val, double max_val, int n_bins) {

    // Ensure that min_val is positive and non-zero for logarithmic binning.
    if (min_val <= 0) {
        Logger::error("The minimum value for logarithmic binning must be greater than zero.");
        return std::vector<Double_t>();
    }

    std::vector<Double_t> bin_edges;
    bin_edges.reserve(n_bins + 1);

    // Calculate the logarithmic step
    double log_min = TMath::Log10(min_val);
    double log_max = TMath::Log10(max_val);
    double log_step = (log_max - log_min) / n_bins;

    // Populate the vector with bin edges
    for (int i = 0; i <= n_bins; ++i) {
        double current_log_val = log_min + i * log_step;
        bin_edges.push_back(TMath::Power(10, current_log_val));
    }
    
    return bin_edges;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Generates linearly spaced bin edges.
 *
 * @param min_val The overall minimum value for the histogram range.
 * @param max_val The overall maximum value for the histogram range.
 * @param n_bins The desired number of bins.
 * @return A vector of doubles containing the correctly ordered bin edges.
 */
std::vector<Double_t> GetLinBins(double min_val, double max_val, int n_bins) {

    std::vector<Double_t> bin_edges;
    bin_edges.reserve(n_bins + 1);

    // Calculate the linear step
    double lin_step = (max_val - min_val) / n_bins;

    // Populate the vector with bin edges
    for (int i = 0; i <= n_bins; ++i) {
        double current_lin_val = min_val + i * lin_step;
        bin_edges.push_back(current_lin_val);
    }
    
    return bin_edges;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////