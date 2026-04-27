#pragma once

#include <string>
#include <vector>

class TH1D;

// Fixed 3D (Q^2, x_pom, beta) binning used to map events to a global
// "k index" for relative-resolution studies in the DDIS skim.
//
// The edge vectors are populated at startup by the skim's main() from the
// YAML bin file; `kRelResN*` are cached sizes for O(1) bin indexing.

extern std::vector<double> kRelResQ2Bins;
extern std::vector<double> kRelResBetaBins;
extern std::vector<double> kRelResXpomBins;
extern int kRelResNQ2;
extern int kRelResNBeta;
extern int kRelResNXpom;
extern int kRelResNBins;

int FindBinIndex(const std::vector<double>& edges, double value);
int GetRelResGlobalBin(double q2, double xpom, double beta);

// Running moments over the global-k space. sum/sumsq/count are parallel
// arrays indexed by k; RMS / RMSError are computed on demand.
struct ResAccum {
    std::vector<double> sum;
    std::vector<double> sumsq;
    std::vector<int>    count;

    explicit ResAccum(int nBins = 0)
        : sum(nBins, 0.0), sumsq(nBins, 0.0), count(nBins, 0) {}

    void   Fill(int k, double value);
    double RMS(int k) const;
    double RMSError(int k) const;
};

TH1D* BuildRelResVsKHist(const std::string& name,
                         const std::string& title,
                         const ResAccum& acc);
