#include "ResolutionBinning.hpp"

#include <TH1D.h>

#include <cmath>

std::vector<double> kRelResQ2Bins;
std::vector<double> kRelResBetaBins;
std::vector<double> kRelResXpomBins;
int kRelResNQ2 = 0;
int kRelResNBeta = 0;
int kRelResNXpom = 0;
int kRelResNBins = 0;

int FindBinIndex(const std::vector<double>& edges, double value) {
    if (edges.size() < 2) return -1;
    for (size_t i = 0; i + 1 < edges.size(); i++) {
        const double lo = edges[i];
        const double hi = edges[i + 1];
        const bool last = (i + 1 == edges.size() - 1);
        if (value >= lo && (value < hi || (last && value <= hi))) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

int GetRelResGlobalBin(double q2, double xpom, double beta) {
    const int iQ2 = FindBinIndex(kRelResQ2Bins, q2);
    const int iXpom = FindBinIndex(kRelResXpomBins, xpom);
    const int iBeta = FindBinIndex(kRelResBetaBins, beta);
    if (iQ2 < 0 || iXpom < 0 || iBeta < 0) return -1;
    return iBeta + iXpom * kRelResNBeta + iQ2 * (kRelResNBeta * kRelResNXpom);
}

void ResAccum::Fill(int k, double value) {
    if (k < 0 || k >= static_cast<int>(sum.size())) return;
    sum[k]   += value;
    sumsq[k] += value * value;
    count[k] += 1;
}

double ResAccum::RMS(int k) const {
    if (k < 0 || k >= static_cast<int>(sum.size()) || count[k] == 0) return 0.0;
    const double mean = sum[k] / static_cast<double>(count[k]);
    const double var  = sumsq[k] / static_cast<double>(count[k]) - mean * mean;
    return (var > 0.0) ? std::sqrt(var) : 0.0;
}

double ResAccum::RMSError(int k) const {
    if (k < 0 || k >= static_cast<int>(sum.size()) || count[k] < 2) return 0.0;
    const double rms = RMS(k);
    return rms / std::sqrt(2.0 * (count[k] - 1.0));
}

TH1D* BuildRelResVsKHist(const std::string& name,
                         const std::string& title,
                         const ResAccum& acc) {
    TH1D* h = new TH1D(name.c_str(), title.c_str(), kRelResNBins, 0.5, kRelResNBins + 0.5);
    for (int k = 0; k < kRelResNBins; k++) {
        if (acc.count[k] <= 0) continue;
        h->SetBinContent(k + 1, acc.RMS(k));
        h->SetBinError(k + 1, acc.RMSError(k));
    }
    return h;
}
