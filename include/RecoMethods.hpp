#ifndef RECO_METHODS_HPP
#define RECO_METHODS_HPP

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "Math/Vector4D.h"
#include <vector>
#include <string>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <cctype>
#include <cmath>

using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;

// Structure to hold histograms for a reconstruction method
struct MethodHistograms {
    std::string name;
    TH1D* h_t_reco;
    TH2D* h_t_corr;
    TH1D* h_t_res;
    TH1D* h_theta;
    TH1D* h_xL;
    TH2D* h_xL_corr;
    TH1D* h_xL_res;
    TH1D* h_xpom;
    TH2D* h_xpom_corr;
    TH1D* h_xpom_res;
    TH1D* h_MX2;
    TH2D* h_MX2_corr;
    std::vector<double> t_truth_vec;
    std::vector<double> t_reco_vec;
    
    MethodHistograms(const std::string& method_name, const std::vector<Double_t>& t_bins);
    TGraph* MakeCorrelationGraph();
    void FillCorrelation(double t_truth, double t_reco);
    void Write();
};

// Structure to hold beam information
struct BeamInfo {
    P3MVector e_beam;
    P3MVector p_beam;
    double fMass_electron = 0.000511;
    double fMass_proton   = 0.938272;
};

// Calculate |t| from momentum transfer: BABE method
inline Double_t CalcT(const P3MVector& p_initial, const P3MVector& p_final) {
    return (p_final - p_initial).M2();
}

// Calculate |t| for eX method: t = (q - PX)^2
inline Double_t CalcT_eX(const P3MVector& q_gamma, const P3MVector& X_system) {
    return (q_gamma - X_system).M2();
}

// Calculate x_L = Pz_proton / E_proton
inline double CalcXL(const P3MVector& proton){
    return proton.Pz() / proton.E();
}

inline double CalcMX2(const P3MVector& X_system){
    return X_system.M2();
}

inline double CalcMX2_LPS(double xL, double xBj, double W){
    return (1.0 - xL*(1.0 + xBj)) * W * W;
}

////////////////////////////////////////////////////////////////////////
// Find and return scattered electron using EICRecon's identification
struct ScatteredElectronInfo {
    P3MVector p4;
    bool found;
    unsigned int index;
};

ScatteredElectronInfo GetScatteredElectron(
    TTreeReaderArray<int>& electron_index_array,
    TTreeReaderArray<float>& rpf_px, TTreeReaderArray<float>& rpf_py,
    TTreeReaderArray<float>& rpf_pz, double electron_mass
);

inline bool ParseBeamEnergiesFromFilename(const std::string& filePath,
                                          double& eBeamGeV,
                                          double& pBeamGeV) {
    auto isDelim = [](char c) {
        return c == '_' || c == '-' || c == '.';
    };

    const size_t slashPos = filePath.find_last_of("/\\");
    const std::string fileName = (slashPos == std::string::npos) ? filePath : filePath.substr(slashPos + 1);
    if (fileName.empty()) return false;

    for (size_t start = 0; start < fileName.size(); ++start) {
        const unsigned char c0 = static_cast<unsigned char>(fileName[start]);
        if (!std::isdigit(c0)) continue;

        size_t mid = start;
        bool dotSeen = false;
        while (mid < fileName.size()) {
            const char c = fileName[mid];
            if (std::isdigit(static_cast<unsigned char>(c))) {
                ++mid;
                continue;
            }
            if (c == '.' && !dotSeen) {
                dotSeen = true;
                ++mid;
                continue;
            }
            break;
        }
        if (mid <= start || mid >= fileName.size()) continue;
        if (fileName[mid] != 'x' && fileName[mid] != 'X') continue;

        const size_t rightStart = mid + 1;
        if (rightStart >= fileName.size()) continue;
        if (!std::isdigit(static_cast<unsigned char>(fileName[rightStart]))) continue;

        size_t end = rightStart;
        dotSeen = false;
        while (end < fileName.size()) {
            const char c = fileName[end];
            if (std::isdigit(static_cast<unsigned char>(c))) {
                ++end;
                continue;
            }
            if (c == '.' && !dotSeen) {
                dotSeen = true;
                ++end;
                continue;
            }
            break;
        }
        if (end <= rightStart) continue;

        if (start > 0 && !isDelim(fileName[start - 1])) continue;
        if (end < fileName.size() && !isDelim(fileName[end])) continue;

        try {
            const std::string eStr = fileName.substr(start, mid - start);
            const std::string pStr = fileName.substr(rightStart, end - rightStart);
            const double eCandidate = std::stod(eStr);
            const double pCandidate = std::stod(pStr);
            if (eCandidate > 0.0 && pCandidate > 0.0) {
                eBeamGeV = eCandidate;
                pBeamGeV = pCandidate;
                return true;
            }
        } catch (...) {
            // Continue scanning other potential tokens.
        }
    }
    return false;
}

inline bool InferBeamEnergiesFromFileList(const std::vector<std::string>& files,
                                          double& eBeamGeV,
                                          double& pBeamGeV,
                                          std::string* firstMatchedFile = nullptr,
                                          std::string* mismatchFile = nullptr) {
    bool found = false;
    double eRef = 0.0;
    double pRef = 0.0;

    for (const auto& file : files) {
        double eCand = 0.0;
        double pCand = 0.0;
        if (!ParseBeamEnergiesFromFilename(file, eCand, pCand)) continue;

        if (!found) {
            eRef = eCand;
            pRef = pCand;
            found = true;
            if (firstMatchedFile) *firstMatchedFile = file;
            continue;
        }

        const bool sameE = std::fabs(eCand - eRef) < 1e-9;
        const bool sameP = std::fabs(pCand - pRef) < 1e-9;
        if (!sameE || !sameP) {
            if (mismatchFile) *mismatchFile = file;
            return false;
        }
    }

    if (!found) return false;
    eBeamGeV = eRef;
    pBeamGeV = pRef;
    return true;
}

#endif // RECO_METHODS_HPP
