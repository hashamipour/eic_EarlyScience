// Inclusive DIS Skimmer: Q2/x/y/W2 and scattered-electron observables
// g++ DDIS_Skim_Final.cpp -o DDIS_Skim_Final $(root-config --cflags --glibs)
// ./DDIS_Skim_Final filelist.txt

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <TProfile2D.h>
#include <TParameter.h>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <random>
#include <cctype>
#include <sstream>
#include <limits>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include "Utility.hpp"
#include "RecoMethods.hpp"

// These are the crucial headers for the ROOT::Math objects
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"

using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

struct BinDef {
    int bin_id = -1;
    double Q2_min = -1.0;
    double Q2_max = -1.0;
    double beta_min = -1.0;
    double beta_max = -1.0;
    double xpom_min = -1.0;
    double xpom_max = -1.0;
};

static bool StartsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static std::string TrimWS(const std::string& s) {
    const char* ws = " \t\r\n";
    const size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

static std::vector<double> ParseInlineList(const std::string& line) {
    std::vector<double> values;
    const auto lbr = line.find('[');
    const auto rbr = line.find(']');
    if (lbr == std::string::npos || rbr == std::string::npos || rbr <= lbr) return values;
    std::string body = line.substr(lbr + 1, rbr - lbr - 1);
    for (char& c : body) {
        if (c == ',') c = ' ';
    }
    std::istringstream iss(body);
    double v = 0.0;
    while (iss >> v) {
        values.push_back(v);
    }
    return values;
}

static std::vector<double> ReadInlineListFromYAML(const std::string& path, const std::string& key) {
    std::ifstream in(path);
    if (!in.is_open()) return {};
    std::string line;
    const std::string prefix = key + ":";
    while (std::getline(in, line)) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = TrimWS(line);
        if (line.empty()) continue;
        if (StartsWith(line, prefix)) {
            return ParseInlineList(line);
        }
    }
    return {};
}

static std::vector<BinDef> ReadBinsFromYAML(const std::string& path) {
    std::vector<BinDef> bins;
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "ERROR: cannot open YAML file " << path << std::endl;
        return bins;
    }
    bool in_bins = false;
    BinDef cur;
    bool have_cur = false;
    std::string line;
    while (std::getline(in, line)) {
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = TrimWS(line);
        if (line.empty()) continue;
        if (StartsWith(line, "3D_bins:")) {
            in_bins = true;
            continue;
        }
        if (!in_bins) continue;
        if (StartsWith(line, "t_bins:")) continue;
        if (StartsWith(line, "W2_bins:")) continue;
        if (StartsWith(line, "-")) {
            if (have_cur) bins.push_back(cur);
            cur = BinDef();
            have_cur = true;
            const size_t pos = line.find("bin_id");
            if (pos != std::string::npos) {
                const size_t colon = line.find(':', pos);
                if (colon != std::string::npos) {
                    cur.bin_id = std::stoi(TrimWS(line.substr(colon + 1)));
                }
            }
            continue;
        }
        const size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        const std::string key = TrimWS(line.substr(0, colon));
        const std::string val_str = TrimWS(line.substr(colon + 1));
        if (val_str.empty()) continue;
        const double val = std::stod(val_str);
        if (key == "bin_id") cur.bin_id = static_cast<int>(val);
        else if (key == "Q2_min") cur.Q2_min = val;
        else if (key == "Q2_max") cur.Q2_max = val;
        else if (key == "beta_min") cur.beta_min = val;
        else if (key == "beta_max") cur.beta_max = val;
        else if (key == "x_min" || key == "x_pom_min" || key == "xpom_min") cur.xpom_min = val;
        else if (key == "x_max" || key == "x_pom_max" || key == "xpom_max") cur.xpom_max = val;
    }
    if (have_cur) bins.push_back(cur);
    return bins;
}

static std::vector<double> UniqueSorted(std::vector<double> vals) {
    std::sort(vals.begin(), vals.end());
    const double eps = 1e-12;
    vals.erase(std::unique(vals.begin(), vals.end(), [eps](double a, double b) {
        return std::fabs(a - b) < eps;
    }), vals.end());
    return vals;
}

struct CutflowConfig {
    bool apply_dis_cuts = true;// all cuts
    bool apply_proton_tag_cuts = true;
    bool use_em_for_dis_selection = false; // y-method switch can override later

    double q2_min = 1.0;
    double q2_max = std::numeric_limits<double>::infinity();
    double y_min = 0.01;
    double y_max = 0.95;

    bool apply_epz_cut = true;
    bool apply_epz_cut_to_truth = false;
    double epz_min = 16.0;// this is for 18x275 GeV2 , for 10x100 GeV2 it could be 16-24 GeV
    double epz_max = 24.0;

    double theta_rp_max_mrad = 5.0;
    double theta_b0_min_mrad = 5.5;
    double theta_b0_max_mrad = 20.0;
};

static bool PassDISKinematicCuts(const CutflowConfig& cfg,
                                 const bool validQ2,
                                 const bool validY,
                                 const double q2,
                                 const double y,
                                 const double epz,
                                 const bool isTruth) {
    if (!validQ2 || !validY) return false;
    if (!cfg.apply_dis_cuts) return true;

    if (!(q2 > cfg.q2_min && q2 < cfg.q2_max)) return false;
    if (!(y > cfg.y_min && y < cfg.y_max)) return false;

    const bool applyEPz = cfg.apply_epz_cut && (!isTruth || cfg.apply_epz_cut_to_truth);
    if (applyEPz) {
        if (!(std::isfinite(epz) && epz > cfg.epz_min && epz < cfg.epz_max)) return false;
    }

    return true;
}

static bool PassRPAcceptance(const CutflowConfig& cfg, const double theta_mrad) {
    if (!std::isfinite(theta_mrad)) return false;
    if (!cfg.apply_proton_tag_cuts) return true;
    return theta_mrad <= cfg.theta_rp_max_mrad;
}

static bool PassB0Acceptance(const CutflowConfig& cfg, const double theta_mrad) {
    if (!std::isfinite(theta_mrad)) return false;
    if (!cfg.apply_proton_tag_cuts) return true;
    return theta_mrad > cfg.theta_b0_min_mrad && theta_mrad < cfg.theta_b0_max_mrad;
}

static void CollectEdges(const std::vector<BinDef>& bins,
                         std::vector<double>& q2_edges,
                         std::vector<double>& beta_edges,
                         std::vector<double>& xpom_edges) {
    for (const auto& b : bins) {
        if (b.Q2_min > 0) q2_edges.push_back(b.Q2_min);
        if (b.Q2_max > 0) q2_edges.push_back(b.Q2_max);
        if (b.beta_min >= 0) beta_edges.push_back(b.beta_min);
        if (b.beta_max >= 0) beta_edges.push_back(b.beta_max);
        if (b.xpom_min > 0) xpom_edges.push_back(b.xpom_min);
        if (b.xpom_max > 0) xpom_edges.push_back(b.xpom_max);
    }
    q2_edges = UniqueSorted(q2_edges);
    beta_edges = UniqueSorted(beta_edges);
    xpom_edges = UniqueSorted(xpom_edges);
}

// Smooth DA->EM blend for W^2 reconstruction:
// alpha=1 at low W^2 (DA-dominated), alpha->0 at high W^2 (EM-dominated).
static constexpr double kW2BestW0 = 320.0;
static constexpr double kW2BestGamma = 4.0;

static double ComputeW2BestAlpha(const double w2_proxy) {
    if (!std::isfinite(w2_proxy) || w2_proxy <= 0.0) return 0.5;
    const double ratio = w2_proxy / kW2BestW0;
    const double powTerm = std::pow(ratio, kW2BestGamma);
    const double alpha = 1.0 / (1.0 + powTerm);
    return std::clamp(alpha, 0.0, 1.0);
}

static bool ComputeW2Best(const double w2_em,
                          const double w2_da,
                          double& w2_best,
                          double& alpha,
                          double& w2_proxy) {
    const bool validEM = std::isfinite(w2_em) && w2_em > 0.0;
    const bool validDA = std::isfinite(w2_da) && w2_da > 0.0;

    if (!validEM && !validDA) return false;

    if (validEM && validDA) {
        w2_proxy = std::sqrt(w2_em * w2_da);
        alpha = ComputeW2BestAlpha(w2_proxy);
        w2_best = alpha * w2_da + (1.0 - alpha) * w2_em;
        return std::isfinite(w2_best) && w2_best > 0.0;
    }

    if (validDA) {
        w2_proxy = w2_da;
        alpha = 1.0;
        w2_best = w2_da;
        return true;
    }

    w2_proxy = w2_em;
    alpha = 0.0;
    w2_best = w2_em;
    return true;
}

//==============================================================
// Fixed 3D binning for relative resolution vs global bin index
//==============================================================
static std::vector<double> kRelResQ2Bins;
static std::vector<double> kRelResBetaBins;
static std::vector<double> kRelResXpomBins;
static int kRelResNQ2 = 0;
static int kRelResNBeta = 0;
static int kRelResNXpom = 0;
static int kRelResNBins = 0;

static int FindBinIndex(const std::vector<double>& edges, double value) {
    if(edges.size() < 2) return -1;
    for(size_t i = 0; i + 1 < edges.size(); i++) {
        const double lo = edges[i];
        const double hi = edges[i + 1];
        const bool last = (i + 1 == edges.size() - 1);
        if(value >= lo && (value < hi || (last && value <= hi))) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

static int GetRelResGlobalBin(double q2, double xpom, double beta) {
    const int iQ2 = FindBinIndex(kRelResQ2Bins, q2);
    const int iXpom = FindBinIndex(kRelResXpomBins, xpom);
    const int iBeta = FindBinIndex(kRelResBetaBins, beta);
    if(iQ2 < 0 || iXpom < 0 || iBeta < 0) return -1;
    return iBeta + iXpom * kRelResNBeta + iQ2 * (kRelResNBeta * kRelResNXpom);
}

struct ResAccum {
    std::vector<double> sum;
    std::vector<double> sumsq;
    std::vector<int> count;

    explicit ResAccum(int nBins = 0)
        : sum(nBins, 0.0), sumsq(nBins, 0.0), count(nBins, 0) {}

    void Fill(int k, double value) {
        if(k < 0 || k >= static_cast<int>(sum.size())) return;
        sum[k] += value;
        sumsq[k] += value * value;
        count[k] += 1;
    }

    double RMS(int k) const {
        if(k < 0 || k >= static_cast<int>(sum.size()) || count[k] == 0) return 0.0;
        const double mean = sum[k] / static_cast<double>(count[k]);
        const double var = sumsq[k] / static_cast<double>(count[k]) - mean * mean;
        return (var > 0.0) ? std::sqrt(var) : 0.0;
    }

    double RMSError(int k) const {
        if(k < 0 || k >= static_cast<int>(sum.size()) || count[k] < 2) return 0.0;
        const double rms = RMS(k);
        return rms / std::sqrt(2.0 * (count[k] - 1.0));
    }
};

static TH1D* BuildRelResVsKHist(const std::string& name,
                                const std::string& title,
                                const ResAccum& acc) {
    TH1D* h = new TH1D(name.c_str(), title.c_str(), kRelResNBins, 0.5, kRelResNBins + 0.5);
    for(int k = 0; k < kRelResNBins; k++) {
        if(acc.count[k] <= 0) continue;
        h->SetBinContent(k + 1, acc.RMS(k));
        h->SetBinError(k + 1, acc.RMSError(k));
    }
    return h;
}

static bool WriteBinsTSV(const std::string& path, const std::vector<BinDef>& bins) {
    std::ofstream out(path);
    if (!out.is_open()) return false;
    out << "# bin_id\tQ2_min\tQ2_max\tbeta_min\tbeta_max\tx_pom_min\tx_pom_max\n";
    out << std::setprecision(8);
    for (const auto& b : bins) {
        out << b.bin_id << "\t"
            << b.Q2_min << "\t" << b.Q2_max << "\t"
            << b.beta_min << "\t" << b.beta_max << "\t"
            << b.xpom_min << "\t" << b.xpom_max << "\n";
    }
    return true;
}

static bool EdgesMatchWithinTolerance(double a, double b) {
    const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
    return std::abs(a - b) <= 1e-12 * scale;
}

static int FindLowerEdgeIndex(const std::vector<double>& edges, double value) {
    if (edges.size() < 2) return -1;
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        if (EdgesMatchWithinTolerance(edges[i], value)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

static int GetGlobalBinFromBinDef(const BinDef& b) {
    const int iQ2 = FindLowerEdgeIndex(kRelResQ2Bins, b.Q2_min);
    const int iXpom = FindLowerEdgeIndex(kRelResXpomBins, b.xpom_min);
    const int iBeta = FindLowerEdgeIndex(kRelResBetaBins, b.beta_min);
    if (iQ2 < 0 || iXpom < 0 || iBeta < 0) return -1;

    if (iQ2 + 1 >= static_cast<int>(kRelResQ2Bins.size()) ||
        iXpom + 1 >= static_cast<int>(kRelResXpomBins.size()) ||
        iBeta + 1 >= static_cast<int>(kRelResBetaBins.size())) {
        return -1;
    }

    if (!EdgesMatchWithinTolerance(kRelResQ2Bins[iQ2 + 1], b.Q2_max) ||
        !EdgesMatchWithinTolerance(kRelResXpomBins[iXpom + 1], b.xpom_max) ||
        !EdgesMatchWithinTolerance(kRelResBetaBins[iBeta + 1], b.beta_max)) {
        return -1;
    }

    return iBeta + iXpom * kRelResNBeta + iQ2 * (kRelResNBeta * kRelResNXpom);
}

static bool WriteBinOccupancyYAML(const std::string& path,
                                  const std::vector<BinDef>& bins,
                                  const std::vector<int>& occupancy) {
    std::ofstream out(path);
    if (!out.is_open()) return false;

    long long totalAssigned = 0;
    for (const int c : occupancy) totalAssigned += c;

    out << "metadata:\n";
    out << "  description: \"Per-bin occupancy from truth global bin assignment (Q2_truth, x_pom_truth, beta_truth).\"\n";
    out << "  n_bins: " << bins.size() << "\n";
    out << "  total_assigned_events: " << totalAssigned << "\n";
    out << "  k_index_convention: \"one_based\"\n";
    out << "bins:\n";
    out << std::setprecision(8);

    for (const auto& b : bins) {
        const int kZeroBased = GetGlobalBinFromBinDef(b);
        const int occ = (kZeroBased >= 0 && kZeroBased < static_cast<int>(occupancy.size()))
                            ? occupancy[kZeroBased]
                            : 0;
        out << "  - bin_id: " << b.bin_id << "\n";
        out << "    k_index: " << (kZeroBased >= 0 ? kZeroBased + 1 : -1) << "\n";
        out << "    occupancy: " << occ << "\n";
        out << "    Q2_min: " << b.Q2_min << "\n";
        out << "    Q2_max: " << b.Q2_max << "\n";
        out << "    x_pom_min: " << b.xpom_min << "\n";
        out << "    x_pom_max: " << b.xpom_max << "\n";
        out << "    beta_min: " << b.beta_min << "\n";
        out << "    beta_max: " << b.beta_max << "\n";
    }
    return true;
}

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};
double MASS_PROTON   = fMass_proton;
double MASS_ELECTRON = fMass_electron;

// Global afterburner correction parameters
Float_t fXAngle{-0.025};
RotationX rotAboutX;
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

void undoAfterburnAndCalc(P3MVector& p, P3MVector& k){
    P3MVector p_beam(fXAngle*p.E(), 0., p.E(), p.M());
    P3MVector e_beam(0., 0., -k.E(), k.M());
    
    P3MVector CoM_boost = p_beam + e_beam;
    vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
    
    p_beam = boost(p_beam, vBoostToCoM);
    e_beam = boost(e_beam, vBoostToCoM);
    
    Float_t fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
    Float_t fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());
    
    rotAboutY = RotationY(fRotY);
    rotAboutX = RotationX(fRotX);
    
    p_beam = rotAboutY(p_beam);
    p_beam = rotAboutX(p_beam);
    e_beam = rotAboutY(e_beam);
    e_beam = rotAboutX(e_beam);
    
    P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
    vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());
    
    p_beam = boost(p_beam, vBoostToHoF);
    e_beam = boost(e_beam, vBoostToHoF);
    
    p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
    k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}

void undoAfterburn(P3MVector& a){
    a = boost(a, vBoostToCoM);
    a = rotAboutY(a);
    a = rotAboutX(a);
    a = boost(a, vBoostToHoF);
}

// E - pz calculation for matched reco/truth particles
static void CalculateSumEPz_Matched(
    TTreeReaderArray<float>& re_px,
    TTreeReaderArray<float>& re_py,
    TTreeReaderArray<float>& re_pz,
    TTreeReaderArray<float>& re_energy,
    TTreeReaderArray<double>& mc_px,
    TTreeReaderArray<double>& mc_py,
    TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass,
    TTreeReaderArray<unsigned int>& assoc_rec_id,
    TTreeReaderArray<unsigned int>& assoc_sim_id,
    double& sumEPz_truth,
    double& sumEPz_reco
) {
    sumEPz_truth = 0.0;
    sumEPz_reco = 0.0;

    for (unsigned int i = 0; i < re_energy.GetSize(); i++) {
        int mc_idx = -1;
        for (unsigned int j = 0; j < assoc_rec_id.GetSize(); j++) {
            if (assoc_rec_id[j] == i) {
                mc_idx = static_cast<int>(assoc_sim_id[j]);
                break;
            }
        }
        if (mc_idx < 0 || mc_idx >= static_cast<int>(mc_px.GetSize())) continue;

        const double E_reco = re_energy[i];
        const double pz_reco = re_pz[i];
        sumEPz_reco += (E_reco - pz_reco);

        const double px_mc = mc_px[mc_idx];
        const double py_mc = mc_py[mc_idx];
        const double pz_mc = mc_pz[mc_idx];
        const double m_mc = mc_mass[mc_idx];
        const double E_mc = TMath::Sqrt(px_mc * px_mc + py_mc * py_mc + pz_mc * pz_mc + m_mc * m_mc);
        sumEPz_truth += (E_mc - pz_mc);
    }
}

static float GetArrayValue(const TTreeReaderArray<float>& arr, float fallback = -999.0f) {
    return (arr.GetSize() > 0) ? arr[0] : fallback;
}

static float GetArrayValue(const TTreeReaderArray<float>* arr, float fallback = -999.0f) {
    if (!arr || arr->GetSize() == 0) return fallback;
    return (*arr)[0];
}

static bool ComputeTruthKinematics(const BeamInfo& beams,
                                   const TTreeReaderArray<double>& mc_px,
                                   const TTreeReaderArray<double>& mc_py,
                                   const TTreeReaderArray<double>& mc_pz,
                                   const TTreeReaderArray<double>& mc_mass,
                                   int scat_idx,
                                   float& Q2,
                                   float& x,
                                   float& y,
                                   float& W) {
    if (scat_idx < 0 || scat_idx >= static_cast<int>(mc_px.GetSize())) return false;

    const double kx = beams.e_beam.Px();
    const double ky = beams.e_beam.Py();
    const double kz = beams.e_beam.Pz();
    const double kE = beams.e_beam.E();

    const double px = beams.p_beam.Px();
    const double py = beams.p_beam.Py();
    const double pz = beams.p_beam.Pz();
    const double pE = beams.p_beam.E();

    const double kpx = mc_px[scat_idx];
    const double kpy = mc_py[scat_idx];
    const double kpz = mc_pz[scat_idx];
    const double km  = mc_mass[scat_idx];
    const double kEprime = std::sqrt(kpx*kpx + kpy*kpy + kpz*kpz + km*km);

    const double qx = kx - kpx;
    const double qy = ky - kpy;
    const double qz = kz - kpz;
    const double qE = kE - kEprime;

    const double q2 = qE*qE - (qx*qx + qy*qy + qz*qz);
    const double Q2calc = -q2;
    if (!std::isfinite(Q2calc)) return false;

    const double pDotq = pE*qE - (px*qx + py*qy + pz*qz);
    const double pDotk = pE*kE - (px*kx + py*ky + pz*kz);
    if (pDotq <= 0.0 || pDotk <= 0.0) return false;

    const double xcalc = Q2calc / (2.0 * pDotq);
    const double ycalc = pDotq / pDotk;

    const double W2calc = (pE + qE)*(pE + qE)
                        - ((px + qx)*(px + qx) + (py + qy)*(py + qy) + (pz + qz)*(pz + qz));
    const double Wcalc = (W2calc > 0.0) ? std::sqrt(W2calc) : -1.0;

    Q2 = static_cast<float>(Q2calc);
    x  = static_cast<float>(xcalc);
    y  = static_cast<float>(ycalc);
    W  = static_cast<float>(Wcalc);
    return std::isfinite(x) && std::isfinite(y);
}

static std::string TrimLine(const std::string& line) {
    const auto start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    const auto end = line.find_last_not_of(" \t\r\n");
    return line.substr(start, end - start + 1);
}

static std::vector<std::string> ReadFileList(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::string> files;
    if (!in.is_open()) {
        std::cerr << "Error: Could not open file list " << path << std::endl;
        return files;
    }
    std::string line;
    while (std::getline(in, line)) {
        line = TrimLine(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        files.push_back(line);
    }
    return files;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fileList.txt> [N] <bins.yaml>" << std::endl;
        return 1;
    }
    TString fileList = argv[1];
    int sampleN = -1;
    std::string yamlPath;
    if (argc == 3) {
        yamlPath = argv[2];
    } else {
        char* end = nullptr;
        long parsed = std::strtol(argv[2], &end, 10);
        if (!(end && *end == '\0' && parsed > 0)) {
            std::cerr << "ERROR: Expected integer N as second argument when providing 3 args." << std::endl;
            return 1;
        }
        sampleN = static_cast<int>(parsed);
        yamlPath = argv[3];
    }
    if (yamlPath.empty()) {
        std::cerr << "ERROR: bins.yaml is mandatory." << std::endl;
        return 1;
    }

    std::cout<< " __ __ __ __ __ __ __ __ __ __" <<std::endl;
    std::cout<< "|                             |"<<std::endl;
    std::cout<< "|   Inclusive DIS Skim        |"<<std::endl;
    std::cout<< "|   Q2/x/y/W2 + e- vars       |"<<std::endl;
    std::cout<< "|__ __ __ __ __ __ __ __ __ __|"<<std::endl;
    std::cout<< "\nInput filelist: " << fileList <<std::endl;
    if (!yamlPath.empty()) {
        std::cout << "Using YAML bins: " << yamlPath << std::endl;
    }

    std::vector<BinDef> yaml_bins = ReadBinsFromYAML(yamlPath);
    if (yaml_bins.empty()) {
        std::cerr << "ERROR: No bins found in YAML: " << yamlPath << std::endl;
        return 1;
    }
    std::vector<double> yaml_q2_edges;
    std::vector<double> yaml_beta_edges;
    std::vector<double> yaml_xpom_edges;
    CollectEdges(yaml_bins, yaml_q2_edges, yaml_beta_edges, yaml_xpom_edges);
    if (yaml_q2_edges.size() < 2 || yaml_beta_edges.size() < 2 || yaml_xpom_edges.size() < 2) {
        std::cerr << "ERROR: Invalid YAML bin edges in " << yamlPath << std::endl;
        return 1;
    }
    kRelResQ2Bins = yaml_q2_edges;
    kRelResBetaBins = yaml_beta_edges;
    kRelResXpomBins = yaml_xpom_edges;
    kRelResNQ2 = static_cast<int>(kRelResQ2Bins.size()) - 1;
    kRelResNBeta = static_cast<int>(kRelResBetaBins.size()) - 1;
    kRelResNXpom = static_cast<int>(kRelResXpomBins.size()) - 1;
    kRelResNBins = kRelResNQ2 * kRelResNXpom * kRelResNBeta;
    std::vector<double> yaml_t_edges = ReadInlineListFromYAML(yamlPath, "t_bins");
    std::vector<double> yaml_w2_edges = ReadInlineListFromYAML(yamlPath, "W2_bins");

    std::vector<std::string> allFiles = ReadFileList(fileList.Data());
    if (allFiles.empty()) {
        std::cerr << "Error: file list is empty after filtering." << std::endl;
        return 1;
    }

    std::vector<std::string> selectedFiles = allFiles;
    if (sampleN > 0 && sampleN < static_cast<int>(allFiles.size())) {
        std::mt19937 rng(std::random_device{}());
        std::shuffle(selectedFiles.begin(), selectedFiles.end(), rng);
        selectedFiles.resize(sampleN);
        std::cout << "Sampling " << sampleN << " of " << allFiles.size() << " files." << std::endl;
    } else if (sampleN > 0) {
        std::cout << "Requested N=" << sampleN << " >= available files (" << allFiles.size()
                  << "); using full list." << std::endl;
    }

    //---------------------------------------------------------
    // CREATE TCHAIN AND OUTPUT ROOT FILE
    //---------------------------------------------------------
    TChain* events = new TChain("events");
    Int_t nFiles{0};
    std::vector<std::string> addedFiles;
    for (const auto& fileName : selectedFiles) {
        TString tmp = fileName;
        // Skip existence check for XRootD - TFile::Open handles remote files
        if(tmp.BeginsWith("root://")) {
            events->Add(tmp);
            nFiles++;
            addedFiles.push_back(fileName);
            continue;
        }
        // Local file check
        if (!std::filesystem::exists(tmp.Data())) {
            std::cerr << "\nWarning: File does not exist: " << fileName << std::endl;
            continue;
        }
        events->Add(tmp);
        nFiles++;
        addedFiles.push_back(fileName);
    }
    std::cout<<"\nNo. of files: "<<nFiles<<"; no. of events: "<<events->GetEntries()<<std::endl;
    
    // Create output file
    TFile* outputFile = new TFile("DDIS_Combined_output.root", "RECREATE");

    // Create TTree for event-level data
    TTree* tree = new TTree("Q2_tree", "Q2 Kinematics Data");
    float out_Q2_truth, out_Q2_EM, out_Q2_DA, out_Q2_Sigma;
    float out_x_truth, out_x_EM, out_x_DA, out_x_Sigma;
    float out_y_truth, out_y_EM, out_y_DA, out_y_Sigma;
    tree->Branch("Q2_truth",   &out_Q2_truth,  "Q2_truth/F");
    tree->Branch("Q2_EM",      &out_Q2_EM,     "Q2_EM/F");
    tree->Branch("Q2_DA",      &out_Q2_DA,     "Q2_DA/F");
    tree->Branch("Q2_Sigma",   &out_Q2_Sigma,  "Q2_Sigma/F");
    tree->Branch("x_truth",    &out_x_truth,   "x_truth/F");
    tree->Branch("x_EM",       &out_x_EM,      "x_EM/F");
    tree->Branch("x_DA",       &out_x_DA,      "x_DA/F");
    tree->Branch("x_Sigma",    &out_x_Sigma,   "x_Sigma/F");
    tree->Branch("y_truth",    &out_y_truth,   "y_truth/F");
    tree->Branch("y_EM",       &out_y_EM,      "y_EM/F");
    tree->Branch("y_DA",       &out_y_DA,      "y_DA/F");
    tree->Branch("y_Sigma",    &out_y_Sigma,   "y_Sigma/F");

    //---------------------------------------------------------
    // DECLARE OUTPUT HISTOGRAMS (Q2/xy part)
    //---------------------------------------------------------

    int n_bins = 10;
    std::vector<Double_t> bin_edges_Q2;
    if (!yaml_q2_edges.empty()) {
        bin_edges_Q2.assign(yaml_q2_edges.begin(), yaml_q2_edges.end());
    } else {
        bin_edges_Q2 = GetRoundedLogBins(3.4, 150.0, n_bins);
    }
    n_bins = bin_edges_Q2.size()-1;
    
    // x histograms (0..1)
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, 25);
    TH1D* h_x_EM     = new TH1D("x_EM",     "electron method;x_{Bj}", x_bins.size()-1, x_bins.data());
    TH1D* h_x_DA     = new TH1D("x_DA",     "DA method;x_{Bj}",       x_bins.size()-1, x_bins.data());
    TH1D* h_x_Sigma  = new TH1D("x_Sigma",  "Sigma method;x_{Bj}",    x_bins.size()-1, x_bins.data());
    TH1D* h_x_truth  = new TH1D("x_truth",  "truth;x_{Bj}",           x_bins.size()-1, x_bins.data());

    // y histograms (0..1)
    int n_y_bins = 10;
    TH1D* h_y_EM     = new TH1D("y_EM",     "electron method;y", n_y_bins, 0.0, 1.0);
    TH1D* h_y_DA     = new TH1D("y_DA",     "DA method;y",       n_y_bins, 0.0, 1.0);
    TH1D* h_y_Sigma  = new TH1D("y_Sigma",  "Sigma method;y",    n_y_bins, 0.0, 1.0);
    TH1D* h_y_truth  = new TH1D("y_truth",  "truth;y",           n_y_bins, 0.0, 1.0);

    // Event density in (x_Bj, Q2) with finer binning for phase space coverage
    const int n_x_bins_density = 120;
    const int n_q2_bins_density = 120;
    const double q2_min_density = 1.0;
    const double q2_max_density = 1000.0;
    std::vector<Double_t> x_bins_density = GetLogBins(1.0e-4, 1.0, n_x_bins_density);
    std::vector<Double_t> q2_bins_density = GetLogBins(q2_min_density, q2_max_density, n_q2_bins_density);
    TH2D* h_xQ2_truth = new TH2D("xQ2_truth",
                                "Event Density;x_{Bj};Q^{2} [GeV^{2}]",
                                x_bins_density.size()-1, x_bins_density.data(),
                                q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_xQ2_reco = new TH2D("xQ2_reco",
                               "Event Density (Reco);x_{Bj};Q^{2} [GeV^{2}]",
                               x_bins_density.size()-1, x_bins_density.data(),
                               q2_bins_density.size()-1, q2_bins_density.data());

    // Relative resolution histograms
    TH1D* h_RelRes_x_EM     = new TH1D("x_RelRes_EM",     "electron method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_x_DA     = new TH1D("x_RelRes_DA",     "DA method;#frac{x(Reco)-x(MC)}{x(MC)}"      , 41, -0.5, 0.5);
    TH1D* h_RelRes_x_Sigma  = new TH1D("x_RelRes_Sigma",  "Sigma method;#frac{x(Reco)-x(MC)}{x(MC)}"   , 41, -0.5, 0.5);

    TH1D* h_RelRes_y_EM     = new TH1D("y_RelRes_EM",     "electron method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_y_DA     = new TH1D("y_RelRes_DA",     "DA method;#frac{y(Reco)-y(MC)}{y(MC)}"      , 51, -0.5, 0.5);
    TH1D* h_RelRes_y_Sigma  = new TH1D("y_RelRes_Sigma", "Sigma method;#frac{y(Reco)-y(MC)}{y(MC)}"    , 51, -0.5, 0.5);

    // 2D binned relres vs truth
    int n_binned = 51;
    TH2D* h_RelRes_x_binned_EM     = new TH2D("x_RelRes_binned_EM",     "x: truth vs rel. res (EM);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_DA     = new TH2D("x_RelRes_binned_DA",     "x: truth vs rel. res (DA);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_Sigma = new TH2D("x_RelRes_binned_Sigma", "x: truth vs rel. res (Sigma);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",   
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);

    TH2D* h_RelRes_y_binned_EM     = new TH2D("y_RelRes_binned_EM",     "y: truth vs rel. res (EM);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_DA     = new TH2D("y_RelRes_binned_DA",     "y: truth vs rel. res (DA);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.3, 0.3);
    TH2D* h_RelRes_y_binned_Sigma = new TH2D("y_RelRes_binned_Sigma", "y: truth vs rel. res (Sigma);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",   n_y_bins, 0.0, 1.0, n_binned-20, -0.6, 0.3);

    // Legacy 2D correlation matrices (truth vs reco) for EM/DA/Sigma
    TH2D* h_Corr_x_EM = new TH2D("x_Corr_EM", "x correlation (EM);x_{truth};x_{EM}",
                                 x_bins.size() - 1, x_bins.data(), x_bins.size() - 1, x_bins.data());
    TH2D* h_Corr_x_DA = new TH2D("x_Corr_DA", "x correlation (DA);x_{truth};x_{DA}",
                                 x_bins.size() - 1, x_bins.data(), x_bins.size() - 1, x_bins.data());
    TH2D* h_Corr_x_Sigma = new TH2D("x_Corr_Sigma", "x correlation (Sigma);x_{truth};x_{Sigma}",
                                    x_bins.size() - 1, x_bins.data(), x_bins.size() - 1, x_bins.data());
    TH2D* h_Corr_y_EM = new TH2D("y_Corr_EM", "y correlation (EM);y_{truth};y_{EM}",
                                 n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_DA = new TH2D("y_Corr_DA", "y correlation (DA);y_{truth};y_{DA}",
                                 n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_Sigma = new TH2D("y_Corr_Sigma", "y correlation (Sigma);y_{truth};y_{Sigma}",
                                    n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);

    // TProfile2D for Q2 resolution - BINNED in (x, Q2) and DISPLAYED in (x, Q2)
    TProfile2D* h_Q2_RelRes_vs_xy_EM = new TProfile2D("Q2_RelRes_vs_xy_EM",
        "Q^{2} Rel. Res. binned in (x, Q^{2}) (EM);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_Q2_RelRes_vs_xy_DA = new TProfile2D("Q2_RelRes_vs_xy_DA",
        "Q^{2} Rel. Res. binned in (x, Q^{2}) (DA);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_Q2_RelRes_vs_xy_Sigma = new TProfile2D("Q2_RelRes_vs_xy_Sigma",
        "Q^{2} Rel. Res. binned in (x, Q^{2}) (Sigma);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");

    // TProfile2D for x resolution - BINNED in (x, Q2) and DISPLAYED in (x, Q2)
    TProfile2D* h_x_RelRes_vs_xQ2_EM = new TProfile2D("x_RelRes_vs_xQ2_EM",
        "x Rel. Res. binned in (x, Q^{2}) (EM);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_x_RelRes_vs_xQ2_DA = new TProfile2D("x_RelRes_vs_xQ2_DA",
        "x Rel. Res. binned in (x, Q^{2}) (DA);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_x_RelRes_vs_xQ2_Sigma = new TProfile2D("x_RelRes_vs_xQ2_Sigma",
        "x Rel. Res. binned in (x, Q^{2}) (Sigma);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");

    // TProfile2D for y resolution - BINNED in (x, Q2) and DISPLAYED in (x, Q2)
    TProfile2D* h_y_RelRes_vs_xQ2_EM = new TProfile2D("y_RelRes_vs_xQ2_EM",
        "y Rel. Res. binned in (x, Q^{2}) (EM);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_y_RelRes_vs_xQ2_DA = new TProfile2D("y_RelRes_vs_xQ2_DA",
        "y Rel. Res. binned in (x, Q^{2}) (DA);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_y_RelRes_vs_xQ2_Sigma = new TProfile2D("y_RelRes_vs_xQ2_Sigma",
        "y Rel. Res. binned in (x, Q^{2}) (Sigma);x;Q^{2} [GeV^{2}]",
        x_bins.size()-1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");

    TH1D* h_RelRes_Q2_EM = new TH1D("Q2_RelRes_EM","electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}",71,-0.15,0.15);
    TH1D* h_RelRes_Q2_DA = new TH1D("Q2_RelRes_DA","DA method;#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",71,-0.15,0.15);
    TH1D* h_RelRes_Q2_Sigma = new TH1D("Q2_RelRes_Sigma","Sigma method;#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",43,-0.4,0.15);
    
    TH2D* h_RelRes_Q2_binned_EM = new TH2D("Q2_RelRes_binned_EM",";Q^{2} [GeV^{2}];#frac{Q^{2}(EM)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_DA = new TH2D("Q2_RelRes_binned_DA",";Q^{2} [GeV^{2}];#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_Sigma = new TH2D("Q2_RelRes_binned_Sigma",";Q^{2} [GeV^{2}];#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_Corr_Q2_EM = new TH2D("Corr_Q2_EM", ";Q^{2}_{MC};Q^{2}_{EM}",
                                  n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_DA = new TH2D("Corr_Q2_DA", ";Q^{2}_{MC};Q^{2}_{DA}",
                                  n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_Sigma = new TH2D("Corr_Q2_Sigma", ";Q^{2}_{MC};Q^{2}_{Sigma}",
                                     n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());

    // Relative resolution vs global 3D bin index (Q2, x_pom, beta)
    ResAccum res_Q2_EM_k(kRelResNBins);
    ResAccum res_Q2_DA_k(kRelResNBins);
    ResAccum res_Q2_Sigma_k(kRelResNBins);
    ResAccum res_xpom_EM_B0_k(kRelResNBins);
    ResAccum res_xpom_DA_B0_k(kRelResNBins);
    ResAccum res_xpom_Sigma_B0_k(kRelResNBins);
    ResAccum res_xpom_EM_RP_k(kRelResNBins);
    ResAccum res_xpom_DA_RP_k(kRelResNBins);
    ResAccum res_xpom_Sigma_RP_k(kRelResNBins);
    ResAccum res_beta_EM_B0_k(kRelResNBins);
    ResAccum res_beta_DA_B0_k(kRelResNBins);
    ResAccum res_beta_Sigma_B0_k(kRelResNBins);
    ResAccum res_beta_EM_RP_k(kRelResNBins);
    ResAccum res_beta_DA_RP_k(kRelResNBins);
    ResAccum res_beta_Sigma_RP_k(kRelResNBins);
    std::vector<int> occupancy_truth_k(kRelResNBins, 0);
    TH1D* h_phase_bin_gen = new TH1D("phase_bin_gen",
                                     "Generated counts per global bin;Global bin k;Counts",
                                     kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_meas = new TH1D("phase_bin_meas",
                                      "Measured counts per global bin (after selections);Global bin k;Counts",
                                      kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_gen_meas_same = new TH1D("phase_bin_gen_meas_same",
                                               "Generated and measured in same global bin;Global bin k;Counts",
                                               kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_gen_setA = new TH1D("phase_bin_gen_setA",
                                          "Generated counts per global bin (Set A);Global bin k;Counts",
                                          kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_meas_setA = new TH1D("phase_bin_meas_setA",
                                           "Measured counts per global bin (Set A);Global bin k;Counts",
                                           kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_gen_meas_same_setA = new TH1D("phase_bin_gen_meas_same_setA",
                                                    "Generated and measured in same global bin (Set A);Global bin k;Counts",
                                                    kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_gen_setB = new TH1D("phase_bin_gen_setB",
                                          "Generated counts per global bin (Set B);Global bin k;Counts",
                                          kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_meas_setB = new TH1D("phase_bin_meas_setB",
                                           "Measured counts per global bin (Set B);Global bin k;Counts",
                                           kRelResNBins, 0.5, kRelResNBins + 0.5);
    TH1D* h_phase_bin_gen_meas_same_setB = new TH1D("phase_bin_gen_meas_same_setB",
                                                    "Generated and measured in same global bin (Set B);Global bin k;Counts",
                                                    kRelResNBins, 0.5, kRelResNBins + 0.5);

    
    TH1D* h_Q2_truth    = new TH1D("h_Q2_truth","Q^2;# of events",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_EM       = new TH1D("h_Q2_EM",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_DA       = new TH1D("h_Q2_DA",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_Sigma   = new TH1D("h_Q2_Sigma",";Q^{2}",n_bins, bin_edges_Q2.data());

    //---------------------------------------------------------
    // NEW INCLUSIVE DIS HISTOGRAMS
    //---------------------------------------------------------

    // W² histograms
    std::vector<double> w2_bins_yaml = yaml_w2_edges;
    std::vector<Double_t> w2_bins;
    if (w2_bins_yaml.size() >= 2) {
        w2_bins.assign(w2_bins_yaml.begin(), w2_bins_yaml.end());
    } else {
        const int n_w2_bins = 60;
        const double w2_min = 10.0;
        const double w2_max = 1.0e4;
        w2_bins = GetLogBins(w2_min, w2_max, n_w2_bins);
    }
    TH1D* h_W2_EM    = new TH1D("W2_EM",    "Electron method;W^{2} [GeV^{2}]", w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_DA    = new TH1D("W2_DA",    "DA method;W^{2} [GeV^{2}]",       w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_Best  = new TH1D("W2_Best",  "Best blend (DA#rightarrowEM);W^{2} [GeV^{2}]", w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_Sigma = new TH1D("W2_Sigma", "Sigma method;W^{2} [GeV^{2}]",    w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_truth = new TH1D("W2_truth", "Truth;W^{2} [GeV^{2}]",           w2_bins.size()-1, w2_bins.data());
    std::vector<double> w2_relres_em_edges;
    {
        // Piecewise test binning:
        // coarse [-0.2, -0.02], tight [-0.02, 0.02], coarse [0.02, 0.2].
        const int n_coarse_left = 18;   // width 0.01
        const int n_tight_mid = 20;     // 
        const int n_coarse_right = 18;  // width 0.01

        auto append_linear_segment = [&w2_relres_em_edges](double x0, double x1, int nBins, bool includeStart) {
            for (int i = includeStart ? 0 : 1; i <= nBins; ++i) {
                w2_relres_em_edges.push_back(x0 + (x1 - x0) * static_cast<double>(i) / nBins);
            }
        };

        append_linear_segment(-0.1, -0.005, n_coarse_left, true);
        append_linear_segment(-0.005, 0.005, n_tight_mid, false);
        append_linear_segment(0.005, 0.1, n_coarse_right, false);
    }
    TH1D* h_RelRes_W2_EM = new TH1D(
        "W2_RelRes_EM",
        "electron method;#frac{W^{2}_{EM}-W^{2}_{MC}}{W^{2}_{MC}}",
        static_cast<int>(w2_relres_em_edges.size()) - 1,
        w2_relres_em_edges.data()
    );
    TH1D* h_RelRes_W2_DA = new TH1D("W2_RelRes_DA", "DA method;#frac{W^{2}_{DA}-W^{2}_{MC}}{W^{2}_{MC}}", 99, -1.0, 1.0);
    TH1D* h_RelRes_W2_Best = new TH1D("W2_RelRes_Best", "Best blend;#frac{W^{2}_{best}-W^{2}_{MC}}{W^{2}_{MC}}", 101, -1.0, 1.0);
    TH1D* h_RelRes_W2_Sigma = new TH1D("W2_RelRes_Sigma", "Sigma method;#frac{W^{2}_{#Sigma}-W^{2}_{MC}}{W^{2}_{MC}}", 101, -1.0, 1.0);

    TH2D* h_RelRes_W2_binned_EM = new TH2D("W2_RelRes_binned_EM",
        ";W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{EM}-W^{2}_{MC}}{W^{2}_{MC}}",
        w2_bins.size()-1, w2_bins.data(), n_binned, -1.0, 1.0);
    TH2D* h_RelRes_W2_binned_DA = new TH2D("W2_RelRes_binned_DA",
        ";W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{DA}-W^{2}_{MC}}{W^{2}_{MC}}",
        w2_bins.size()-1, w2_bins.data(), n_binned, -1.0, 1.0);
    TH2D* h_RelRes_W2_binned_Best = new TH2D("W2_RelRes_binned_Best",
        ";W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{best}-W^{2}_{MC}}{W^{2}_{MC}}",
        w2_bins.size()-1, w2_bins.data(), n_binned, -1.0, 1.0);
    TH2D* h_RelRes_W2_binned_Sigma = new TH2D("W2_RelRes_binned_Sigma",
        ";W^{2}_{MC} [GeV^{2}];#frac{W^{2}_{#Sigma}-W^{2}_{MC}}{W^{2}_{MC}}",
        w2_bins.size()-1, w2_bins.data(), n_binned, -1.0, 1.0);

    // M_X^2 histograms (hadronic invariant mass squared, excluding scattered e- and leading proton)
    const int n_mx2_bins = 120;
    const double mx2_min = 1.0e-3;
    const double mx2_max = 1000.0;
    std::vector<Double_t> mx2_bins = GetLogBins(mx2_min, mx2_max, n_mx2_bins);
    TH1D* h_MX2_truth = new TH1D("MX2_truth", "Truth M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts",
                                mx2_bins.size() - 1, mx2_bins.data());
    TH1D* h_MX2_reco  = new TH1D("MX2_reco",  "Reco M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts",
                                mx2_bins.size() - 1, mx2_bins.data());
    TH2D* h_MX2_corr  = new TH2D("MX2_corr",
                                "M_{X}^{2} Truth vs Reco;M_{X,truth}^{2} [GeV^{2}];M_{X,reco}^{2} [GeV^{2}]",
                                mx2_bins.size() - 1, mx2_bins.data(),
                                mx2_bins.size() - 1, mx2_bins.data());
    std::vector<double> mx2_relres_edges;
    {
        // Piecewise binning with dense bins around zero:
        // coarse [-1.0, -0.02], tight [-0.02, 0.02], coarse [0.02, 1.0].
        const int n_coarse_left = 18;
        const int n_tight_mid = 20;
        const int n_coarse_right = 18;

        auto append_linear_segment = [&mx2_relres_edges](double x0, double x1, int nBins, bool includeStart) {
            for (int i = includeStart ? 0 : 1; i <= nBins; ++i) {
                mx2_relres_edges.push_back(x0 + (x1 - x0) * static_cast<double>(i) / nBins);
            }
        };

        append_linear_segment(-1.0, -0.05, n_coarse_left, true);
        append_linear_segment(-0.05, 0.05, n_tight_mid, false);
        append_linear_segment(0.05, 1.0, n_coarse_right, false);
    }
    TH1D* h_MX2_RelRes = new TH1D("MX2_RelRes",
                                  "M_{X}^{2} relative resolution;#frac{M_{X,reco}^{2}-M_{X,truth}^{2}}{M_{X,truth}^{2}};Counts",
                                  static_cast<int>(mx2_relres_edges.size()) - 1, mx2_relres_edges.data());
    TH2D* h_MX2_RelRes_binned = new TH2D("MX2_RelRes_binned",
                                         "M_{X}^{2} relative resolution;M_{X,truth}^{2} [GeV^{2}];#frac{M_{X,reco}^{2}-M_{X,truth}^{2}}{M_{X,truth}^{2}}",
                                         mx2_bins.size() - 1, mx2_bins.data(), n_binned, -1.0, 1.0);
    TProfile2D* h_MX2_RelRes_vs_MX2Q2 = new TProfile2D("MX2_RelRes_vs_MX2Q2",
                                                       "M_{X}^{2} Rel. Res. binned in (M_{X}^{2}, Q^{2});M_{X}^{2} [GeV^{2}];Q^{2} [GeV^{2}]",
                                                       mx2_bins.size() - 1, mx2_bins.data(),
                                                       n_bins, bin_edges_Q2.data(),
                                                       "s");

    // E-pz and eta_{max} distributions (legacy plots)
    TH1D* h_EPz_truth = new TH1D("h_EPz_truth", "MC Truth Sum(E-p_{z}) - Matched Particles;#Sigma(E-p_{z}) [GeV];Counts", 80, 0.0, 40.0);
    TH1D* h_EPz = new TH1D("h_EPz", "Reco Sum(E-p_{z}) - Matched Particles;#Sigma(E-p_{z}) [GeV];Counts", 80, 0.0, 40.0);
    TH2D* h_EPz_2D = new TH2D("h_EPz_2D",
                              "E-p_{z} Truth vs Reco;#Sigma(E-p_{z})_{truth} [GeV];#Sigma(E-p_{z})_{reco} [GeV]",
                              120, 0.0, 40.0, 120, 0.0, 40.0);
    TH1D* h_eta_max = new TH1D("h_eta_max", "Maximum Pseudorapidity per Event (Reco);#eta_{max};Counts", 60, -5.0, 7.0);
    TH1D* h_eta_max_truth = new TH1D("h_eta_max_truth", "Maximum Pseudorapidity per Event (Truth);#eta_{max};Counts", 60, -5.0, 7.0);

    // Scattered electron quantities (reco)
    TH1D* h_Ep_e     = new TH1D("Ep_e",     "Scattered electron energy;E'_{e} [GeV]",  50, 0, 20);
    TH1D* h_phi_e    = new TH1D("phi_e",    "Scattered electron #phi;#phi_{e} [rad]",  50, -TMath::Pi(), TMath::Pi());
    TH1D* h_pT_e     = new TH1D("pT_e",     "Scattered electron p_{T};p_{T}^{e} [GeV]", 50, 0, 10);

    // Truth electron quantities
    TH1D* h_Ep_e_truth    = new TH1D("Ep_e_truth",    "Truth E'_{e};E'_{e} [GeV]",       50, 0, 20);
    TH1D* h_phi_e_truth   = new TH1D("phi_e_truth",   "Truth #phi_{e};#phi_{e} [rad]",   50, -TMath::Pi(), TMath::Pi());
    TH1D* h_pT_e_truth    = new TH1D("pT_e_truth",    "Truth p_{T}^{e};p_{T}^{e} [GeV]", 50, 0, 10);

    // Correlation plots for electron quantities are stored as TGraphs (unbinned)

    // Response matrices (for unfolding) - using Electron method
    TH2D* h_Response_Q2 = new TH2D("Response_Q2", "Q^{2} Response Matrix;Q^{2}_{truth} [GeV^{2}];Q^{2}_{reco} [GeV^{2}]",
                                    n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());
    TH2D* h_Response_x  = new TH2D("Response_x",  "x Response Matrix;x_{truth};x_{reco}",
                                    x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Response_y  = new TH2D("Response_y",  "y Response Matrix;y_{truth};y_{reco}",
                                    n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Response_W2_EM = new TH2D("Response_W2_EM",
                                      "W^{2} Response Matrix (EM);W^{2}_{truth} [GeV^{2}];W^{2}_{reco} [GeV^{2}]",
                                      w2_bins.size()-1, w2_bins.data(), w2_bins.size()-1, w2_bins.data());
    TH2D* h_Response_W2_DA = new TH2D("Response_W2_DA",
                                      "W^{2} Response Matrix (DA);W^{2}_{truth} [GeV^{2}];W^{2}_{reco} [GeV^{2}]",
                                      w2_bins.size()-1, w2_bins.data(), w2_bins.size()-1, w2_bins.data());
    TH2D* h_Response_W2_Best = new TH2D("Response_W2_Best",
                                        "W^{2} Response Matrix (best);W^{2}_{truth} [GeV^{2}];W^{2}_{reco} [GeV^{2}]",
                                        w2_bins.size()-1, w2_bins.data(), w2_bins.size()-1, w2_bins.data());
    TH2D* h_Response_W2_Sigma = new TH2D("Response_W2_Sigma",
                                         "W^{2} Response Matrix (#Sigma);W^{2}_{truth} [GeV^{2}];W^{2}_{reco} [GeV^{2}]",
                                         w2_bins.size()-1, w2_bins.data(), w2_bins.size()-1, w2_bins.data());

    // Unbinned correlation graphs (reco vs truth)
    TGraph* g_Q2_EM = new TGraph(); g_Q2_EM->SetName("g_Q2_EM");
    TGraph* g_Q2_DA = new TGraph(); g_Q2_DA->SetName("g_Q2_DA");
    TGraph* g_Q2_Sigma = new TGraph(); g_Q2_Sigma->SetName("g_Q2_Sigma");

    TGraph* g_x_EM = new TGraph(); g_x_EM->SetName("g_x_EM");
    TGraph* g_x_DA = new TGraph(); g_x_DA->SetName("g_x_DA");
    TGraph* g_x_Sigma = new TGraph(); g_x_Sigma->SetName("g_x_Sigma");

    TGraph* g_y_EM = new TGraph(); g_y_EM->SetName("g_y_EM");
    TGraph* g_y_DA = new TGraph(); g_y_DA->SetName("g_y_DA");
    TGraph* g_y_Sigma = new TGraph(); g_y_Sigma->SetName("g_y_Sigma");

    TGraph* g_W2_EM = new TGraph(); g_W2_EM->SetName("g_W2_EM");
    TGraph* g_W2_DA = new TGraph(); g_W2_DA->SetName("g_W2_DA");
    TGraph* g_W2_Best = new TGraph(); g_W2_Best->SetName("g_W2_Best");
    TGraph* g_W2_Sigma = new TGraph(); g_W2_Sigma->SetName("g_W2_Sigma");

    TGraph* g_Ep_e = new TGraph(); g_Ep_e->SetName("g_Ep_e");
    TGraph* g_phi_e = new TGraph(); g_phi_e->SetName("g_phi_e");
    TGraph* g_pT_e = new TGraph(); g_pT_e->SetName("g_pT_e");

    int n_g_Q2_EM = 0, n_g_Q2_DA = 0, n_g_Q2_Sigma = 0;
    int n_g_x_EM = 0, n_g_x_DA = 0, n_g_x_Sigma = 0;
    int n_g_y_EM = 0, n_g_y_DA = 0, n_g_y_Sigma = 0;
    int n_g_W2_EM = 0, n_g_W2_DA = 0, n_g_W2_Best = 0, n_g_W2_Sigma = 0;
    int n_g_Ep_e = 0, n_g_phi_e = 0, n_g_pT_e = 0;

    //---------------------------------------------------------
    // DIFFRACTIVE: Mandelstam t
    //---------------------------------------------------------
    std::vector<Double_t> t_bins;
    if (yaml_t_edges.size() >= 2) {
        t_bins.assign(yaml_t_edges.begin(), yaml_t_edges.end());
    } else {
        std::vector<Double_t> t_bins_low = GetLogBins(1e-3, 0.5, 20);
        std::vector<Double_t> t_bins_high = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.6, 1.8, 2.0};
        t_bins = t_bins_low;
        for(size_t i = 1; i < t_bins_high.size(); i++){
            t_bins.push_back(t_bins_high[i]);
        }
    }

    TH1D* h_t_MC = new TH1D("t_MC", "Truth Mandelstam t;|t| [GeV^{2}];Counts",
                            t_bins.size()-1, t_bins.data());
    TH1D* h_t_B0 = new TH1D("t_B0", "B0 Reco Mandelstam t;|t| [GeV^{2}];Counts",
                            t_bins.size()-1, t_bins.data());
    TH1D* h_t_RP_histo = new TH1D("t_RP_histo", "RP Reco Mandelstam t;|t| [GeV^{2}];Counts",
                                 t_bins.size()-1, t_bins.data());

    TH1D* h_dsigma_dt_MC = new TH1D("dsigma_dt_MC", "Truth d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_dsigma_dt_B0 = new TH1D("dsigma_dt_B0", "B0 Reco d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_dsigma_dt_RP = new TH1D("dsigma_dt_RP", "RP Reco d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_dsigma_dt_Sum = nullptr;

    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_all_TS = new TH1D("theta_all_TS", "All Truth-Seeded Proton Angles;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_RP = new TH1D("theta_RP", "RP Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);

    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_B0 = new TH1D("xL_B0", "B0 Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_RP = new TH1D("xL_RP", "RP Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_res_B0 = new TH1D("xL_res_B0", "B0 x_{L} Resolution;(x_{L,reco}-x_{L,truth})/x_{L,truth};Counts", 31, -0.5, 0.5);
    TH1D* h_xL_res_RP = new TH1D("xL_res_RP", "RP x_{L} Resolution;(x_{L,reco}-x_{L,truth})/x_{L,truth};Counts", 31, -0.1, 0.1);
    TH2D* h_xL_corr_B0 = new TH2D("xL_corr_B0", "B0 x_{L} Correlation;Truth x_{L};Reco x_{L}", 100, 0.7, 1.1, 100, 0.7, 1.1);
    TH2D* h_xL_corr_RP = new TH2D("xL_corr_RP", "RP x_{L} Correlation;Truth x_{L};Reco x_{L}", 100, 0.7, 1.1, 100, 0.7, 1.1);
    TH2D* h_xL_RelRes_binned_B0 = new TH2D("xL_RelRes_binned_B0",
                                           "x_{L} relative resolution (B0);x_{L,truth};#frac{x_{L,reco}-x_{L,truth}}{x_{L,truth}}",
                                           30, 0.75, 1.05, 100, -2.0, 2.0);
    TH2D* h_xL_RelRes_binned_RP = new TH2D("xL_RelRes_binned_RP",
                                           "x_{L} relative resolution (RP);x_{L,truth};#frac{x_{L,reco}-x_{L,truth}}{x_{L,truth}}",
                                           30, 0.75, 1.05, 100, -2.0, 2.0);

    // x_pom histograms from definition: x_pom = (Q^2 + M_X^2 - t)/(Q^2 + W^2 - m_p^2)
    std::vector<Double_t> xpom_bins;
    xpom_bins.assign(yaml_xpom_edges.begin(), yaml_xpom_edges.end());
    std::vector<Double_t> xpom_centers;
    xpom_centers.reserve(xpom_bins.size() - 1);
    for (size_t i = 0; i + 1 < xpom_bins.size(); ++i) {
        xpom_centers.push_back(std::sqrt(xpom_bins[i] * xpom_bins[i + 1]));
    }

    std::cout << "\n[x_pom] bin edges:";
    for (size_t i = 0; i < xpom_bins.size(); ++i) {
        std::cout << " " << std::scientific << xpom_bins[i];
    }
    std::cout << "\n[x_pom] bin centers:";
    for (size_t i = 0; i < xpom_centers.size(); ++i) {
        std::cout << " " << std::scientific << xpom_centers[i];
    }
    std::cout << "\n" << std::defaultfloat;

    std::vector<Double_t> beta_bins;
    beta_bins.assign(yaml_beta_edges.begin(), yaml_beta_edges.end());
    TH1D* h_xpom_truth_all = new TH1D("xpom_truth_all", "Truth x_{pom} (All);x_{pom};Counts",
                                     xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_truth_B0 = new TH1D("xpom_truth_B0", "Truth x_{pom} (B0);x_{pom};Counts",
                                    xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_truth_RP = new TH1D("xpom_truth_RP", "Truth x_{pom} (RP);x_{pom};Counts",
                                    xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_EM_all = new TH1D("xpom_reco_EM_all", "Reco x_{pom} EM (All);x_{pom};Counts",
                                       xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_EM_B0 = new TH1D("xpom_reco_EM_B0", "Reco x_{pom} EM (B0);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_EM_RP = new TH1D("xpom_reco_EM_RP", "Reco x_{pom} EM (RP);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_DA_all = new TH1D("xpom_reco_DA_all", "Reco x_{pom} DA (All);x_{pom};Counts",
                                       xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_DA_B0 = new TH1D("xpom_reco_DA_B0", "Reco x_{pom} DA (B0);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_DA_RP = new TH1D("xpom_reco_DA_RP", "Reco x_{pom} DA (RP);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_Sigma_all = new TH1D("xpom_reco_Sigma_all", "Reco x_{pom} #Sigma (All);x_{pom};Counts",
                                          xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_Sigma_B0 = new TH1D("xpom_reco_Sigma_B0", "Reco x_{pom} #Sigma (B0);x_{pom};Counts",
                                         xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_Sigma_RP = new TH1D("xpom_reco_Sigma_RP", "Reco x_{pom} #Sigma (RP);x_{pom};Counts",
                                         xpom_bins.size()-1, xpom_bins.data());
    // Legacy x_{pom} definitions: from x_{L} and from definition, plus correlations
    TH1D* h_xpom_MC = new TH1D("xpom_MC", "Truth x_{pom} (from x_{L});x_{pom};Counts", xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_B0 = new TH1D("xpom_B0", "B0 Reco x_{pom} (from x_{L});x_{pom};Counts", xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_RP = new TH1D("xpom_RP", "RP Reco x_{pom} (from x_{L});x_{pom};Counts", xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_def_MC = new TH1D("xpom_def_MC", "Truth x_{pom} (from definition);x_{pom};Counts", xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_def_B0 = new TH1D("xpom_def_B0", "B0 Reco x_{pom} (from definition);x_{pom};Counts", xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_def_RP = new TH1D("xpom_def_RP", "RP Reco x_{pom} (from definition);x_{pom};Counts", xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_xpom_comp_MC = new TH2D("xpom_comp_MC", "Truth: x_{pom}(1-x_{L}) vs x_{pom}(def);x_{pom} (1-x_{L});x_{pom} (definition)",
                                    xpom_bins.size()-1, xpom_bins.data(), xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_xpom_comp_B0 = new TH2D("xpom_comp_B0", "B0: x_{pom}(1-x_{L}) vs x_{pom}(def);x_{pom} (1-x_{L});x_{pom} (definition)",
                                    xpom_bins.size()-1, xpom_bins.data(), xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_xpom_comp_RP = new TH2D("xpom_comp_RP", "RP: x_{pom}(1-x_{L}) vs x_{pom}(def);x_{pom} (1-x_{L});x_{pom} (definition)",
                                    xpom_bins.size()-1, xpom_bins.data(), xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_xpom_corr_B0 = new TH2D("xpom_corr_B0", "B0 x_{pom} Correlation (from x_{L});Truth x_{pom};Reco x_{pom}",
                                    xpom_bins.size()-1, xpom_bins.data(), xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_xpom_corr_RP = new TH2D("xpom_corr_RP", "RP x_{pom} Correlation (from x_{L});Truth x_{pom};Reco x_{pom}",
                                    xpom_bins.size()-1, xpom_bins.data(), xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_res_B0 = new TH1D("xpom_res_B0", "B0 x_{pom} Resolution (from x_{L});(x_{pom,reco}-x_{pom,truth})/x_{pom,truth};Counts", 30, -1.5, 1.5);
    TH1D* h_xpom_res_RP = new TH1D("xpom_res_RP", "RP x_{pom} Resolution (from x_{L});(x_{pom,reco}-x_{pom,truth})/x_{pom,truth};Counts", 60, -1.5, 1.0);
    TH2D* h_xpom_RelRes_binned_B0 = new TH2D("xpom_RelRes_binned_B0",
                                             "x_{pom} relative resolution (B0);x_{pom,truth};#frac{x_{pom,reco}-x_{pom,truth}}{x_{pom,truth}}",
                                             xpom_bins.size()-1, xpom_bins.data(), 21, -2.0, 2.0);
    TH2D* h_xpom_RelRes_binned_RP = new TH2D("xpom_RelRes_binned_RP",
                                             "x_{pom} relative resolution (RP);x_{pom,truth};#frac{x_{pom,reco}-x_{pom,truth}}{x_{pom,truth}}",
                                             xpom_bins.size()-1, xpom_bins.data(), 21, -2.0, 2.0);
    TH2D* h_xpom_RelRes_binned_W2Best_B0 = new TH2D("xpom_RelRes_binned_W2Best_B0",
                                                    "x_{pom} relative resolution (W^{2}_{best}, B0);x_{pom,truth};#frac{x_{pom,reco}-x_{pom,truth}}{x_{pom,truth}}",
                                                    xpom_bins.size()-1, xpom_bins.data(), 21, -2.0, 2.0);
    TH2D* h_xpom_RelRes_binned_W2Best_RP = new TH2D("xpom_RelRes_binned_W2Best_RP",
                                                    "x_{pom} relative resolution (W^{2}_{best}, RP);x_{pom,truth};#frac{x_{pom,reco}-x_{pom,truth}}{x_{pom,truth}}",
                                                    xpom_bins.size()-1, xpom_bins.data(), 21, -2.0, 2.0);

    TH2D* h_Response_xpom_EM_B0 = new TH2D("Response_xpom_EM_B0", "x_{pom} Response EM (B0);Truth x_{pom};Reco x_{pom}",
                                          xpom_bins.size()-1, xpom_bins.data(),
                                          xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_EM_RP = new TH2D("Response_xpom_EM_RP", "x_{pom} Response EM (RP);Truth x_{pom};Reco x_{pom}",
                                          xpom_bins.size()-1, xpom_bins.data(),
                                          xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_DA_B0 = new TH2D("Response_xpom_DA_B0", "x_{pom} Response DA (B0);Truth x_{pom};Reco x_{pom}",
                                          xpom_bins.size()-1, xpom_bins.data(),
                                          xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_DA_RP = new TH2D("Response_xpom_DA_RP", "x_{pom} Response DA (RP);Truth x_{pom};Reco x_{pom}",
                                          xpom_bins.size()-1, xpom_bins.data(),
                                          xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_Sigma_B0 = new TH2D("Response_xpom_Sigma_B0", "x_{pom} Response #Sigma (B0);Truth x_{pom};Reco x_{pom}",
                                             xpom_bins.size()-1, xpom_bins.data(),
                                             xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_Sigma_RP = new TH2D("Response_xpom_Sigma_RP", "x_{pom} Response #Sigma (RP);Truth x_{pom};Reco x_{pom}",
                                             xpom_bins.size()-1, xpom_bins.data(),
                                             xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_W2Best_B0 = new TH2D("Response_xpom_W2Best_B0", "x_{pom} Response W^{2}_{best} (B0);Truth x_{pom};Reco x_{pom}",
                                              xpom_bins.size()-1, xpom_bins.data(),
                                              xpom_bins.size()-1, xpom_bins.data());
    TH2D* h_Response_xpom_W2Best_RP = new TH2D("Response_xpom_W2Best_RP", "x_{pom} Response W^{2}_{best} (RP);Truth x_{pom};Reco x_{pom}",
                                              xpom_bins.size()-1, xpom_bins.data(),
                                              xpom_bins.size()-1, xpom_bins.data());

    const int n_beta_resp_bins = static_cast<int>(beta_bins.size()) - 1;
    TH2D* h_Response_beta_EM_B0 = new TH2D("Response_beta_EM_B0", "#beta Response EM (B0);Truth #beta;Reco #beta",
                                          n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_EM_RP = new TH2D("Response_beta_EM_RP", "#beta Response EM (RP);Truth #beta;Reco #beta",
                                          n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_DA_B0 = new TH2D("Response_beta_DA_B0", "#beta Response DA (B0);Truth #beta;Reco #beta",
                                          n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_DA_RP = new TH2D("Response_beta_DA_RP", "#beta Response DA (RP);Truth #beta;Reco #beta",
                                          n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_Sigma_B0 = new TH2D("Response_beta_Sigma_B0", "#beta Response #Sigma (B0);Truth #beta;Reco #beta",
                                             n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_Sigma_RP = new TH2D("Response_beta_Sigma_RP", "#beta Response #Sigma (RP);Truth #beta;Reco #beta",
                                             n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_W2Best_B0 = new TH2D("Response_beta_W2Best_B0", "#beta Response W^{2}_{best} (B0);Truth #beta;Reco #beta",
                                              n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());
    TH2D* h_Response_beta_W2Best_RP = new TH2D("Response_beta_W2Best_RP", "#beta Response W^{2}_{best} (RP);Truth #beta;Reco #beta",
                                              n_beta_resp_bins, beta_bins.data(), n_beta_resp_bins, beta_bins.data());

    TH1D* h_beta_truth_all = new TH1D("beta_truth_all", "Truth #beta (All);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_truth_B0 = new TH1D("beta_truth_B0", "Truth #beta (B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_truth_RP = new TH1D("beta_truth_RP", "Truth #beta (RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_EM_all = new TH1D("beta_reco_EM_all", "Reco #beta EM (All);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_EM_B0 = new TH1D("beta_reco_EM_B0", "Reco #beta EM (B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_EM_RP = new TH1D("beta_reco_EM_RP", "Reco #beta EM (RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_DA_all = new TH1D("beta_reco_DA_all", "Reco #beta DA (All);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_DA_B0 = new TH1D("beta_reco_DA_B0", "Reco #beta DA (B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_DA_RP = new TH1D("beta_reco_DA_RP", "Reco #beta DA (RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_Sigma_all = new TH1D("beta_reco_Sigma_all", "Reco #beta #Sigma (All);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_Sigma_B0 = new TH1D("beta_reco_Sigma_B0", "Reco #beta #Sigma (B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_Sigma_RP = new TH1D("beta_reco_Sigma_RP", "Reco #beta #Sigma (RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    // Legacy beta histograms/correlations
    TH1D* h_beta_MC = new TH1D("beta_MC", "Truth #beta (from x_{pom} def);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_B0 = new TH1D("beta_B0", "B0 Reco #beta (EM,x_{pom} def);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_RP = new TH1D("beta_RP", "RP Reco #beta (EM,x_{pom} def);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_res_B0 = new TH1D("beta_res_B0", "B0 #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts", 40, -0.4, 0.4);
    TH1D* h_beta_res_RP = new TH1D("beta_res_RP", "RP #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts", 40, -0.5, 0.5);
    TH2D* h_beta_corr_B0 = new TH2D("beta_corr_B0", "B0 #beta Correlation;Truth #beta;Reco #beta",
                                    beta_bins.size()-1, beta_bins.data(), beta_bins.size()-1, beta_bins.data());
    TH2D* h_beta_corr_RP = new TH2D("beta_corr_RP", "RP #beta Correlation;Truth #beta;Reco #beta",
                                    beta_bins.size()-1, beta_bins.data(), beta_bins.size()-1, beta_bins.data());
    TH2D* h_beta_RelRes_binned_B0 = new TH2D("beta_RelRes_binned_B0",
                                             "#beta relative resolution (B0);#beta_{truth};#frac{#beta_{reco}-#beta_{truth}}{#beta_{truth}}",
                                             beta_bins.size()-1, beta_bins.data(), 21, -1.0, 1.0);
    TH2D* h_beta_RelRes_binned_RP = new TH2D("beta_RelRes_binned_RP",
                                             "#beta relative resolution (RP);#beta_{truth};#frac{#beta_{reco}-#beta_{truth}}{#beta_{truth}}",
                                             beta_bins.size()-1, beta_bins.data(), 21, -1.0, 1.0);
    TH2D* h_beta_RelRes_binned_W2Best_B0 = new TH2D("beta_RelRes_binned_W2Best_B0",
                                                    "#beta relative resolution (W^{2}_{best}, B0);#beta_{truth};#frac{#beta_{reco}-#beta_{truth}}{#beta_{truth}}",
                                                    beta_bins.size()-1, beta_bins.data(), 21, -1.0, 1.0);
    TH2D* h_beta_RelRes_binned_W2Best_RP = new TH2D("beta_RelRes_binned_W2Best_RP",
                                                    "#beta relative resolution (W^{2}_{best}, RP);#beta_{truth};#frac{#beta_{reco}-#beta_{truth}}{#beta_{truth}}",
                                                    beta_bins.size()-1, beta_bins.data(), 21, -1.0, 1.0);
    TH2D* h_beta_vs_Q2 = new TH2D("beta_vs_Q2", "#beta vs Q^{2};Q^{2} [GeV^{2}];#beta",
                                  n_bins, bin_edges_Q2.data(), beta_bins.size()-1, beta_bins.data());
    TH2D* h_beta_vs_xpom = new TH2D("beta_vs_xpom", "#beta vs x_{pom};x_{pom};#beta",
                                    xpom_bins.size()-1, xpom_bins.data(), beta_bins.size()-1, beta_bins.data());
    TH2D* h_beta_vs_t = new TH2D("beta_vs_t", "#beta vs |t|;|t| [GeV^{2}];#beta",
                                 t_bins.size()-1, t_bins.data(), beta_bins.size()-1, beta_bins.data());

    // Diffractive 2D profile maps for circle plots
    std::vector<Double_t> xL_bins_profile(31);
    for (size_t i = 0; i < xL_bins_profile.size(); ++i) {
        xL_bins_profile[i] = 0.75 + (1.05 - 0.75) * static_cast<double>(i) / static_cast<double>(xL_bins_profile.size() - 1);
    }
    TProfile2D* h_t_RelRes_vs_xpomQ2_B0 = new TProfile2D("t_RelRes_vs_xpomQ2_B0",
                                                         "|t| rel. res. vs (x_{pom},Q^{2}) B0;x_{pom};Q^{2} [GeV^{2}]",
                                                         xpom_bins.size()-1, xpom_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_t_RelRes_vs_xpomQ2_RP = new TProfile2D("t_RelRes_vs_xpomQ2_RP",
                                                         "|t| rel. res. vs (x_{pom},Q^{2}) RP;x_{pom};Q^{2} [GeV^{2}]",
                                                         xpom_bins.size()-1, xpom_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_xpom_RelRes_vs_xpomQ2_B0 = new TProfile2D("xpom_RelRes_vs_xpomQ2_B0",
                                                            "x_{pom} rel. res. vs (x_{pom},Q^{2}) B0;x_{pom};Q^{2} [GeV^{2}]",
                                                            xpom_bins.size()-1, xpom_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_xpom_RelRes_vs_xpomQ2_RP = new TProfile2D("xpom_RelRes_vs_xpomQ2_RP",
                                                            "x_{pom} rel. res. vs (x_{pom},Q^{2}) RP;x_{pom};Q^{2} [GeV^{2}]",
                                                            xpom_bins.size()-1, xpom_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_beta_RelRes_vs_betaQ2_B0 = new TProfile2D("beta_RelRes_vs_betaQ2_B0",
                                                            "#beta rel. res. vs (#beta,Q^{2}) B0;#beta;Q^{2} [GeV^{2}]",
                                                            beta_bins.size()-1, beta_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_beta_RelRes_vs_betaQ2_RP = new TProfile2D("beta_RelRes_vs_betaQ2_RP",
                                                            "#beta rel. res. vs (#beta,Q^{2}) RP;#beta;Q^{2} [GeV^{2}]",
                                                            beta_bins.size()-1, beta_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_xL_RelRes_vs_xLQ2_B0 = new TProfile2D("xL_RelRes_vs_xLQ2_B0",
                                                        "x_{L} rel. res. vs (x_{L},Q^{2}) B0;x_{L};Q^{2} [GeV^{2}]",
                                                        xL_bins_profile.size()-1, xL_bins_profile.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_xL_RelRes_vs_xLQ2_RP = new TProfile2D("xL_RelRes_vs_xLQ2_RP",
                                                        "x_{L} rel. res. vs (x_{L},Q^{2}) RP;x_{L};Q^{2} [GeV^{2}]",
                                                        xL_bins_profile.size()-1, xL_bins_profile.data(), n_bins, bin_edges_Q2.data(), "s");

    // Triple differential cross sections d^3sigma / (dQ^2 d#beta dx_{pom})
    const int n_beta_3d_bins = 5;
    double beta_3d_bins[n_beta_3d_bins + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    const int n_Q2_3d_bins = 10;
    std::vector<Double_t> Q2_3d_bins = GetLogBins(1.0, 100.0, n_Q2_3d_bins);
    const int n_xpom_3d_bins = 8;
    std::vector<Double_t> xpom_3d_bins = GetLogBins(1.0e-3, 0.1, n_xpom_3d_bins);
    TH3D* h_d3sigma_MC = new TH3D("d3sigma_dQ2dbeta_dxpom_MC",
                                  "Truth d^{3}#sigma/(dQ^{2}d#betadx_{pom});Q^{2} [GeV^{2}];#beta;x_{pom}",
                                  n_Q2_3d_bins, Q2_3d_bins.data(),
                                  n_beta_3d_bins, beta_3d_bins,
                                  n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_B0 = new TH3D("d3sigma_dQ2dbeta_dxpom_B0",
                                  "B0 Reco d^{3}#sigma/(dQ^{2}d#betadx_{pom});Q^{2} [GeV^{2}];#beta;x_{pom}",
                                  n_Q2_3d_bins, Q2_3d_bins.data(),
                                  n_beta_3d_bins, beta_3d_bins,
                                  n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_RP = new TH3D("d3sigma_dQ2dbeta_dxpom_RP",
                                  "RP Reco d^{3}#sigma/(dQ^{2}d#betadx_{pom});Q^{2} [GeV^{2}];#beta;x_{pom}",
                                  n_Q2_3d_bins, Q2_3d_bins.data(),
                                  n_beta_3d_bins, beta_3d_bins,
                                  n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_truth_mc_B0 = new TH3D("d3sigma_dQ2dbeta_dxpom_truth_mc_B0",
                                           "Truth d^{3}#sigma inputs (MC subsample, B0);Q^{2} [GeV^{2}];#beta;x_{pom}",
                                           n_Q2_3d_bins, Q2_3d_bins.data(),
                                           n_beta_3d_bins, beta_3d_bins,
                                           n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_reco_mc_B0 = new TH3D("d3sigma_dQ2dbeta_dxpom_reco_mc_B0",
                                          "Reco d^{3}#sigma inputs (MC subsample, B0);Q^{2} [GeV^{2}];#beta;x_{pom}",
                                          n_Q2_3d_bins, Q2_3d_bins.data(),
                                          n_beta_3d_bins, beta_3d_bins,
                                          n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_same_mc_B0 = new TH3D("d3sigma_dQ2dbeta_dxpom_same_mc_B0",
                                          "Generated and measured in same d^{3}#sigma bin (MC subsample, B0);Q^{2} [GeV^{2}];#beta;x_{pom}",
                                          n_Q2_3d_bins, Q2_3d_bins.data(),
                                          n_beta_3d_bins, beta_3d_bins,
                                          n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_truth_mc_RP = new TH3D("d3sigma_dQ2dbeta_dxpom_truth_mc_RP",
                                           "Truth d^{3}#sigma inputs (MC subsample, RP);Q^{2} [GeV^{2}];#beta;x_{pom}",
                                           n_Q2_3d_bins, Q2_3d_bins.data(),
                                           n_beta_3d_bins, beta_3d_bins,
                                           n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_reco_mc_RP = new TH3D("d3sigma_dQ2dbeta_dxpom_reco_mc_RP",
                                          "Reco d^{3}#sigma inputs (MC subsample, RP);Q^{2} [GeV^{2}];#beta;x_{pom}",
                                          n_Q2_3d_bins, Q2_3d_bins.data(),
                                          n_beta_3d_bins, beta_3d_bins,
                                          n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_same_mc_RP = new TH3D("d3sigma_dQ2dbeta_dxpom_same_mc_RP",
                                          "Generated and measured in same d^{3}#sigma bin (MC subsample, RP);Q^{2} [GeV^{2}];#beta;x_{pom}",
                                          n_Q2_3d_bins, Q2_3d_bins.data(),
                                          n_beta_3d_bins, beta_3d_bins,
                                          n_xpom_3d_bins, xpom_3d_bins.data());
    TH3D* h_d3sigma_Sum = nullptr;

    // --------------------------------------------------------
    // Efficiency-correction inputs (MC subsample + pseudo-data)
    // --------------------------------------------------------
    TH1D* h_t_truth_mc = new TH1D("t_truth_mc", "Truth |t| (MC);|t| [GeV^{2}];Counts",
                                 t_bins.size()-1, t_bins.data());
    TH1D* h_t_reco_mc = new TH1D("t_reco_mc", "Reco |t| (MC);|t| [GeV^{2}];Counts",
                                t_bins.size()-1, t_bins.data());
    TH1D* h_t_reco_pdata = new TH1D("t_reco_pdata", "Reco |t| (Pseudo-data);|t| [GeV^{2}];Counts",
                                   t_bins.size()-1, t_bins.data());
    TH1D* h_t_truth_pdata_all = new TH1D("t_truth_pdata_all", "Truth |t| (Pseudo-data);|t| [GeV^{2}];Counts",
                                        t_bins.size()-1, t_bins.data());

    TH1D* h_beta_truth_mc = new TH1D("beta_truth_mc", "Truth #beta (MC);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_mc = new TH1D("beta_reco_mc", "Reco #beta (MC);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_pdata = new TH1D("beta_reco_pdata", "Reco #beta (Pseudo-data);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_truth_pdata_all = new TH1D("beta_truth_pdata_all", "Truth #beta (Pseudo-data);#beta;Counts", beta_bins.size()-1, beta_bins.data());

    TH1D* h_xpom_truth_mc = new TH1D("xpom_truth_mc", "Truth x_{pom} (MC);x_{pom};Counts",
                                    xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_mc = new TH1D("xpom_reco_mc", "Reco x_{pom} (MC);x_{pom};Counts",
                                   xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_pdata = new TH1D("xpom_reco_pdata", "Reco x_{pom} (Pseudo-data);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_truth_pdata_all = new TH1D("xpom_truth_pdata_all", "Truth x_{pom} (Pseudo-data);x_{pom};Counts",
                                           xpom_bins.size()-1, xpom_bins.data());

    TH1D* h_Q2_truth_mc = new TH1D("Q2_truth_mc", "Truth Q^{2} (MC);Q^{2} [GeV^{2}];Counts",
                                  n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_reco_mc = new TH1D("Q2_reco_mc", "Reco Q^{2} (MC);Q^{2} [GeV^{2}];Counts",
                                 n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_reco_pdata = new TH1D("Q2_reco_pdata", "Reco Q^{2} (Pseudo-data);Q^{2} [GeV^{2}];Counts",
                                    n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_truth_pdata = new TH1D("Q2_truth_pdata", "Truth Q^{2} (Pseudo-data);Q^{2} [GeV^{2}];Counts",
                                     n_bins, bin_edges_Q2.data());

    TH1D* h_MX2_truth_mc = new TH1D("MX2_truth_mc", "Truth M_{X}^{2} (MC);M_{X}^{2} [GeV^{2}];Counts",
                                   mx2_bins.size() - 1, mx2_bins.data());
    TH1D* h_MX2_reco_mc = new TH1D("MX2_reco_mc", "Reco M_{X}^{2} (MC);M_{X}^{2} [GeV^{2}];Counts",
                                  mx2_bins.size() - 1, mx2_bins.data());
    TH1D* h_MX2_reco_pdata = new TH1D("MX2_reco_pdata", "Reco M_{X}^{2} (Pseudo-data);M_{X}^{2} [GeV^{2}];Counts",
                                     mx2_bins.size() - 1, mx2_bins.data());
    TH1D* h_MX2_truth_pdata = new TH1D("MX2_truth_pdata", "Truth M_{X}^{2} (Pseudo-data);M_{X}^{2} [GeV^{2}];Counts",
                                      mx2_bins.size() - 1, mx2_bins.data());

    TH1D* h_x_truth_mc = new TH1D("x_truth_mc", "Truth x_{Bj} (MC);x_{Bj};Counts",
                                 x_bins.size()-1, x_bins.data());
    TH1D* h_x_reco_mc = new TH1D("x_reco_mc", "Reco x_{Bj} (MC);x_{Bj};Counts",
                                x_bins.size()-1, x_bins.data());
    TH1D* h_x_reco_pdata = new TH1D("x_reco_pdata", "Reco x_{Bj} (Pseudo-data);x_{Bj};Counts",
                                   x_bins.size()-1, x_bins.data());
    TH1D* h_x_truth_pdata = new TH1D("x_truth_pdata", "Truth x_{Bj} (Pseudo-data);x_{Bj};Counts",
                                    x_bins.size()-1, x_bins.data());

    TH1D* h_y_truth_mc = new TH1D("y_truth_mc", "Truth y (MC);y;Counts", n_y_bins, 0.0, 1.0);
    TH1D* h_y_reco_mc = new TH1D("y_reco_mc", "Reco y (MC);y;Counts", n_y_bins, 0.0, 1.0);
    TH1D* h_y_reco_pdata = new TH1D("y_reco_pdata", "Reco y (Pseudo-data);y;Counts", n_y_bins, 0.0, 1.0);
    TH1D* h_y_truth_pdata = new TH1D("y_truth_pdata", "Truth y (Pseudo-data);y;Counts", n_y_bins, 0.0, 1.0);

    TH1D* h_W2_truth_mc = new TH1D("W2_truth_mc", "Truth W^{2} (MC);W^{2} [GeV^{2}];Counts",
                                  w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_reco_mc = new TH1D("W2_reco_mc", "Reco W^{2} (MC, best DA#rightarrowEM blend);W^{2} [GeV^{2}];Counts",
                                 w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_reco_pdata = new TH1D("W2_reco_pdata", "Reco W^{2} (Pseudo-data, best DA#rightarrowEM blend);W^{2} [GeV^{2}];Counts",
                                    w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_truth_pdata = new TH1D("W2_truth_pdata", "Truth W^{2} (Pseudo-data);W^{2} [GeV^{2}];Counts",
                                     w2_bins.size()-1, w2_bins.data());

    TH1D* h_t_truth_mc_B0 = new TH1D("t_truth_mc_B0", "Truth |t| (MC, B0);|t| [GeV^{2}];Counts",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_t_reco_mc_B0 = new TH1D("t_reco_mc_B0", "Reco |t| (MC, B0);|t| [GeV^{2}];Counts",
                                   t_bins.size()-1, t_bins.data());
    TH1D* h_t_same_mc_B0 = new TH1D("t_same_mc_B0", "Generated and measured in same |t| bin (MC, B0);|t| [GeV^{2}];Counts",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_t_reco_pdata_B0 = new TH1D("t_reco_pdata_B0", "Reco |t| (Pseudo-data, B0);|t| [GeV^{2}];Counts",
                                      t_bins.size()-1, t_bins.data());
    TH1D* h_t_truth_pdata_B0 = new TH1D("t_truth_pdata_B0", "Truth |t| (Pseudo-data, B0);|t| [GeV^{2}];Counts",
                                       t_bins.size()-1, t_bins.data());

    TH1D* h_t_truth_mc_RP = new TH1D("t_truth_mc_RP", "Truth |t| (MC, RP);|t| [GeV^{2}];Counts",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_t_reco_mc_RP = new TH1D("t_reco_mc_RP", "Reco |t| (MC, RP);|t| [GeV^{2}];Counts",
                                   t_bins.size()-1, t_bins.data());
    TH1D* h_t_same_mc_RP = new TH1D("t_same_mc_RP", "Generated and measured in same |t| bin (MC, RP);|t| [GeV^{2}];Counts",
                                    t_bins.size()-1, t_bins.data());
    TH1D* h_t_reco_pdata_RP = new TH1D("t_reco_pdata_RP", "Reco |t| (Pseudo-data, RP);|t| [GeV^{2}];Counts",
                                      t_bins.size()-1, t_bins.data());
    TH1D* h_t_truth_pdata_RP = new TH1D("t_truth_pdata_RP", "Truth |t| (Pseudo-data, RP);|t| [GeV^{2}];Counts",
                                       t_bins.size()-1, t_bins.data());

    TH1D* h_beta_truth_mc_B0 = new TH1D("beta_truth_mc_B0", "Truth #beta (MC, B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_mc_B0 = new TH1D("beta_reco_mc_B0", "Reco #beta (MC, B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_pdata_B0 = new TH1D("beta_reco_pdata_B0", "Reco #beta (Pseudo-data, B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_truth_pdata_B0 = new TH1D("beta_truth_pdata_B0", "Truth #beta (Pseudo-data, B0);#beta;Counts", beta_bins.size()-1, beta_bins.data());

    TH1D* h_beta_truth_mc_RP = new TH1D("beta_truth_mc_RP", "Truth #beta (MC, RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_mc_RP = new TH1D("beta_reco_mc_RP", "Reco #beta (MC, RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_reco_pdata_RP = new TH1D("beta_reco_pdata_RP", "Reco #beta (Pseudo-data, RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());
    TH1D* h_beta_truth_pdata_RP = new TH1D("beta_truth_pdata_RP", "Truth #beta (Pseudo-data, RP);#beta;Counts", beta_bins.size()-1, beta_bins.data());

    TH1D* h_xpom_truth_mc_B0 = new TH1D("xpom_truth_mc_B0", "Truth x_{pom} (MC, B0);x_{pom};Counts",
                                       xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_mc_B0 = new TH1D("xpom_reco_mc_B0", "Reco x_{pom} (MC, B0);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_pdata_B0 = new TH1D("xpom_reco_pdata_B0", "Reco x_{pom} (Pseudo-data, B0);x_{pom};Counts",
                                         xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_truth_pdata_B0 = new TH1D("xpom_truth_pdata_B0", "Truth x_{pom} (Pseudo-data, B0);x_{pom};Counts",
                                          xpom_bins.size()-1, xpom_bins.data());

    TH1D* h_xpom_truth_mc_RP = new TH1D("xpom_truth_mc_RP", "Truth x_{pom} (MC, RP);x_{pom};Counts",
                                       xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_mc_RP = new TH1D("xpom_reco_mc_RP", "Reco x_{pom} (MC, RP);x_{pom};Counts",
                                      xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_reco_pdata_RP = new TH1D("xpom_reco_pdata_RP", "Reco x_{pom} (Pseudo-data, RP);x_{pom};Counts",
                                         xpom_bins.size()-1, xpom_bins.data());
    TH1D* h_xpom_truth_pdata_RP = new TH1D("xpom_truth_pdata_RP", "Truth x_{pom} (Pseudo-data, RP);x_{pom};Counts",
                                          xpom_bins.size()-1, xpom_bins.data());

    TH1D* h_theta_truth_mc_B0 = new TH1D("theta_truth_mc_B0", "Truth #theta (MC, B0);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_reco_mc_B0 = new TH1D("theta_reco_mc_B0", "Reco #theta (MC, B0);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_reco_pdata_B0 = new TH1D("theta_reco_pdata_B0", "Reco #theta (Pseudo-data, B0);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_truth_pdata_B0 = new TH1D("theta_truth_pdata_B0", "Truth #theta (Pseudo-data, B0);#theta [mrad];Counts", 100, 0.0, 25.0);

    TH1D* h_theta_truth_mc_RP = new TH1D("theta_truth_mc_RP", "Truth #theta (MC, RP);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_reco_mc_RP = new TH1D("theta_reco_mc_RP", "Reco #theta (MC, RP);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_reco_pdata_RP = new TH1D("theta_reco_pdata_RP", "Reco #theta (Pseudo-data, RP);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_truth_pdata_RP = new TH1D("theta_truth_pdata_RP", "Truth #theta (Pseudo-data, RP);#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_truth_pdata_all = new TH1D("theta_truth_pdata_all", "Truth #theta (Pseudo-data);#theta [mrad];Counts", 100, 0.0, 25.0);

    // Phase-space density plots (truth, unbinned-style)
    const int n_beta_bins_density = 120;
    const int n_xpom_bins_density = 120;
    std::vector<Double_t> xpom_bins_density = GetLogBins(1.0e-4, 0.3, n_xpom_bins_density);
    // keep a finer t binning for density plots (independent of YAML binning)
    std::vector<Double_t> t_bins_density = GetLogBins(1e-3, 2.0, 120);
    std::vector<Double_t> beta_bins_density;
    beta_bins_density.reserve(n_beta_bins_density + 1);
    for (int i = 0; i <= n_beta_bins_density; ++i) {
        beta_bins_density.push_back(static_cast<Double_t>(i) / n_beta_bins_density);
    }
    TH2D* h_beta_Q2_truth = new TH2D("beta_Q2_truth", "Event Density;#beta;Q^{2} [GeV^{2}]",
                                    n_beta_bins_density, 0.0, 1.0,
                                    q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_t_Q2_truth = new TH2D("t_Q2_truth", "Event Density;|t| [GeV^{2}];Q^{2} [GeV^{2}]",
                                 t_bins_density.size()-1, t_bins_density.data(),
                                 q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_xpom_Q2_truth = new TH2D("xpom_Q2_truth", "Event Density;x_{pom};Q^{2} [GeV^{2}]",
                                    xpom_bins_density.size()-1, xpom_bins_density.data(),
                                    q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_beta_t_truth = new TH2D("beta_t_truth", "Event Density;#beta;|t| [GeV^{2}]",
                                   n_beta_bins_density, 0.0, 1.0,
                                   t_bins_density.size()-1, t_bins_density.data());
    TH2D* h_xbj_t_truth = new TH2D("xbj_t_truth", "Event Density;x_{Bj};|t| [GeV^{2}]",
                                  x_bins_density.size()-1, x_bins_density.data(),
                                  t_bins_density.size()-1, t_bins_density.data());
    TH2D* h_xpom_t_truth = new TH2D("xpom_t_truth", "Event Density;x_{pom};|t| [GeV^{2}]",
                                   xpom_bins_density.size()-1, xpom_bins_density.data(),
                                   t_bins_density.size()-1, t_bins_density.data());
    TH2D* h_xpom_beta_truth = new TH2D("xpom_beta_truth", "Event Density;x_{pom};#beta",
                                      xpom_bins_density.size()-1, xpom_bins_density.data(),
                                      n_beta_bins_density, 0.0, 1.0);
    TH2D* h_xbj_beta_truth = new TH2D("xbj_beta_truth", "Event Density;x_{Bj};#beta",
                                     x_bins_density.size()-1, x_bins_density.data(),
                                     n_beta_bins_density, 0.0, 1.0);
    TH2D* h_xbj_xpom_truth = new TH2D("xbj_xpom_truth", "Event Density;x_{Bj};x_{pom}",
                                     x_bins_density.size()-1, x_bins_density.data(),
                                     xpom_bins_density.size()-1, xpom_bins_density.data());

    TH2D* h_beta_Q2_reco = new TH2D("beta_Q2_reco", "Event Density (Reco);#beta;Q^{2} [GeV^{2}]",
                                   n_beta_bins_density, 0.0, 1.0,
                                   q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_t_Q2_reco = new TH2D("t_Q2_reco", "Event Density (Reco);|t| [GeV^{2}];Q^{2} [GeV^{2}]",
                                t_bins_density.size()-1, t_bins_density.data(),
                                q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_xpom_Q2_reco = new TH2D("xpom_Q2_reco", "Event Density (Reco);x_{pom};Q^{2} [GeV^{2}]",
                                   xpom_bins_density.size()-1, xpom_bins_density.data(),
                                   q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_beta_t_reco = new TH2D("beta_t_reco", "Event Density (Reco);#beta;|t| [GeV^{2}]",
                                  n_beta_bins_density, 0.0, 1.0,
                                  t_bins_density.size()-1, t_bins_density.data());
    TH2D* h_xbj_t_reco = new TH2D("xbj_t_reco", "Event Density (Reco);x_{Bj};|t| [GeV^{2}]",
                                 x_bins_density.size()-1, x_bins_density.data(),
                                 t_bins_density.size()-1, t_bins_density.data());
    TH2D* h_xpom_t_reco = new TH2D("xpom_t_reco", "Event Density (Reco);x_{pom};|t| [GeV^{2}]",
                                  xpom_bins_density.size()-1, xpom_bins_density.data(),
                                  t_bins_density.size()-1, t_bins_density.data());
    TH2D* h_xpom_beta_reco = new TH2D("xpom_beta_reco", "Event Density (Reco);x_{pom};#beta",
                                     xpom_bins_density.size()-1, xpom_bins_density.data(),
                                     n_beta_bins_density, 0.0, 1.0);
    TH2D* h_xbj_beta_reco = new TH2D("xbj_beta_reco", "Event Density (Reco);x_{Bj};#beta",
                                    x_bins_density.size()-1, x_bins_density.data(),
                                    n_beta_bins_density, 0.0, 1.0);
    TH2D* h_xbj_xpom_reco = new TH2D("xbj_xpom_reco", "Event Density (Reco);x_{Bj};x_{pom}",
                                    x_bins_density.size()-1, x_bins_density.data(),
                                    xpom_bins_density.size()-1, xpom_bins_density.data());

    TH3D* h_phase3D_reco = new TH3D("phase3D_reco",
                                    "Reco Phase Space;Q^{2} [GeV^{2}];x_{pom};#beta",
                                    q2_bins_density.size()-1, q2_bins_density.data(),
                                    xpom_bins_density.size()-1, xpom_bins_density.data(),
                                    beta_bins_density.size()-1, beta_bins_density.data());

    TH3D* h_phase3D_reco_yaml = nullptr;
    h_phase3D_reco_yaml = new TH3D("phase3D_reco_yaml",
                                   "Reco Phase Space (YAML);Q^{2} [GeV^{2}];x_{pom};#beta",
                                   yaml_q2_edges.size()-1, yaml_q2_edges.data(),
                                   yaml_xpom_edges.size()-1, yaml_xpom_edges.data(),
                                   yaml_beta_edges.size()-1, yaml_beta_edges.data());
    if (!WriteBinsTSV("/data/bins_output.txt", yaml_bins)) {
        if (WriteBinsTSV("data/bins_output.txt", yaml_bins)) {
            std::cerr << "Warning: /data not writable; wrote data/bins_output.txt instead." << std::endl;
        } else {
            std::cerr << "Warning: failed to write bins_output.txt to /data and ./data." << std::endl;
        }
    }

    TH1D* h_t_res_B0 = new TH1D("t_res_B0", "B0 t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 51, -0.3, 0.3);
    TH1D* h_t_res_RP = new TH1D("t_res_RP", "RP t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 41, -0.5, 0.5);

    TH2D* h_t_RelRes_binned_B0 = new TH2D("t_RelRes_binned_B0",
                                         "|t|: truth vs rel. res (B0);|t|_{truth} [GeV^{2}];(#t_{reco}-#t_{truth})/#t_{truth}",
                                         t_bins.size()-1, t_bins.data(), 100, -2.0, 2.0);
    TH2D* h_t_RelRes_binned_RP = new TH2D("t_RelRes_binned_RP",
                                         "|t|: truth vs rel. res (RP);|t|_{truth} [GeV^{2}];(#t_{reco}-#t_{truth})/#t_{truth}",
                                         t_bins.size()-1, t_bins.data(), 100, -2.0, 2.0);

    TH2D* h_t_corr_B0 = new TH2D("t_corr_B0", "B0 |t| Response;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                 t_bins.size()-1, t_bins.data(), t_bins.size()-1, t_bins.data());
    TH2D* h_t_corr_RP = new TH2D("t_corr_RP", "RP |t| Response;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                 t_bins.size()-1, t_bins.data(), t_bins.size()-1, t_bins.data());

    TH2D* h_MX2_t_truth = new TH2D("MX2_t_truth", "M_{X}^{2} vs |t| (Truth);M_{X}^{2} [GeV^{2}];|t| [GeV^{2}]",
                                   mx2_bins.size() - 1, mx2_bins.data(),
                                   t_bins.size()-1, t_bins.data());
    TH2D* h_MX2_t_B0 = new TH2D("MX2_t_B0", "M_{X}^{2} vs |t| (B0 Reco);M_{X}^{2} [GeV^{2}];|t| [GeV^{2}]",
                                mx2_bins.size() - 1, mx2_bins.data(),
                                t_bins.size()-1, t_bins.data());
    TH2D* h_MX2_t_RP = new TH2D("MX2_t_RP", "M_{X}^{2} vs |t| (RP Reco);M_{X}^{2} [GeV^{2}];|t| [GeV^{2}]",
                                mx2_bins.size() - 1, mx2_bins.data(),
                                t_bins.size()-1, t_bins.data());

    TGraph* g_t_B0 = new TGraph(); g_t_B0->SetName("g_t_B0");
    TGraph* g_t_RP = new TGraph(); g_t_RP->SetName("g_t_RP");
    int n_g_t_B0 = 0, n_g_t_RP = 0;

    TGraph* g_xL_B0 = new TGraph(); g_xL_B0->SetName("g_xL_B0");
    TGraph* g_xL_RP = new TGraph(); g_xL_RP->SetName("g_xL_RP");
    int n_g_xL_B0 = 0, n_g_xL_RP = 0;

    TGraph* g_xpom_EM_B0 = new TGraph(); g_xpom_EM_B0->SetName("g_xpom_EM_B0");
    TGraph* g_xpom_EM_RP = new TGraph(); g_xpom_EM_RP->SetName("g_xpom_EM_RP");
    TGraph* g_xpom_DA_B0 = new TGraph(); g_xpom_DA_B0->SetName("g_xpom_DA_B0");
    TGraph* g_xpom_DA_RP = new TGraph(); g_xpom_DA_RP->SetName("g_xpom_DA_RP");
    TGraph* g_xpom_Sigma_B0 = new TGraph(); g_xpom_Sigma_B0->SetName("g_xpom_Sigma_B0");
    TGraph* g_xpom_Sigma_RP = new TGraph(); g_xpom_Sigma_RP->SetName("g_xpom_Sigma_RP");
    TGraph* g_beta_EM_B0 = new TGraph(); g_beta_EM_B0->SetName("g_beta_EM_B0");
    TGraph* g_beta_EM_RP = new TGraph(); g_beta_EM_RP->SetName("g_beta_EM_RP");
    TGraph* g_beta_DA_B0 = new TGraph(); g_beta_DA_B0->SetName("g_beta_DA_B0");
    TGraph* g_beta_DA_RP = new TGraph(); g_beta_DA_RP->SetName("g_beta_DA_RP");
    TGraph* g_beta_Sigma_B0 = new TGraph(); g_beta_Sigma_B0->SetName("g_beta_Sigma_B0");
    TGraph* g_beta_Sigma_RP = new TGraph(); g_beta_Sigma_RP->SetName("g_beta_Sigma_RP");
    TGraph* g_MX2 = new TGraph(); g_MX2->SetName("g_MX2");
    int n_g_xpom_em_b0 = 0, n_g_xpom_em_rp = 0;
    int n_g_xpom_da_b0 = 0, n_g_xpom_da_rp = 0;
    int n_g_xpom_sigma_b0 = 0, n_g_xpom_sigma_rp = 0;
    int n_g_beta_em_b0 = 0, n_g_beta_em_rp = 0;
    int n_g_beta_da_b0 = 0, n_g_beta_da_rp = 0;
    int n_g_beta_sigma_b0 = 0, n_g_beta_sigma_rp = 0;
    int n_g_MX2 = 0;

    //---------------------------------------------------------
    // DECLARE TTREEREADER AND BRANCHES TO USE
    //---------------------------------------------------------
    TTreeReader tree_reader(events);
    TTreeReaderArray<double>  mc_px_array         = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double>  mc_py_array         = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double>  mc_pz_array         = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array        = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int>    mc_genStatus_array   = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int>    mc_pdg_array         = {tree_reader, "MCParticles.PDG"};
    TTreeReaderArray<unsigned int> assoc_rec_id   = {tree_reader, "ReconstructedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> assoc_sim_id   = {tree_reader, "ReconstructedParticleAssociations.simID"};
    TTreeReaderArray<float>  re_px_array          = {tree_reader, "ReconstructedParticles.momentum.x"};
    TTreeReaderArray<float>  re_py_array          = {tree_reader, "ReconstructedParticles.momentum.y"};
    TTreeReaderArray<float>  re_pz_array          = {tree_reader, "ReconstructedParticles.momentum.z"};
    TTreeReaderArray<float>  re_energy_array      = {tree_reader, "ReconstructedParticles.energy"};
    TTreeReaderArray<int>    re_pdg_array         = {tree_reader, "ReconstructedParticles.PDG"};
    TTreeReaderArray<int>    electron_scat_index  = {tree_reader, "ScatteredElectronsTruth_objIdx.index"};

    const bool hasTSProtons = events->GetBranch("ReconstructedTruthSeededChargedParticles.momentum.x") &&
                              events->GetBranch("ReconstructedTruthSeededChargedParticles.momentum.y") &&
                              events->GetBranch("ReconstructedTruthSeededChargedParticles.momentum.z") &&
                              events->GetBranch("ReconstructedTruthSeededChargedParticleAssociations.recID") &&
                              events->GetBranch("ReconstructedTruthSeededChargedParticleAssociations.simID");
    if (!hasTSProtons) {
        std::cerr << "WARNING: Truth-seeded charged particle branches missing; B0 t reconstruction will be skipped." << std::endl;
    }

    const bool hasRPProtons = events->GetBranch("ForwardRomanPotRecParticles.momentum.x") &&
                              events->GetBranch("ForwardRomanPotRecParticles.momentum.y") &&
                              events->GetBranch("ForwardRomanPotRecParticles.momentum.z") &&
                              events->GetBranch("ForwardRomanPotRecParticles.mass") &&
                              events->GetBranch("ForwardRomanPotRecParticles.PDG");
    if (!hasRPProtons) {
        std::cerr << "WARNING: Roman Pot branches missing; RP t reconstruction will be skipped." << std::endl;
    }

    std::unique_ptr<TTreeReaderArray<unsigned int>> tsassoc_rec_id;
    std::unique_ptr<TTreeReaderArray<unsigned int>> tsassoc_sim_id;
    std::unique_ptr<TTreeReaderArray<float>> tsre_px_array;
    std::unique_ptr<TTreeReaderArray<float>> tsre_py_array;
    std::unique_ptr<TTreeReaderArray<float>> tsre_pz_array;
    if (hasTSProtons) {
        tsassoc_rec_id = std::make_unique<TTreeReaderArray<unsigned int>>(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID");
        tsassoc_sim_id = std::make_unique<TTreeReaderArray<unsigned int>>(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID");
        tsre_px_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
        tsre_py_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
        tsre_pz_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    }

    std::unique_ptr<TTreeReaderArray<float>> rp_px_array;
    std::unique_ptr<TTreeReaderArray<float>> rp_py_array;
    std::unique_ptr<TTreeReaderArray<float>> rp_pz_array;
    std::unique_ptr<TTreeReaderArray<float>> rp_mass_array;
    std::unique_ptr<TTreeReaderArray<int>> rp_pdg_array;
    if (hasRPProtons) {
        rp_px_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
        rp_py_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
        rp_pz_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ForwardRomanPotRecParticles.momentum.z");
        rp_mass_array = std::make_unique<TTreeReaderArray<float>>(tree_reader, "ForwardRomanPotRecParticles.mass");
        rp_pdg_array = std::make_unique<TTreeReaderArray<int>>(tree_reader, "ForwardRomanPotRecParticles.PDG");
    }

    auto requireBranch = [&](const char* name) {
        if (!events->GetBranch(name)) {
            std::cerr << "ERROR: Missing required branch " << name << std::endl;
            return false;
        }
        return true;
    };

    if (!requireBranch("InclusiveKinematicsElectron.Q2") ||
        !requireBranch("InclusiveKinematicsElectron.x") ||
        !requireBranch("InclusiveKinematicsElectron.y") ||
        !requireBranch("InclusiveKinematicsElectron.W") ||
        !requireBranch("InclusiveKinematicsDA.Q2") ||
        !requireBranch("InclusiveKinematicsDA.x") ||
        !requireBranch("InclusiveKinematicsDA.y") ||
        !requireBranch("InclusiveKinematicsDA.W")) {
        return 1;
    }

    const bool hasSigmaKin = events->GetBranch("InclusiveKinematicsSigma.Q2") != nullptr;
    const bool hasESigmaKin = events->GetBranch("InclusiveKinematicsESigma.Q2") != nullptr;
    const bool sigmaAvailable = hasSigmaKin || hasESigmaKin;
    if (hasSigmaKin) {
        std::cout << "Using InclusiveKinematicsSigma for Sigma method." << std::endl;
    } else if (hasESigmaKin) {
        std::cout << "Using InclusiveKinematicsESigma for Sigma method." << std::endl;
    } else {
        std::cerr << "WARNING: No Sigma/ESigma branches found; Sigma method will be skipped." << std::endl;
    }

    const bool hasTruthKin = events->GetBranch("InclusiveKinematicsTruth.Q2") &&
                             events->GetBranch("InclusiveKinematicsTruth.x") &&
                             events->GetBranch("InclusiveKinematicsTruth.y") &&
                             events->GetBranch("InclusiveKinematicsTruth.W");
    if (!hasTruthKin) {
        std::cerr << "WARNING: InclusiveKinematicsTruth not found; truth kinematics will be derived from MC." << std::endl;
    }

    TTreeReaderArray<float> kin_Q2_EM(tree_reader, "InclusiveKinematicsElectron.Q2");
    TTreeReaderArray<float> kin_x_EM(tree_reader, "InclusiveKinematicsElectron.x");
    TTreeReaderArray<float> kin_y_EM(tree_reader, "InclusiveKinematicsElectron.y");
    TTreeReaderArray<float> kin_W_EM(tree_reader, "InclusiveKinematicsElectron.W");

    TTreeReaderArray<float> kin_Q2_DA(tree_reader, "InclusiveKinematicsDA.Q2");
    TTreeReaderArray<float> kin_x_DA(tree_reader, "InclusiveKinematicsDA.x");
    TTreeReaderArray<float> kin_y_DA(tree_reader, "InclusiveKinematicsDA.y");
    TTreeReaderArray<float> kin_W_DA(tree_reader, "InclusiveKinematicsDA.W");

    std::unique_ptr<TTreeReaderArray<float>> kin_Q2_Sigma;
    std::unique_ptr<TTreeReaderArray<float>> kin_x_Sigma;
    std::unique_ptr<TTreeReaderArray<float>> kin_y_Sigma;
    std::unique_ptr<TTreeReaderArray<float>> kin_W_Sigma;
    if (hasSigmaKin) {
        kin_Q2_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsSigma.Q2");
        kin_x_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsSigma.x");
        kin_y_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsSigma.y");
        kin_W_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsSigma.W");
    } else if (hasESigmaKin) {
        kin_Q2_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsESigma.Q2");
        kin_x_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsESigma.x");
        kin_y_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsESigma.y");
        kin_W_Sigma = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsESigma.W");
    }

    std::unique_ptr<TTreeReaderArray<float>> kin_Q2_truth;
    std::unique_ptr<TTreeReaderArray<float>> kin_x_truth;
    std::unique_ptr<TTreeReaderArray<float>> kin_y_truth;
    std::unique_ptr<TTreeReaderArray<float>> kin_W_truth;
    if (hasTruthKin) {
        kin_Q2_truth = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsTruth.Q2");
        kin_x_truth = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsTruth.x");
        kin_y_truth = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsTruth.y");
        kin_W_truth = std::make_unique<TTreeReaderArray<float>>(tree_reader, "InclusiveKinematicsTruth.W");
    }

    //---------------------------------------------------------
    // FIND BEAM PARTICLES
    //---------------------------------------------------------
    BeamInfo beams;
    double eBeamGeV = 0.0;
    double pBeamGeV = 0.0;
    std::string firstMatchedBeamFile;
    std::string mismatchBeamFile;
    const bool parsedFromFilename = InferBeamEnergiesFromFileList(
        addedFiles, eBeamGeV, pBeamGeV, &firstMatchedBeamFile, &mismatchBeamFile
    );

    bool usingFilenameBeams = false;
    if (parsedFromFilename && eBeamGeV > beams.fMass_electron && pBeamGeV > beams.fMass_proton) {
        const double ePz = -std::sqrt(std::max(0.0, eBeamGeV * eBeamGeV -
                                                     beams.fMass_electron * beams.fMass_electron));
        const double pPz = std::sqrt(std::max(0.0, pBeamGeV * pBeamGeV -
                                                    beams.fMass_proton * beams.fMass_proton));
        beams.e_beam.SetCoordinates(0.0, 0.0, ePz, beams.fMass_electron);
        beams.p_beam.SetCoordinates(0.0, 0.0, pPz, beams.fMass_proton);
        usingFilenameBeams = true;
        std::cout << "Using beam energies from filename tag (" << eBeamGeV << "x"
                  << pBeamGeV << " GeV), first match: " << firstMatchedBeamFile << std::endl;
    }

    if (!usingFilenameBeams) {
        if (!mismatchBeamFile.empty()) {
            std::cerr << "WARNING: Inconsistent beam-energy tags in file names (first match: "
                      << firstMatchedBeamFile << ", mismatch: " << mismatchBeamFile
                      << "). Falling back to beam-particle scan." << std::endl;
        } else {
            std::cout << "Beam-energy tag not found in file names; falling back to beam-particle scan." << std::endl;
        }

        std::cout << "Finding beam particles..." << std::endl;
        P3MVector beame4_acc(0,0,0,0), beamp4_acc(0,0,0,0);

        while(tree_reader.Next()){
            for(int i = 0; i < mc_px_array.GetSize(); i++){
                if(mc_genStatus_array[i] != 4) continue;

                if(mc_pdg_array[i] == 2212){
                    P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], beams.fMass_proton);
                    beamp4_acc += p;
                }
                else if(mc_pdg_array[i] == 11){
                    P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], beams.fMass_electron);
                    beame4_acc += p;
                }
            }
        }

        auto nEntries = std::max<Long64_t>(1, events->GetEntries());
        beams.e_beam.SetCoordinates(
            beame4_acc.X()/nEntries,
            beame4_acc.Y()/nEntries,
            beame4_acc.Z()/nEntries,
            beams.fMass_electron
        );
        beams.p_beam.SetCoordinates(
            beamp4_acc.X()/nEntries,
            beamp4_acc.Y()/nEntries,
            beamp4_acc.Z()/nEntries,
            beams.fMass_proton
        );

        std::cout << "Found beam energies " << beams.e_beam.E() << "x" << beams.p_beam.E() << " GeV" << std::endl;
    }
    
    undoAfterburnAndCalc(beams.p_beam, beams.e_beam);
    
    tree_reader.Restart();

    //---------------------------------------------------------
    // RUN OVER TTREEREADER
    //---------------------------------------------------------
    Long64_t nentries = events->GetEntries();
    Long64_t nProcessed = 0;
    Long64_t nTruthValid = 0;
    Long64_t nEMValid = 0;
    Long64_t nDAValid = 0;
    Long64_t nSigmaValid = 0;
    Long64_t nTruthFromMC = 0;

    Long64_t nOutQ2Truth = 0, nOutXTruth = 0, nOutYTruth = 0;
    Long64_t nOutQ2EM = 0, nOutXEM = 0, nOutYEM = 0;
    Long64_t nOutQ2DA = 0, nOutXDA = 0, nOutYDA = 0;
    Long64_t nOutQ2Sigma = 0, nOutXSigma = 0, nOutYSigma = 0;
    Long64_t nPassTruthDIS = 0, nPassEMDIS = 0, nPassDADIS = 0, nPassSigmaDIS = 0;

    const CutflowConfig cutCfg{};
    auto fillSameBin1D = [](TH1D* refHist, TH1D* sameHist, double truthVal, double recoVal) {
        if (!refHist || !sameHist) return;
        const int nBins = refHist->GetNbinsX();
        const int truthBin = refHist->GetXaxis()->FindBin(truthVal);
        const int recoBin = refHist->GetXaxis()->FindBin(recoVal);
        if (truthBin < 1 || truthBin > nBins) return;
        if (recoBin < 1 || recoBin > nBins) return;
        if (truthBin == recoBin) sameHist->Fill(recoVal);
    };
    auto fillSameBin3D = [](TH3D* refHist,
                            TH3D* sameHist,
                            double q2Truth,
                            double betaTruth,
                            double xpomTruth,
                            double q2Reco,
                            double betaReco,
                            double xpomReco) {
        if (!refHist || !sameHist) return;
        const int nX = refHist->GetNbinsX();
        const int nY = refHist->GetNbinsY();
        const int nZ = refHist->GetNbinsZ();
        const int bxt = refHist->GetXaxis()->FindBin(q2Truth);
        const int byt = refHist->GetYaxis()->FindBin(betaTruth);
        const int bzt = refHist->GetZaxis()->FindBin(xpomTruth);
        const int bxr = refHist->GetXaxis()->FindBin(q2Reco);
        const int byr = refHist->GetYaxis()->FindBin(betaReco);
        const int bzr = refHist->GetZaxis()->FindBin(xpomReco);
        if (bxt < 1 || bxt > nX || byt < 1 || byt > nY || bzt < 1 || bzt > nZ) return;
        if (bxr < 1 || bxr > nX || byr < 1 || byr > nY || bzr < 1 || bzr > nZ) return;
        if (bxt == bxr && byt == byr && bzt == bzr) sameHist->Fill(q2Reco, betaReco, xpomReco);
    };
    for (Long64_t i = 0; i < nentries; i++) {
        // Update the counter every 1000 events
        if (i % 1000 == 0) {
        printf("\rProcessing event %lld of %lld; %.2f percent done.", i, nentries, 100.0*i/nentries);
        fflush(stdout);
        }
        if (!tree_reader.Next()) break;
        nProcessed++;
        const bool is_mc_subsample = ((i % 2) == 0);

        float electron_Q2_EM = GetArrayValue(kin_Q2_EM);
        float electron_x_EM = GetArrayValue(kin_x_EM);
        float electron_y_EM = GetArrayValue(kin_y_EM);
        float electron_W_EM = GetArrayValue(kin_W_EM);

        float electron_Q2_DA = GetArrayValue(kin_Q2_DA);
        float electron_x_DA = GetArrayValue(kin_x_DA);
        float electron_y_DA = GetArrayValue(kin_y_DA);
        float electron_W_DA = GetArrayValue(kin_W_DA);

        const bool hasSigmaMethod = (kin_Q2_Sigma != nullptr);
        float electron_Q2_Sigma = GetArrayValue(kin_Q2_Sigma.get());
        float electron_x_Sigma = GetArrayValue(kin_x_Sigma.get());
        float electron_y_Sigma = GetArrayValue(kin_y_Sigma.get());
        float electron_W_Sigma = GetArrayValue(kin_W_Sigma.get());

        float electron_Q2_truth = GetArrayValue(kin_Q2_truth.get());
        float electron_x_truth = GetArrayValue(kin_x_truth.get());
        float electron_y_truth = GetArrayValue(kin_y_truth.get());
        float electron_W_truth = GetArrayValue(kin_W_truth.get());

        int scat_mc_idx = (electron_scat_index.GetSize() > 0) ? electron_scat_index[0] : -1;
        if (scat_mc_idx < 0 || scat_mc_idx >= (int)mc_px_array.GetSize() || mc_pdg_array[scat_mc_idx] != 11) {
            double bestE = -1.0;
            int bestIdx = -1;
            for (int j = 0; j < mc_px_array.GetSize(); j++) {
                if (mc_genStatus_array[j] != 1) continue;
                if (mc_pdg_array[j] != 11) continue;
                const double px = mc_px_array[j];
                const double py = mc_py_array[j];
                const double pz = mc_pz_array[j];
                const double m = mc_mass_array[j];
                const double E = std::sqrt(px*px + py*py + pz*pz + m*m);
                if (E > bestE) {
                    bestE = E;
                    bestIdx = j;
                }
            }
            scat_mc_idx = bestIdx;
        }

        const bool truth_from_arrays = hasTruthKin && kin_Q2_truth && kin_x_truth && kin_y_truth && kin_W_truth &&
                                       kin_Q2_truth->GetSize() > 0 &&
                                       kin_x_truth->GetSize() > 0 &&
                                       kin_y_truth->GetSize() > 0 &&
                                       kin_W_truth->GetSize() > 0;
        if (!truth_from_arrays) {
            float Q2_calc = -999.0f, x_calc = -999.0f, y_calc = -999.0f, W_calc = -999.0f;
            if (ComputeTruthKinematics(beams, mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
                                       scat_mc_idx, Q2_calc, x_calc, y_calc, W_calc)) {
                electron_Q2_truth = Q2_calc;
                electron_x_truth = x_calc;
                electron_y_truth = y_calc;
                electron_W_truth = W_calc;
                nTruthFromMC++;
            }
        }

        const auto isValidQ2 = [](float v) { return std::isfinite(v) && v > 0.0f; };
        const auto isValidX = [](float v) { return std::isfinite(v) && v > 0.0f && v < 1.0f; };
        const auto isValidY = [](float v) { return std::isfinite(v) && v > 0.0f && v < 1.0f; };

        const bool rawTruthQ2 = isValidQ2(electron_Q2_truth);
        const bool rawTruthX = isValidX(electron_x_truth);
        const bool rawTruthY = isValidY(electron_y_truth);
        if (!rawTruthQ2) nOutQ2Truth++;
        if (!rawTruthX) nOutXTruth++;
        if (!rawTruthY) nOutYTruth++;

        const bool rawEMQ2 = isValidQ2(electron_Q2_EM);
        const bool rawEMX = isValidX(electron_x_EM);
        const bool rawEMY = isValidY(electron_y_EM);
        if (!rawEMQ2) nOutQ2EM++;
        if (!rawEMX) nOutXEM++;
        if (!rawEMY) nOutYEM++;

        const bool rawDAQ2 = isValidQ2(electron_Q2_DA);
        const bool rawDAX = isValidX(electron_x_DA);
        const bool rawDAY = isValidY(electron_y_DA);
        if (!rawDAQ2) nOutQ2DA++;
        if (!rawDAX) nOutXDA++;
        if (!rawDAY) nOutYDA++;

        const bool rawSigmaQ2 = isValidQ2(electron_Q2_Sigma);
        const bool rawSigmaX = isValidX(electron_x_Sigma);
        const bool rawSigmaY = isValidY(electron_y_Sigma);
        if (hasSigmaMethod) {
            if (!rawSigmaQ2) nOutQ2Sigma++;
            if (!rawSigmaX) nOutXSigma++;
            if (!rawSigmaY) nOutYSigma++;
        }

        // Reco/truth E-pz (currently matched-particle based), used by configurable DIS cuts.
        double sumEPz_truth_matched = 0.0;
        double sumEPz_reco_matched = 0.0;
        CalculateSumEPz_Matched(
            re_px_array, re_py_array, re_pz_array, re_energy_array,
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            assoc_rec_id, assoc_sim_id,
            sumEPz_truth_matched, sumEPz_reco_matched
        );
        h_EPz_truth->Fill(sumEPz_truth_matched);
        h_EPz->Fill(sumEPz_reco_matched);
        h_EPz_2D->Fill(sumEPz_truth_matched, sumEPz_reco_matched);

        const bool passTruthDIS = PassDISKinematicCuts(
            cutCfg, rawTruthQ2, rawTruthY, electron_Q2_truth, electron_y_truth, sumEPz_truth_matched, true);
        const bool passEMDIS = PassDISKinematicCuts(
            cutCfg, rawEMQ2, rawEMY, electron_Q2_EM, electron_y_EM, sumEPz_reco_matched, false);
        const bool passDADIS = PassDISKinematicCuts(
            cutCfg, rawDAQ2, rawDAY, electron_Q2_DA, electron_y_DA, sumEPz_reco_matched, false);
        const bool passSigmaDIS = hasSigmaMethod ? PassDISKinematicCuts(
            cutCfg, rawSigmaQ2, rawSigmaY, electron_Q2_Sigma, electron_y_Sigma, sumEPz_reco_matched, false) : false;

        if (passTruthDIS) nPassTruthDIS++;
        if (passEMDIS) nPassEMDIS++;
        if (passDADIS) nPassDADIS++;
        if (hasSigmaMethod && passSigmaDIS) nPassSigmaDIS++;

        // Centralized analysis-valid flags used everywhere below.
        const bool validTruthQ2 = rawTruthQ2 && passTruthDIS;
        const bool validTruthX = rawTruthX && passTruthDIS;
        const bool validTruthY = rawTruthY && passTruthDIS;

        const bool validEMQ2 = rawEMQ2 && passEMDIS;
        const bool validEMX = rawEMX && passEMDIS;
        const bool validEMY = rawEMY && passEMDIS;

        const bool validDAQ2 = rawDAQ2 && passDADIS;
        const bool validDAX = rawDAX && passDADIS;
        const bool validDAY = rawDAY && passDADIS;

        const bool validSigmaQ2 = rawSigmaQ2 && passSigmaDIS;
        const bool validSigmaX = rawSigmaX && passSigmaDIS;
        const bool validSigmaY = rawSigmaY && passSigmaDIS;

        const bool passRecoDIS = cutCfg.use_em_for_dis_selection ? passEMDIS : (passEMDIS || passDADIS || passSigmaDIS);

        if (validTruthQ2 && validTruthX && validTruthY) nTruthValid++;
        if (validEMQ2 && validEMX && validEMY) nEMValid++;
        if (validDAQ2 && validDAX && validDAY) nDAValid++;
        if (hasSigmaMethod && validSigmaQ2 && validSigmaX && validSigmaY) nSigmaValid++;

        if (validTruthX) h_x_truth->Fill(electron_x_truth);
        if (validEMX) h_x_EM->Fill(electron_x_EM);
        if (validDAX) h_x_DA->Fill(electron_x_DA);
        if (hasSigmaMethod && validSigmaX) h_x_Sigma->Fill(electron_x_Sigma);

        if (validTruthY) h_y_truth->Fill(electron_y_truth);
        if (validEMY) h_y_EM->Fill(electron_y_EM);
        if (validDAY) h_y_DA->Fill(electron_y_DA);
        if (hasSigmaMethod && validSigmaY) h_y_Sigma->Fill(electron_y_Sigma);

        // Relative resolutions vs truth (Set A only)
        if (is_mc_subsample) {
        if (validTruthX && validEMX) {
            h_RelRes_x_EM->Fill((electron_x_EM - electron_x_truth) / electron_x_truth);
            h_RelRes_x_binned_EM->Fill(electron_x_truth, (electron_x_EM - electron_x_truth) / electron_x_truth);
        }
        if (validTruthX && validDAX) {
            h_RelRes_x_DA->Fill((electron_x_DA - electron_x_truth) / electron_x_truth);
            h_RelRes_x_binned_DA->Fill(electron_x_truth, (electron_x_DA - electron_x_truth) / electron_x_truth);
        }
        if (hasSigmaMethod && validTruthX && validSigmaX) {
            h_RelRes_x_Sigma->Fill((electron_x_Sigma - electron_x_truth) / electron_x_truth);
            h_RelRes_x_binned_Sigma->Fill(electron_x_truth, (electron_x_Sigma - electron_x_truth) / electron_x_truth);
        }
        if (validTruthY && validEMY) {
            h_RelRes_y_EM->Fill((electron_y_EM - electron_y_truth) / electron_y_truth);
            h_RelRes_y_binned_EM->Fill(electron_y_truth, (electron_y_EM - electron_y_truth) / electron_y_truth);
        }
        if (validTruthY && validDAY) {
            h_RelRes_y_DA->Fill((electron_y_DA - electron_y_truth) / electron_y_truth);
            h_RelRes_y_binned_DA->Fill(electron_y_truth, (electron_y_DA - electron_y_truth) / electron_y_truth);
        }
        if (hasSigmaMethod && validTruthY && validSigmaY) {
            h_RelRes_y_Sigma->Fill((electron_y_Sigma - electron_y_truth) / electron_y_truth);
            h_RelRes_y_binned_Sigma->Fill(electron_y_truth, (electron_y_Sigma - electron_y_truth) / electron_y_truth);
        }

        if (validTruthX && validEMX) h_Corr_x_EM->Fill(electron_x_truth, electron_x_EM);
        if (validTruthX && validDAX) h_Corr_x_DA->Fill(electron_x_truth, electron_x_DA);
        if (hasSigmaMethod && validTruthX && validSigmaX) h_Corr_x_Sigma->Fill(electron_x_truth, electron_x_Sigma);
        if (validTruthY && validEMY) h_Corr_y_EM->Fill(electron_y_truth, electron_y_EM);
        if (validTruthY && validDAY) h_Corr_y_DA->Fill(electron_y_truth, electron_y_DA);
        if (hasSigmaMethod && validTruthY && validSigmaY) h_Corr_y_Sigma->Fill(electron_y_truth, electron_y_Sigma);

        // Correlations
        if (validTruthX && validEMX) g_x_EM->SetPoint(n_g_x_EM++, electron_x_truth, electron_x_EM);
        if (validTruthX && validDAX) g_x_DA->SetPoint(n_g_x_DA++, electron_x_truth, electron_x_DA);
        if (hasSigmaMethod && validTruthX && validSigmaX) g_x_Sigma->SetPoint(n_g_x_Sigma++, electron_x_truth, electron_x_Sigma);

        if (validTruthY && validEMY) g_y_EM->SetPoint(n_g_y_EM++, electron_y_truth, electron_y_EM);
        if (validTruthY && validDAY) g_y_DA->SetPoint(n_g_y_DA++, electron_y_truth, electron_y_DA);
        if (hasSigmaMethod && validTruthY && validSigmaY) g_y_Sigma->SetPoint(n_g_y_Sigma++, electron_y_truth, electron_y_Sigma);
        } // end is_mc_subsample gate for x/y resolution

        if (validEMQ2) h_Q2_EM->Fill(electron_Q2_EM);
        if (validDAQ2) h_Q2_DA->Fill(electron_Q2_DA);
        if (hasSigmaMethod && validSigmaQ2) h_Q2_Sigma->Fill(electron_Q2_Sigma);

        // Also fill the output tree for x and y
        if (validTruthX && validTruthY) {
            out_Q2_truth = electron_Q2_truth;
            out_Q2_EM    = electron_Q2_EM;
            out_Q2_DA    = electron_Q2_DA;
            out_Q2_Sigma = electron_Q2_Sigma;
            out_x_truth  = electron_x_truth;
            out_x_EM     = electron_x_EM;
            out_x_DA     = electron_x_DA;
            out_x_Sigma = electron_x_Sigma;
            out_y_truth  = electron_y_truth;
            out_y_EM     = electron_y_EM;
            out_y_DA     = electron_y_DA;
            out_y_Sigma = electron_y_Sigma;
            tree->Fill();
        }
        
        // Fill histograms for truth
        if (validTruthQ2) h_Q2_truth->Fill(electron_Q2_truth);
        if (validTruthQ2 && validTruthX) {
            h_xQ2_truth->Fill(electron_x_truth, electron_Q2_truth);
        }
        if (validEMQ2 && validEMX) {
            h_xQ2_reco->Fill(electron_x_EM, electron_Q2_EM);
        }

        // Q2/x/y resolution, correlations, and profile maps (Set A only)
        if (is_mc_subsample) {
        if (validTruthQ2 && validEMQ2) {
            h_RelRes_Q2_EM->Fill((electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth);
            h_RelRes_Q2_binned_EM->Fill(electron_Q2_truth, (electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth);
            h_Corr_Q2_EM->Fill(electron_Q2_truth, electron_Q2_EM);
        }
        if (validTruthQ2 && validDAQ2) {
            h_RelRes_Q2_DA->Fill((electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth);
            h_RelRes_Q2_binned_DA->Fill(electron_Q2_truth, (electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth);
            h_Corr_Q2_DA->Fill(electron_Q2_truth, electron_Q2_DA);
        }
        if (hasSigmaMethod && validTruthQ2 && validSigmaQ2) {
            h_RelRes_Q2_Sigma->Fill((electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth);
            h_RelRes_Q2_binned_Sigma->Fill(electron_Q2_truth, (electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth);
            h_Corr_Q2_Sigma->Fill(electron_Q2_truth, electron_Q2_Sigma);
        }

        // Fill Q2 resolution profile histograms - binned in (x, Q2)
        if (validTruthQ2 && validTruthX) {
            if (validEMQ2) {
                h_Q2_RelRes_vs_xy_EM->Fill(electron_x_truth, electron_Q2_truth,
                                           (electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth);
            }
            if (validDAQ2) {
                h_Q2_RelRes_vs_xy_DA->Fill(electron_x_truth, electron_Q2_truth,
                                           (electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth);
            }
            if (hasSigmaMethod && validSigmaQ2) {
                h_Q2_RelRes_vs_xy_Sigma->Fill(electron_x_truth, electron_Q2_truth,
                                              (electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth);
            }

            // Fill x resolution profile histograms - binned in (x, Q2)
            if (validEMX) {
                h_x_RelRes_vs_xQ2_EM->Fill(electron_x_truth, electron_Q2_truth,
                                           (electron_x_EM - electron_x_truth) / electron_x_truth);
            }
            if (validDAX) {
                h_x_RelRes_vs_xQ2_DA->Fill(electron_x_truth, electron_Q2_truth,
                                           (electron_x_DA - electron_x_truth) / electron_x_truth);
            }
            if (hasSigmaMethod && validSigmaX) {
                h_x_RelRes_vs_xQ2_Sigma->Fill(electron_x_truth, electron_Q2_truth,
                                              (electron_x_Sigma - electron_x_truth) / electron_x_truth);
            }
        }

        // Fill y resolution profile histograms - binned in (x, Q2)
        if (validTruthQ2 && validTruthX && validTruthY) {
            if (validEMY) {
                h_y_RelRes_vs_xQ2_EM->Fill(electron_x_truth, electron_Q2_truth,
                                           (electron_y_EM - electron_y_truth) / electron_y_truth);
            }
            if (validDAY) {
                h_y_RelRes_vs_xQ2_DA->Fill(electron_x_truth, electron_Q2_truth,
                                           (electron_y_DA - electron_y_truth) / electron_y_truth);
            }
            if (hasSigmaMethod && validSigmaY) {
                h_y_RelRes_vs_xQ2_Sigma->Fill(electron_x_truth, electron_Q2_truth,
                                              (electron_y_Sigma - electron_y_truth) / electron_y_truth);
            }
        }

        // Fill correlation histograms
        if (validTruthQ2 && validEMQ2) g_Q2_EM->SetPoint(n_g_Q2_EM++, electron_Q2_truth, electron_Q2_EM);
        if (validTruthQ2 && validDAQ2) g_Q2_DA->SetPoint(n_g_Q2_DA++, electron_Q2_truth, electron_Q2_DA);
        if (hasSigmaMethod && validTruthQ2 && validSigmaQ2) g_Q2_Sigma->SetPoint(n_g_Q2_Sigma++, electron_Q2_truth, electron_Q2_Sigma);
        } // end is_mc_subsample gate for Q2/x/y resolution

        //=================================================================
        // NEW INCLUSIVE DIS FILLS
        //=================================================================

        // Fill W² histograms
        const bool validWTruth = std::isfinite(electron_W_truth) && electron_W_truth > 0.0f && passTruthDIS;
        const bool validWEM = std::isfinite(electron_W_EM) && electron_W_EM > 0.0f && passEMDIS;
        const bool validWDA = std::isfinite(electron_W_DA) && electron_W_DA > 0.0f && passDADIS;
        const bool validWSigma = std::isfinite(electron_W_Sigma) && electron_W_Sigma > 0.0f && passSigmaDIS;

        const double W2_truth = validWTruth ? electron_W_truth * electron_W_truth : -1.0;
        const double W2_EM = validWEM ? electron_W_EM * electron_W_EM : -1.0;
        const double W2_DA = validWDA ? electron_W_DA * electron_W_DA : -1.0;
        const double W2_Sigma = validWSigma ? electron_W_Sigma * electron_W_Sigma : -1.0;
        double W2_best = -1.0;
        double alpha_W2_best = 0.0;
        double W2_proxy = -1.0;
        const bool validWBest = ComputeW2Best(W2_EM, W2_DA, W2_best, alpha_W2_best, W2_proxy);
        (void)alpha_W2_best;
        (void)W2_proxy;

        if (validWEM) h_W2_EM->Fill(W2_EM);
        if (validWDA) h_W2_DA->Fill(W2_DA);
        if (validWBest) h_W2_Best->Fill(W2_best);
        if (hasSigmaMethod && validWSigma) h_W2_Sigma->Fill(W2_Sigma);
        if (validWTruth) h_W2_truth->Fill(W2_truth);

        // W2 correlations, resolution, and response (Set A only)
        if (is_mc_subsample) {
        if (validWTruth && validWEM) g_W2_EM->SetPoint(n_g_W2_EM++, W2_truth, W2_EM);
        if (validWTruth && validWDA) g_W2_DA->SetPoint(n_g_W2_DA++, W2_truth, W2_DA);
        if (validWTruth && validWBest) g_W2_Best->SetPoint(n_g_W2_Best++, W2_truth, W2_best);
        if (hasSigmaMethod && validWTruth && validWSigma) g_W2_Sigma->SetPoint(n_g_W2_Sigma++, W2_truth, W2_Sigma);

        if (validWTruth && validWEM && W2_truth > 0.0) {
            const double rel = (W2_EM - W2_truth) / W2_truth;
            h_RelRes_W2_EM->Fill(rel);
            h_RelRes_W2_binned_EM->Fill(W2_truth, rel);
            h_Response_W2_EM->Fill(W2_truth, W2_EM);
        }
        if (validWTruth && validWDA && W2_truth > 0.0) {
            const double rel = (W2_DA - W2_truth) / W2_truth;
            h_RelRes_W2_DA->Fill(rel);
            h_RelRes_W2_binned_DA->Fill(W2_truth, rel);
            h_Response_W2_DA->Fill(W2_truth, W2_DA);
        }
        if (validWTruth && validWBest && W2_truth > 0.0) {
            const double rel = (W2_best - W2_truth) / W2_truth;
            h_RelRes_W2_Best->Fill(rel);
            h_RelRes_W2_binned_Best->Fill(W2_truth, rel);
            h_Response_W2_Best->Fill(W2_truth, W2_best);
        }
        if (hasSigmaMethod && validWTruth && validWSigma && W2_truth > 0.0) {
            const double rel = (W2_Sigma - W2_truth) / W2_truth;
            h_RelRes_W2_Sigma->Fill(rel);
            h_RelRes_W2_binned_Sigma->Fill(W2_truth, rel);
            h_Response_W2_Sigma->Fill(W2_truth, W2_Sigma);
        }
        } // end is_mc_subsample gate for W2 resolution

        const double m_p_sq = MASS_PROTON * MASS_PROTON;
        const double xpom_denominator_EM    = (validEMQ2 && validWEM)       ? (electron_Q2_EM    + W2_EM    - m_p_sq) : -1.0;
        const double xpom_denominator_truth = (validTruthQ2 && validWTruth) ? (electron_Q2_truth + W2_truth - m_p_sq) : -1.0;
        const double xpom_denominator_DA    = (validDAQ2 && validWDA)       ? (electron_Q2_DA    + W2_DA    - m_p_sq) : -1.0;
        const double xpom_denominator_Sigma = (hasSigmaMethod && validSigmaQ2 && validWSigma) ? (electron_Q2_Sigma + W2_Sigma - m_p_sq): -1.0;
        const double xpom_denominator_W2Best = (validEMQ2 && validWBest) ? (electron_Q2_EM + W2_best - m_p_sq) : -1.0;

        // eta_max (matched reco/truth particles)
        double eta_max_reco = -1.0e9;
        double eta_max_truth = -1.0e9;
        for (unsigned int j = 0; j < re_energy_array.GetSize(); ++j) {
            const P3MVector particle_reco(re_px_array[j], re_py_array[j], re_pz_array[j], 0.0);
            eta_max_reco = std::max(eta_max_reco, particle_reco.Eta());

            int mc_idx = -1;
            for (unsigned int k = 0; k < assoc_rec_id.GetSize(); ++k) {
                if (assoc_rec_id[k] == j) {
                    mc_idx = static_cast<int>(assoc_sim_id[k]);
                    break;
                }
            }
            if (mc_idx >= 0 && mc_idx < static_cast<int>(mc_px_array.GetSize())) {
                P3MVector particle_mc(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
                eta_max_truth = std::max(eta_max_truth, particle_mc.Eta());
            }
        }
        if (eta_max_reco > -1.0e8) h_eta_max->Fill(eta_max_reco);
        if (eta_max_truth > -1.0e8) h_eta_max_truth->Fill(eta_max_truth);

        // Response matrices (Electron method) (Set A only)
        if (is_mc_subsample) {
        if (validTruthQ2 && validEMQ2) h_Response_Q2->Fill(electron_Q2_truth, electron_Q2_EM);
        if (validTruthX && validEMX) h_Response_x->Fill(electron_x_truth, electron_x_EM);
        if (validTruthY && validEMY) h_Response_y->Fill(electron_y_truth, electron_y_EM);
        } // end is_mc_subsample gate for response matrices

        // Scattered electron kinematics
        int scat_reco_idx = -1;
        if (scat_mc_idx >= 0) {
            for (unsigned int j = 0; j < assoc_rec_id.GetSize(); j++) {
                if (assoc_sim_id[j] == (unsigned int)scat_mc_idx) {
                    scat_reco_idx = assoc_rec_id[j];
                    break;
                }
            }
        }

        if (scat_reco_idx < 0 && scat_mc_idx >= 0 &&
            scat_mc_idx < (int)re_px_array.GetSize() && re_pdg_array[scat_mc_idx] == 11) {
            scat_reco_idx = scat_mc_idx;
        }

        bool hasRecoElectron = false;
        double reco_E = 0.0, reco_phi = 0.0, reco_pT = 0.0;
        if (scat_reco_idx >= 0 && scat_reco_idx < (int)re_px_array.GetSize()) {
            const double px = re_px_array[scat_reco_idx];
            const double py = re_py_array[scat_reco_idx];
            const double pz = re_pz_array[scat_reco_idx];
            const double E  = re_energy_array[scat_reco_idx];

            reco_phi = TMath::ATan2(py, px);
            reco_pT = TMath::Sqrt(px*px + py*py);
            reco_E = E;

            if (passRecoDIS) {
                h_Ep_e->Fill(E);
                h_phi_e->Fill(reco_phi);
                h_pT_e->Fill(reco_pT);
                hasRecoElectron = true;
            }
        }

        bool hasTruthElectron = false;
        double truth_E = 0.0, truth_phi = 0.0, truth_pT = 0.0;
        if (scat_mc_idx >= 0 && scat_mc_idx < (int)mc_px_array.GetSize()) {
            const double mc_px = mc_px_array[scat_mc_idx];
            const double mc_py = mc_py_array[scat_mc_idx];
            const double mc_pz = mc_pz_array[scat_mc_idx];
            const double mc_m  = mc_mass_array[scat_mc_idx];
            const double mc_E  = TMath::Sqrt(mc_px*mc_px + mc_py*mc_py + mc_pz*mc_pz + mc_m*mc_m);
            truth_phi   = TMath::ATan2(mc_py, mc_px);
            truth_pT    = TMath::Sqrt(mc_px*mc_px + mc_py*mc_py);
            truth_E = mc_E;

            if (passTruthDIS) {
                h_Ep_e_truth->Fill(mc_E);
                h_phi_e_truth->Fill(truth_phi);
                h_pT_e_truth->Fill(truth_pT);
                hasTruthElectron = true;
            }
        }

        if (hasTruthElectron && hasRecoElectron) {
            g_Ep_e->SetPoint(n_g_Ep_e++, truth_E, reco_E);
            g_phi_e->SetPoint(n_g_phi_e++, truth_phi, reco_phi);
            g_pT_e->SetPoint(n_g_pT_e++, truth_pT, reco_pT);
        }

        //=================================================================
        // DIFFRACTIVE: M_X^2 (hadronic system) excluding scattered e- and leading proton
        //=================================================================
        int scat_reco_idx_mx2 = scat_reco_idx;
        if (scat_reco_idx_mx2 < 0) {
            double bestE = -1.0;
            int bestIdx = -1;
            for (int j = 0; j < re_energy_array.GetSize(); j++) {
                if (re_pdg_array[j] != 11) continue;
                if (re_energy_array[j] > bestE) {
                    bestE = re_energy_array[j];
                    bestIdx = j;
                }
            }
            scat_reco_idx_mx2 = bestIdx;
        }

        int lead_reco_proton_idx = -1;
        double lead_reco_pz = -1.0e12;
        for (int j = 0; j < re_px_array.GetSize(); j++) {
            if (re_pdg_array[j] != 2212) continue;
            if (re_pz_array[j] > lead_reco_pz) {
                lead_reco_pz = re_pz_array[j];
                lead_reco_proton_idx = j;
            }
        }

        P3EVector total_hadrons_reco(0.0, 0.0, 0.0, 0.0);
        for (int j = 0; j < re_energy_array.GetSize(); j++) {
            if (j == scat_reco_idx_mx2) continue;
            if (j == lead_reco_proton_idx) continue;
            P3EVector particle(re_px_array[j], re_py_array[j], re_pz_array[j], re_energy_array[j]);
            total_hadrons_reco += particle;
        }
        double MX2_reco = total_hadrons_reco.M2();

        P3EVector total_hadrons_truth(0.0, 0.0, 0.0, 0.0);
        for (unsigned int j = 0; j < re_energy_array.GetSize(); j++) {
            if ((int)j == scat_reco_idx_mx2) continue;
            if ((int)j == lead_reco_proton_idx) continue;

            int mc_idx = -1;
            for (unsigned int k = 0; k < assoc_rec_id.GetSize(); k++) {
                if (assoc_rec_id[k] == j) {
                    mc_idx = assoc_sim_id[k];
                    break;
                }
            }
            if (mc_idx < 0 || mc_idx >= (int)mc_px_array.GetSize()) continue;
            if (mc_genStatus_array[mc_idx] != 1) continue;
            if (mc_idx == scat_mc_idx) continue;
            const int pdg = mc_pdg_array[mc_idx];
            if (pdg == 12 || pdg == -12 || pdg == 14 || pdg == -14 || pdg == 16 || pdg == -16) continue;

            const double px = mc_px_array[mc_idx];
            const double py = mc_py_array[mc_idx];
            const double pz = mc_pz_array[mc_idx];
            const double m = mc_mass_array[mc_idx];
            const double E = TMath::Sqrt(px*px + py*py + pz*pz + m*m);
            P3EVector particle(px, py, pz, E);
            total_hadrons_truth += particle;
        }
        double MX2_truth = total_hadrons_truth.M2();

        const bool valid_MX2_reco = std::isfinite(MX2_reco) && MX2_reco >= 0.0;
        const bool valid_MX2_truth = std::isfinite(MX2_truth) && MX2_truth >= 0.0;
        if (valid_MX2_reco) h_MX2_reco->Fill(MX2_reco);
        if (valid_MX2_truth) h_MX2_truth->Fill(MX2_truth);
        // MX2 correlations and resolution (Set A only)
        if (is_mc_subsample) {
        if (valid_MX2_reco && valid_MX2_truth) h_MX2_corr->Fill(MX2_truth, MX2_reco);
        if (valid_MX2_reco && valid_MX2_truth) g_MX2->SetPoint(n_g_MX2++, MX2_truth, MX2_reco);
        if (valid_MX2_reco && valid_MX2_truth && MX2_truth > 1e-9) {
            const double rel_mx2 = (MX2_reco - MX2_truth) / MX2_truth;
            h_MX2_RelRes->Fill(rel_mx2);
            h_MX2_RelRes_binned->Fill(MX2_truth, rel_mx2);
            if (validTruthQ2) {
                h_MX2_RelRes_vs_MX2Q2->Fill(MX2_truth, electron_Q2_truth, rel_mx2);
            }
        }
        } // end is_mc_subsample gate for MX2 resolution

        //=================================================================
        // DIFFRACTIVE: Mandelstam t analysis
        //=================================================================
        bool xpom_truth_b0_filled = false;
        bool xpom_truth_rp_filled = false;
        bool xpom_reco_em_b0_filled = false;
        bool xpom_reco_w2best_b0_filled = false;
        bool xpom_reco_da_b0_filled = false;
        bool xpom_reco_sigma_b0_filled = false;
        bool xpom_reco_em_rp_filled = false;
        bool xpom_reco_w2best_rp_filled = false;
        bool xpom_reco_da_rp_filled = false;
        bool xpom_reco_sigma_rp_filled = false;
        bool xpom_reco_em_all_filled = false;
        bool xpom_reco_da_all_filled = false;
        bool xpom_reco_sigma_all_filled = false;
        bool has_truth_lead = false;
        double lead_truth_pz = -1.0;
        P3MVector lead_truth(0.0, 0.0, 0.0, 0.0);
        bool has_reco_proton = false;
        double t_reco_best = -1.0;
        std::vector<P3MVector> truth_protons;
        truth_protons.reserve(4);
        for (int j = 0; j < mc_px_array.GetSize(); j++) {
            if (mc_genStatus_array[j] != 1 || mc_pdg_array[j] != 2212) continue;
            P3MVector p(mc_px_array[j], mc_py_array[j], mc_pz_array[j], mc_mass_array[j]);
            undoAfterburn(p);
            truth_protons.push_back(p);
            if (p.Pz() > lead_truth_pz) {
                lead_truth_pz = p.Pz();
                lead_truth = p;
                has_truth_lead = true;
            }

            const double t_val = TMath::Abs(CalcT(beams.p_beam, p));
            const double theta_mrad = p.Theta() * 1000.0;
            if (passTruthDIS) {
                if (std::isfinite(t_val)) {
                    h_t_MC->Fill(t_val);
                    h_dsigma_dt_MC->Fill(t_val);
                    if (valid_MX2_truth) {
                        h_MX2_t_truth->Fill(MX2_truth, t_val);
                    }
                }
                h_theta_MC->Fill(theta_mrad);
                if (!is_mc_subsample) {
                    h_theta_truth_pdata_all->Fill(theta_mrad);
                }

                const double xL_pz = p.Pz() / beams.p_beam.Pz();
                if (std::isfinite(xL_pz)) {
                    h_xL_MC->Fill(xL_pz);
                    const double xpom_from_xL = 1.0 - xL_pz;
                    if (std::isfinite(xpom_from_xL) && xpom_from_xL > 0.0) {
                        h_xpom_MC->Fill(xpom_from_xL);
                    }
                }
                if (xpom_denominator_truth > 0.0 && valid_MX2_truth && std::isfinite(t_val)) {
                    const double xpom_from_def = (electron_Q2_truth + MX2_truth + t_val) / xpom_denominator_truth;
                    if (std::isfinite(xpom_from_def) && xpom_from_def > 0.0) {
                        h_xpom_def_MC->Fill(xpom_from_def);
                        const double xL_for_comp = p.Pz() / beams.p_beam.Pz();
                        const double xpom_from_xL = 1.0 - xL_for_comp;
                        if (std::isfinite(xpom_from_xL) && xpom_from_xL > 0.0) {
                            h_xpom_comp_MC->Fill(xpom_from_xL, xpom_from_def);
                        }
                        if (validTruthX) {
                            const double beta_truth = electron_x_truth / xpom_from_def;
                            if (std::isfinite(beta_truth) && beta_truth > 0.0 && beta_truth <= 1.0) {
                                h_beta_MC->Fill(beta_truth);
                            }
                        }
                    }
                }
            }
        }

        double t_truth_best = -1.0;
        if (has_truth_lead) {
            t_truth_best = TMath::Abs(CalcT(beams.p_beam, lead_truth));
        }
        double xpom_truth_best = -1.0;
        if (xpom_denominator_truth > 0.0 && valid_MX2_truth && std::isfinite(t_truth_best)) {
            xpom_truth_best = (electron_Q2_truth + MX2_truth + t_truth_best) / xpom_denominator_truth;
            if (!(std::isfinite(xpom_truth_best) && xpom_truth_best > 0.0)) {
                xpom_truth_best = -1.0;
            }
        }
        double beta_truth_best = -1.0;
        if (validTruthX && xpom_truth_best > 0.0) {
            beta_truth_best = electron_x_truth / xpom_truth_best;
            if (!(std::isfinite(beta_truth_best) && beta_truth_best > 0.0 && beta_truth_best <= 1.0)) {
                beta_truth_best = -1.0;
            }
        }

        int k_relres = -1;
        if (validTruthQ2 && xpom_truth_best > 0.0 && beta_truth_best > 0.0) {
            k_relres = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_best, beta_truth_best);
        }
        if (k_relres >= 0) {
            occupancy_truth_k[k_relres] += 1;
            h_phase_bin_gen->Fill(k_relres + 1);
            if (is_mc_subsample) {
                h_phase_bin_gen_setA->Fill(k_relres + 1);
            } else {
                h_phase_bin_gen_setB->Fill(k_relres + 1);
            }
            // ResAccum fills (Set A only)
            if (is_mc_subsample) {
            if (validEMQ2) {
                res_Q2_EM_k.Fill(k_relres, (electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth);
            }
            if (validDAQ2) {
                res_Q2_DA_k.Fill(k_relres, (electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth);
            }
            if (validSigmaQ2) {
                res_Q2_Sigma_k.Fill(k_relres, (electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth);
            }
            } // end is_mc_subsample gate for ResAccum
        }
        if (validTruthQ2 && beta_truth_best > 0.0) {
            h_beta_Q2_truth->Fill(beta_truth_best, electron_Q2_truth);
        }
        if (validTruthQ2 && t_truth_best > 0.0) {
            h_t_Q2_truth->Fill(t_truth_best, electron_Q2_truth);
        }
        if (validTruthQ2 && xpom_truth_best > 0.0) {
            h_xpom_Q2_truth->Fill(xpom_truth_best, electron_Q2_truth);
        }
        if (t_truth_best > 0.0 && beta_truth_best > 0.0) {
            h_beta_t_truth->Fill(beta_truth_best, t_truth_best);
        }
        if (t_truth_best > 0.0 && validTruthX) {
            h_xbj_t_truth->Fill(electron_x_truth, t_truth_best);
        }
        if (t_truth_best > 0.0 && xpom_truth_best > 0.0) {
            h_xpom_t_truth->Fill(xpom_truth_best, t_truth_best);
        }
        if (beta_truth_best > 0.0 && xpom_truth_best > 0.0) {
            h_xpom_beta_truth->Fill(xpom_truth_best, beta_truth_best);
        }
        if (beta_truth_best > 0.0 && validTruthX) {
            h_xbj_beta_truth->Fill(electron_x_truth, beta_truth_best);
        }
        if (xpom_truth_best > 0.0 && validTruthX) {
            h_xbj_xpom_truth->Fill(electron_x_truth, xpom_truth_best);
        }

        if (xpom_truth_best > 0.0) {
            h_xpom_truth_all->Fill(xpom_truth_best);
            if (beta_truth_best > 0.0) {
                h_beta_truth_all->Fill(beta_truth_best);
            }
        }
        if (validTruthQ2 && xpom_truth_best > 0.0 && beta_truth_best > 0.0) {
            h_d3sigma_MC->Fill(electron_Q2_truth, beta_truth_best, xpom_truth_best);
        }

        if (has_truth_lead && passTruthDIS) {
            const double lead_theta_mrad = lead_truth.Theta() * 1000.0;
            const bool lead_in_rp_acceptance = PassRPAcceptance(cutCfg, lead_theta_mrad);
            const bool lead_in_b0_acceptance = PassB0Acceptance(cutCfg, lead_theta_mrad);
            if (is_mc_subsample) {
                if (lead_in_b0_acceptance) {
                    if (t_truth_best > 0.0) h_t_truth_mc_B0->Fill(t_truth_best);
                    if (xpom_truth_best > 0.0) h_xpom_truth_mc_B0->Fill(xpom_truth_best);
                    if (beta_truth_best > 0.0) h_beta_truth_mc_B0->Fill(beta_truth_best);
                    h_theta_truth_mc_B0->Fill(lead_theta_mrad);
                    if (validTruthQ2 && xpom_truth_best > 0.0 && beta_truth_best > 0.0) {
                        h_d3sigma_truth_mc_B0->Fill(electron_Q2_truth, beta_truth_best, xpom_truth_best);
                    }
                }

                if (lead_in_rp_acceptance) {
                    if (t_truth_best > 0.0) h_t_truth_mc_RP->Fill(t_truth_best);
                    if (xpom_truth_best > 0.0) h_xpom_truth_mc_RP->Fill(xpom_truth_best);
                    if (beta_truth_best > 0.0) h_beta_truth_mc_RP->Fill(beta_truth_best);
                    h_theta_truth_mc_RP->Fill(lead_theta_mrad);
                    if (validTruthQ2 && xpom_truth_best > 0.0 && beta_truth_best > 0.0) {
                        h_d3sigma_truth_mc_RP->Fill(electron_Q2_truth, beta_truth_best, xpom_truth_best);
                    }
                }
            } else {
                if (lead_in_b0_acceptance) {
                    if (t_truth_best > 0.0) h_t_truth_pdata_B0->Fill(t_truth_best);
                    if (xpom_truth_best > 0.0) h_xpom_truth_pdata_B0->Fill(xpom_truth_best);
                    if (beta_truth_best > 0.0) h_beta_truth_pdata_B0->Fill(beta_truth_best);
                    h_theta_truth_pdata_B0->Fill(lead_theta_mrad);
                }

                if (lead_in_rp_acceptance) {
                    if (t_truth_best > 0.0) h_t_truth_pdata_RP->Fill(t_truth_best);
                    if (xpom_truth_best > 0.0) h_xpom_truth_pdata_RP->Fill(xpom_truth_best);
                    if (beta_truth_best > 0.0) h_beta_truth_pdata_RP->Fill(beta_truth_best);
                    h_theta_truth_pdata_RP->Fill(lead_theta_mrad);
                }
            }
        }

        if (hasTSProtons && tsassoc_rec_id && tsassoc_sim_id && tsre_px_array && tsre_py_array && tsre_pz_array) {
            for (unsigned int j = 0; j < tsassoc_rec_id->GetSize(); j++) {
                auto mc_idx = (*tsassoc_sim_id)[j];
                if (mc_idx >= (unsigned)mc_pdg_array.GetSize()) continue;
                if (mc_genStatus_array[mc_idx] != 1 || mc_pdg_array[mc_idx] != 2212) continue;

                P3MVector p_reco((*tsre_px_array)[j], (*tsre_py_array)[j], (*tsre_pz_array)[j], mc_mass_array[mc_idx]);
                undoAfterburn(p_reco);
                const double theta_reco_mrad = p_reco.Theta() * 1000.0;
                h_theta_all_TS->Fill(theta_reco_mrad);

                if (!PassB0Acceptance(cutCfg, theta_reco_mrad)) continue;
                if (!passRecoDIS) continue;
                h_theta_B0->Fill(theta_reco_mrad);

                P3MVector p_truth(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
                undoAfterburn(p_truth);

                const double t_reco_abs = TMath::Abs(CalcT(beams.p_beam, p_reco));
                const double t_truth_abs = TMath::Abs(CalcT(beams.p_beam, p_truth));
                if (std::isfinite(t_reco_abs)) {
                    h_t_B0->Fill(t_reco_abs);
                    h_dsigma_dt_B0->Fill(t_reco_abs);
                    if (!has_reco_proton) {
                        t_reco_best = t_reco_abs;
                        has_reco_proton = true;
                    }
                }
                if (is_mc_subsample) {
                    if (std::isfinite(t_reco_abs)) h_t_reco_mc_B0->Fill(t_reco_abs);
                    if (std::isfinite(t_truth_abs) && std::isfinite(t_reco_abs)) {
                        fillSameBin1D(h_t_reco_mc_B0, h_t_same_mc_B0, t_truth_abs, t_reco_abs);
                    }
                } else {
                    if (std::isfinite(t_reco_abs)) h_t_reco_pdata_B0->Fill(t_reco_abs);
                }
                // t correlation and resolution (Set A only)
                if (is_mc_subsample) {
                if (std::isfinite(t_truth_abs) && std::isfinite(t_reco_abs)) {
                    h_t_corr_B0->Fill(t_truth_abs, t_reco_abs);
                    g_t_B0->SetPoint(n_g_t_B0++, t_truth_abs, t_reco_abs);
                    if (t_truth_abs > 1e-6) {
                        const double rel = (t_reco_abs - t_truth_abs) / t_truth_abs;
                        h_t_res_B0->Fill(rel);
                        h_t_RelRes_binned_B0->Fill(t_truth_abs, rel);
                    }
                }
                } // end is_mc_subsample gate for B0 t resolution
                if (std::isfinite(t_reco_abs) && valid_MX2_reco) {
                    h_MX2_t_B0->Fill(MX2_reco, t_reco_abs);
                }

                const bool xpom_truth_valid = (xpom_denominator_truth > 0.0 && valid_MX2_truth && std::isfinite(t_truth_abs));
                double xpom_truth_def = -1.0;
                if (xpom_truth_valid) {
                    xpom_truth_def = (electron_Q2_truth + MX2_truth + t_truth_abs) / xpom_denominator_truth;
                    if (!(std::isfinite(xpom_truth_def) && xpom_truth_def > 0.0)) {
                        xpom_truth_def = -1.0;
                    }
                }

                const bool beta_truth_valid = (validTruthX && xpom_truth_def > 0.0);
                const double beta_truth = beta_truth_valid ? (electron_x_truth / xpom_truth_def) : -1.0;
                if (is_mc_subsample && validTruthQ2 && xpom_truth_def > 0.0 && t_truth_abs > 1e-6 && std::isfinite(t_reco_abs)) {
                    const double rel_t = (t_reco_abs - t_truth_abs) / t_truth_abs;
                    h_t_RelRes_vs_xpomQ2_B0->Fill(xpom_truth_def, electron_Q2_truth, rel_t);
                }

                if (!xpom_truth_b0_filled && xpom_truth_def > 0.0) {
                    h_xpom_truth_B0->Fill(xpom_truth_def);
                    if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                        h_beta_truth_B0->Fill(beta_truth);
                    }
                    xpom_truth_b0_filled = true;
                }

                if (!xpom_reco_em_b0_filled && xpom_denominator_EM > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                    const double xpom_reco_def = (electron_Q2_EM + MX2_reco + t_reco_abs) / xpom_denominator_EM;
                    if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                        h_xpom_reco_EM_B0->Fill(xpom_reco_def);
                        h_xpom_def_B0->Fill(xpom_reco_def);
                        const double xL_reco_for_def = p_reco.Pz() / beams.p_beam.Pz();
                        const double xpom_from_xL_for_def = 1.0 - xL_reco_for_def;
                        if (std::isfinite(xpom_from_xL_for_def) && xpom_from_xL_for_def > 0.0) {
                            h_xpom_comp_B0->Fill(xpom_from_xL_for_def, xpom_reco_def);
                        }
                        xpom_reco_em_b0_filled = true;
                        // xpom EM B0 response/resolution (Set A only)
                        if (is_mc_subsample && xpom_truth_def > 0.0) {
                            h_Response_xpom_EM_B0->Fill(xpom_truth_def, xpom_reco_def);
                            g_xpom_EM_B0->SetPoint(n_g_xpom_em_b0++, xpom_truth_def, xpom_reco_def);
                            if (validTruthQ2) {
                                const double rel_xpom = (xpom_reco_def - xpom_truth_def) / xpom_truth_def;
                                h_xpom_RelRes_vs_xpomQ2_B0->Fill(xpom_truth_def, electron_Q2_truth, rel_xpom);
                            }
                        }
                        int k_rel = -1;
                        if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0 &&
                            xpom_truth_def > 0.0 && validTruthQ2) {
                            k_rel = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_def, beta_truth);
                            if (k_rel >= 0) {
                                res_xpom_EM_B0_k.Fill(k_rel, (xpom_reco_def - xpom_truth_def) / xpom_truth_def);
                            }
                        }
                        if (validEMX) {
                            const double beta_reco = electron_x_EM / xpom_reco_def;
                            if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                h_beta_reco_EM_B0->Fill(beta_reco);
                                h_beta_B0->Fill(beta_reco);
                                if (validEMQ2) h_beta_vs_Q2->Fill(electron_Q2_EM, beta_reco);
                                h_beta_vs_xpom->Fill(xpom_reco_def, beta_reco);
                                if (std::isfinite(t_reco_abs) && t_reco_abs > 0.0) h_beta_vs_t->Fill(t_reco_abs, beta_reco);
                                if (validEMQ2) {
                                    h_d3sigma_B0->Fill(electron_Q2_EM, beta_reco, xpom_reco_def);
                                    if (is_mc_subsample) {
                                        h_d3sigma_reco_mc_B0->Fill(electron_Q2_EM, beta_reco, xpom_reco_def);
                                        if (validTruthQ2 && xpom_truth_best > 0.0 && beta_truth_best > 0.0) {
                                            fillSameBin3D(h_d3sigma_reco_mc_B0,
                                                          h_d3sigma_same_mc_B0,
                                                          electron_Q2_truth,
                                                          beta_truth_best,
                                                          xpom_truth_best,
                                                          electron_Q2_EM,
                                                          beta_reco,
                                                          xpom_reco_def);
                                        }
                                    }
                                }
                                // beta EM B0 response/resolution (Set A only)
                                if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    g_beta_EM_B0->SetPoint(n_g_beta_em_b0++, beta_truth, beta_reco);
                                    h_Response_beta_EM_B0->Fill(beta_truth, beta_reco);
                                    h_beta_corr_B0->Fill(beta_truth, beta_reco);
                                    if (beta_truth > 1e-6) {
                                        const double rel_beta = (beta_reco - beta_truth) / beta_truth;
                                        h_beta_res_B0->Fill(rel_beta);
                                        h_beta_RelRes_binned_B0->Fill(beta_truth, rel_beta);
                                        if (validTruthQ2) {
                                            h_beta_RelRes_vs_betaQ2_B0->Fill(beta_truth, electron_Q2_truth, rel_beta);
                                        }
                                    }
                                }
                            }

                            if (is_mc_subsample && k_rel >= 0 && beta_reco > 0.0 && beta_reco <= 1.0) {
                                res_beta_EM_B0_k.Fill(k_rel, (beta_reco - beta_truth) / beta_truth);
                            }
                        }
                        if (is_mc_subsample) {
                            h_xpom_reco_mc_B0->Fill(xpom_reco_def);
                            if (validEMX) {
                                const double beta_reco = electron_x_EM / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_mc_B0->Fill(beta_reco);
                                }
                            }
                        } else {
                            h_xpom_reco_pdata_B0->Fill(xpom_reco_def);
                            if (validEMX) {
                                const double beta_reco = electron_x_EM / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_pdata_B0->Fill(beta_reco);
                                }
                            }
                        }
                        if (!xpom_reco_em_all_filled) {
                            h_xpom_reco_EM_all->Fill(xpom_reco_def);
                            xpom_reco_em_all_filled = true;
                            if (validEMX) {
                                const double beta_reco = electron_x_EM / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_EM_all->Fill(beta_reco);
                                }
                            }
                        }
                    }
                }

                if (!xpom_reco_w2best_b0_filled && xpom_denominator_W2Best > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                    const double xpom_reco_w2best = (electron_Q2_EM + MX2_reco + t_reco_abs) / xpom_denominator_W2Best;
                    if (std::isfinite(xpom_reco_w2best) && xpom_reco_w2best > 0.0) {
                        xpom_reco_w2best_b0_filled = true;
                        // W2Best xpom/beta response/resolution (Set A only)
                        if (is_mc_subsample && xpom_truth_def > 0.0) {
                            h_Response_xpom_W2Best_B0->Fill(xpom_truth_def, xpom_reco_w2best);
                            if (xpom_truth_def > 1e-9) {
                                const double rel_xpom = (xpom_reco_w2best - xpom_truth_def) / xpom_truth_def;
                                h_xpom_RelRes_binned_W2Best_B0->Fill(xpom_truth_def, rel_xpom);
                            }
                        }
                        if (validEMX) {
                            const double beta_reco_w2best = electron_x_EM / xpom_reco_w2best;
                            if (is_mc_subsample && beta_reco_w2best > 0.0 && beta_reco_w2best <= 1.0 &&
                                beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                h_Response_beta_W2Best_B0->Fill(beta_truth, beta_reco_w2best);
                                if (beta_truth > 1e-6) {
                                    const double rel_beta = (beta_reco_w2best - beta_truth) / beta_truth;
                                    h_beta_RelRes_binned_W2Best_B0->Fill(beta_truth, rel_beta);
                                }
                            }
                        }
                    }
                }

                if (!xpom_reco_da_b0_filled && xpom_denominator_DA > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                    const double xpom_reco_def = (electron_Q2_DA + MX2_reco + t_reco_abs) / xpom_denominator_DA;
                    if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                        h_xpom_reco_DA_B0->Fill(xpom_reco_def);
                        xpom_reco_da_b0_filled = true;
                        // DA B0 response/resolution (Set A only)
                        if (is_mc_subsample && xpom_truth_def > 0.0) {
                            h_Response_xpom_DA_B0->Fill(xpom_truth_def, xpom_reco_def);
                            g_xpom_DA_B0->SetPoint(n_g_xpom_da_b0++, xpom_truth_def, xpom_reco_def);
                        }
                        int k_rel = -1;
                        if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0 &&
                            xpom_truth_def > 0.0 && validTruthQ2) {
                            k_rel = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_def, beta_truth);
                            if (k_rel >= 0) {
                                res_xpom_DA_B0_k.Fill(k_rel, (xpom_reco_def - xpom_truth_def) / xpom_truth_def);
                            }
                        }
                        if (validDAX) {
                            const double beta_reco = electron_x_DA / xpom_reco_def;
                            if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                h_beta_reco_DA_B0->Fill(beta_reco);
                                if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    g_beta_DA_B0->SetPoint(n_g_beta_da_b0++, beta_truth, beta_reco);
                                    h_Response_beta_DA_B0->Fill(beta_truth, beta_reco);
                                }
                            }

                            if (is_mc_subsample && k_rel >= 0 && beta_reco > 0.0 && beta_reco <= 1.0) {
                                res_beta_DA_B0_k.Fill(k_rel, (beta_reco - beta_truth) / beta_truth);
                            }
                        }
                        if (!xpom_reco_da_all_filled) {
                            h_xpom_reco_DA_all->Fill(xpom_reco_def);
                            xpom_reco_da_all_filled = true;
                            if (validDAX) {
                                const double beta_reco = electron_x_DA / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_DA_all->Fill(beta_reco);
                                }
                            }
                        }
                    }
                }

                if (!xpom_reco_sigma_b0_filled && xpom_denominator_Sigma > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                    const double xpom_reco_def = (electron_Q2_Sigma + MX2_reco + t_reco_abs) / xpom_denominator_Sigma;
                    if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                        h_xpom_reco_Sigma_B0->Fill(xpom_reco_def);
                        xpom_reco_sigma_b0_filled = true;
                        // Sigma B0 response/resolution (Set A only)
                        if (is_mc_subsample && xpom_truth_def > 0.0) {
                            h_Response_xpom_Sigma_B0->Fill(xpom_truth_def, xpom_reco_def);
                            g_xpom_Sigma_B0->SetPoint(n_g_xpom_sigma_b0++, xpom_truth_def, xpom_reco_def);
                        }
                        int k_rel = -1;
                        if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0 &&
                            xpom_truth_def > 0.0 && validTruthQ2) {
                            k_rel = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_def, beta_truth);
                            if (k_rel >= 0) {
                                res_xpom_Sigma_B0_k.Fill(k_rel, (xpom_reco_def - xpom_truth_def) / xpom_truth_def);
                            }
                        }
                        if (hasSigmaMethod && validSigmaX) {
                            const double beta_reco = electron_x_Sigma / xpom_reco_def;
                            if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                h_beta_reco_Sigma_B0->Fill(beta_reco);
                                if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    g_beta_Sigma_B0->SetPoint(n_g_beta_sigma_b0++, beta_truth, beta_reco);
                                    h_Response_beta_Sigma_B0->Fill(beta_truth, beta_reco);
                                }
                            }

                            if (is_mc_subsample && k_rel >= 0 && beta_reco > 0.0 && beta_reco <= 1.0) {
                                res_beta_Sigma_B0_k.Fill(k_rel, (beta_reco - beta_truth) / beta_truth);
                            }
                        }
                        if (!xpom_reco_sigma_all_filled) {
                            h_xpom_reco_Sigma_all->Fill(xpom_reco_def);
                            xpom_reco_sigma_all_filled = true;
                            if (hasSigmaMethod && validSigmaX) {
                                const double beta_reco = electron_x_Sigma / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_Sigma_all->Fill(beta_reco);
                                }
                            }
                        }
                    }
                }

                const double xL_reco_pz = p_reco.Pz() / beams.p_beam.Pz();
                const double xL_truth_pz = p_truth.Pz() / beams.p_beam.Pz();
                if (std::isfinite(xL_reco_pz)) {
                    h_xL_B0->Fill(xL_reco_pz);
                }
                if (std::isfinite(xL_truth_pz) && std::isfinite(xL_reco_pz)) {
                    const double xpom_reco_from_xL = 1.0 - xL_reco_pz;
                    if (xpom_reco_from_xL > 0.0 && std::isfinite(xpom_reco_from_xL)) {
                        h_xpom_B0->Fill(xpom_reco_from_xL);
                    }
                    // xL/xpom B0 correlation and resolution (Set A only)
                    if (is_mc_subsample) {
                    g_xL_B0->SetPoint(n_g_xL_B0++, xL_truth_pz, xL_reco_pz);
                    h_xL_corr_B0->Fill(xL_truth_pz, xL_reco_pz);
                    if (xL_truth_pz > 1e-6) {
                        const double rel_xL = (xL_reco_pz - xL_truth_pz) / xL_truth_pz;
                        h_xL_res_B0->Fill(rel_xL);
                        h_xL_RelRes_binned_B0->Fill(xL_truth_pz, rel_xL);
                        if (validTruthQ2) {
                            h_xL_RelRes_vs_xLQ2_B0->Fill(xL_truth_pz, electron_Q2_truth, rel_xL);
                        }
                    }
                    const double xpom_truth_from_xL = 1.0 - xL_truth_pz;
                    if (xpom_truth_from_xL > 1e-9 && std::isfinite(xpom_truth_from_xL) &&
                        xpom_reco_from_xL > 0.0 && std::isfinite(xpom_reco_from_xL)) {
                        h_xpom_corr_B0->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
                        const double rel_xpom = (xpom_reco_from_xL - xpom_truth_from_xL) / xpom_truth_from_xL;
                        h_xpom_res_B0->Fill(rel_xpom);
                        h_xpom_RelRes_binned_B0->Fill(xpom_truth_from_xL, rel_xpom);
                    }
                    } // end is_mc_subsample gate for B0 xL/xpom resolution
                }
                if (is_mc_subsample) {
                    h_theta_reco_mc_B0->Fill(theta_reco_mrad);
                } else {
                    h_theta_reco_pdata_B0->Fill(theta_reco_mrad);
                }
            }
        }

        if (hasRPProtons && rp_px_array && rp_py_array && rp_pz_array && rp_mass_array && rp_pdg_array) {
            for (int j = 0; j < rp_px_array->GetSize(); j++) {
                if ((*rp_pdg_array)[j] != 2212) continue;

                P3MVector p_rp((*rp_px_array)[j], (*rp_py_array)[j], (*rp_pz_array)[j], (*rp_mass_array)[j]);
                const double theta_mrad = p_rp.Theta() * 1000.0;
                if (!PassRPAcceptance(cutCfg, theta_mrad)) continue;
                if (!passRecoDIS) continue;
                h_theta_RP->Fill(theta_mrad);

                P3MVector* best_match = nullptr;
                double best_dr = 0.1;
                for (auto& p_truth : truth_protons) {
                    double dr = TMath::Sqrt(TMath::Power(p_rp.Theta() - p_truth.Theta(), 2) +
                                            TMath::Power(p_rp.Phi() - p_truth.Phi(), 2));
                    if (dr < best_dr) {
                        best_dr = dr;
                        best_match = &p_truth;
                    }
                }

                if (best_match) {
                    const double t_reco_abs = TMath::Abs(CalcT(beams.p_beam, p_rp));
                    const double t_truth_abs = TMath::Abs(CalcT(beams.p_beam, *best_match));
                    if (std::isfinite(t_reco_abs)) {
                        h_t_RP_histo->Fill(t_reco_abs);
                        h_dsigma_dt_RP->Fill(t_reco_abs);
                        if (!has_reco_proton) {
                            t_reco_best = t_reco_abs;
                            has_reco_proton = true;
                        }
                    }
                    if (is_mc_subsample) {
                        if (std::isfinite(t_reco_abs)) h_t_reco_mc_RP->Fill(t_reco_abs);
                        if (std::isfinite(t_truth_abs) && std::isfinite(t_reco_abs)) {
                            fillSameBin1D(h_t_reco_mc_RP, h_t_same_mc_RP, t_truth_abs, t_reco_abs);
                        }
                    } else {
                        if (std::isfinite(t_reco_abs)) h_t_reco_pdata_RP->Fill(t_reco_abs);
                    }
                    // t correlation and resolution (Set A only)
                    if (is_mc_subsample) {
                    if (std::isfinite(t_truth_abs) && std::isfinite(t_reco_abs)) {
                        h_t_corr_RP->Fill(t_truth_abs, t_reco_abs);
                        g_t_RP->SetPoint(n_g_t_RP++, t_truth_abs, t_reco_abs);
                        if (t_truth_abs > 1e-6) {
                            const double rel = (t_reco_abs - t_truth_abs) / t_truth_abs;
                            h_t_res_RP->Fill(rel);
                            h_t_RelRes_binned_RP->Fill(t_truth_abs, rel);
                        }
                    }
                    } // end is_mc_subsample gate for RP t resolution
                    if (std::isfinite(t_reco_abs) && valid_MX2_reco) {
                        h_MX2_t_RP->Fill(MX2_reco, t_reco_abs);
                    }

                    const bool xpom_truth_valid = (xpom_denominator_truth > 0.0 && valid_MX2_truth && std::isfinite(t_truth_abs));
                    double xpom_truth_def = -1.0;
                    if (xpom_truth_valid) {
                        xpom_truth_def = (electron_Q2_truth + MX2_truth + t_truth_abs) / xpom_denominator_truth;
                        if (!(std::isfinite(xpom_truth_def) && xpom_truth_def > 0.0)) {
                            xpom_truth_def = -1.0;
                        }
                    }

                    const bool beta_truth_valid = (validTruthX && xpom_truth_def > 0.0);
                    const double beta_truth = beta_truth_valid ? (electron_x_truth / xpom_truth_def) : -1.0;
                    if (is_mc_subsample && validTruthQ2 && xpom_truth_def > 0.0 && t_truth_abs > 1e-6 && std::isfinite(t_reco_abs)) {
                        const double rel_t = (t_reco_abs - t_truth_abs) / t_truth_abs;
                        h_t_RelRes_vs_xpomQ2_RP->Fill(xpom_truth_def, electron_Q2_truth, rel_t);
                    }

                    if (!xpom_truth_rp_filled && xpom_truth_def > 0.0) {
                        h_xpom_truth_RP->Fill(xpom_truth_def);
                        if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                            h_beta_truth_RP->Fill(beta_truth);
                        }
                        xpom_truth_rp_filled = true;
                    }

                    if (!xpom_reco_em_rp_filled && xpom_denominator_EM > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                        const double xpom_reco_def = (electron_Q2_EM + MX2_reco + t_reco_abs) / xpom_denominator_EM;
                        if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                            h_xpom_reco_EM_RP->Fill(xpom_reco_def);
                            h_xpom_def_RP->Fill(xpom_reco_def);
                            const double xL_reco_for_def = p_rp.Pz() / beams.p_beam.Pz();
                            const double xpom_from_xL_for_def = 1.0 - xL_reco_for_def;
                            if (std::isfinite(xpom_from_xL_for_def) && xpom_from_xL_for_def > 0.0) {
                                h_xpom_comp_RP->Fill(xpom_from_xL_for_def, xpom_reco_def);
                            }
                            xpom_reco_em_rp_filled = true;
                            // xpom EM RP response/resolution (Set A only)
                            if (is_mc_subsample && xpom_truth_def > 0.0) {
                                h_Response_xpom_EM_RP->Fill(xpom_truth_def, xpom_reco_def);
                                g_xpom_EM_RP->SetPoint(n_g_xpom_em_rp++, xpom_truth_def, xpom_reco_def);
                                if (validTruthQ2) {
                                    const double rel_xpom = (xpom_reco_def - xpom_truth_def) / xpom_truth_def;
                                    h_xpom_RelRes_vs_xpomQ2_RP->Fill(xpom_truth_def, electron_Q2_truth, rel_xpom);
                                }
                            }
                            int k_rel = -1;
                            if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0 &&
                                xpom_truth_def > 0.0 && validTruthQ2) {
                                k_rel = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_def, beta_truth);
                                if (k_rel >= 0) {
                                    res_xpom_EM_RP_k.Fill(k_rel, (xpom_reco_def - xpom_truth_def) / xpom_truth_def);
                                }
                            }
                            if (validEMX) {
                                const double beta_reco = electron_x_EM / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_EM_RP->Fill(beta_reco);
                                    h_beta_RP->Fill(beta_reco);
                                    if (validEMQ2) h_beta_vs_Q2->Fill(electron_Q2_EM, beta_reco);
                                    h_beta_vs_xpom->Fill(xpom_reco_def, beta_reco);
                                    if (std::isfinite(t_reco_abs) && t_reco_abs > 0.0) h_beta_vs_t->Fill(t_reco_abs, beta_reco);
                                    if (validEMQ2) {
                                        h_d3sigma_RP->Fill(electron_Q2_EM, beta_reco, xpom_reco_def);
                                        if (is_mc_subsample) {
                                            h_d3sigma_reco_mc_RP->Fill(electron_Q2_EM, beta_reco, xpom_reco_def);
                                            if (validTruthQ2 && xpom_truth_best > 0.0 && beta_truth_best > 0.0) {
                                                fillSameBin3D(h_d3sigma_reco_mc_RP,
                                                              h_d3sigma_same_mc_RP,
                                                              electron_Q2_truth,
                                                              beta_truth_best,
                                                              xpom_truth_best,
                                                              electron_Q2_EM,
                                                              beta_reco,
                                                              xpom_reco_def);
                                            }
                                        }
                                    }
                                    // beta EM RP response/resolution (Set A only)
                                    if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                        g_beta_EM_RP->SetPoint(n_g_beta_em_rp++, beta_truth, beta_reco);
                                        h_Response_beta_EM_RP->Fill(beta_truth, beta_reco);
                                        h_beta_corr_RP->Fill(beta_truth, beta_reco);
                                        if (beta_truth > 1e-6) {
                                            const double rel_beta = (beta_reco - beta_truth) / beta_truth;
                                            h_beta_res_RP->Fill(rel_beta);
                                            h_beta_RelRes_binned_RP->Fill(beta_truth, rel_beta);
                                            if (validTruthQ2) {
                                                h_beta_RelRes_vs_betaQ2_RP->Fill(beta_truth, electron_Q2_truth, rel_beta);
                                            }
                                        }
                                    }
                                }
                                if (is_mc_subsample && k_rel >= 0 && beta_reco > 0.0 && beta_reco <= 1.0) {
                                    res_beta_EM_RP_k.Fill(k_rel, (beta_reco - beta_truth) / beta_truth);
                                }
                            }
                            if (is_mc_subsample) {
                                h_xpom_reco_mc_RP->Fill(xpom_reco_def);
                                if (validEMX) {
                                    const double beta_reco = electron_x_EM / xpom_reco_def;
                                    if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                        h_beta_reco_mc_RP->Fill(beta_reco);
                                    }
                                }
                            } else {
                                h_xpom_reco_pdata_RP->Fill(xpom_reco_def);
                                if (validEMX) {
                                    const double beta_reco = electron_x_EM / xpom_reco_def;
                                    if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                        h_beta_reco_pdata_RP->Fill(beta_reco);
                                    }
                                }
                            }
                            if (!xpom_reco_em_all_filled) {
                                h_xpom_reco_EM_all->Fill(xpom_reco_def);
                                xpom_reco_em_all_filled = true;
                                if (validEMX) {
                                    const double beta_reco = electron_x_EM / xpom_reco_def;
                                    if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                        h_beta_reco_EM_all->Fill(beta_reco);
                                    }
                                }
                            }
                        }
                    }

                    if (!xpom_reco_w2best_rp_filled && xpom_denominator_W2Best > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                        const double xpom_reco_w2best = (electron_Q2_EM + MX2_reco + t_reco_abs) / xpom_denominator_W2Best;
                        if (std::isfinite(xpom_reco_w2best) && xpom_reco_w2best > 0.0) {
                            xpom_reco_w2best_rp_filled = true;
                            // W2Best RP response/resolution (Set A only)
                            if (is_mc_subsample && xpom_truth_def > 0.0) {
                                h_Response_xpom_W2Best_RP->Fill(xpom_truth_def, xpom_reco_w2best);
                                if (xpom_truth_def > 1e-9) {
                                    const double rel_xpom = (xpom_reco_w2best - xpom_truth_def) / xpom_truth_def;
                                    h_xpom_RelRes_binned_W2Best_RP->Fill(xpom_truth_def, rel_xpom);
                                }
                            }
                            if (validEMX) {
                                const double beta_reco_w2best = electron_x_EM / xpom_reco_w2best;
                                if (is_mc_subsample && beta_reco_w2best > 0.0 && beta_reco_w2best <= 1.0 &&
                                    beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    h_Response_beta_W2Best_RP->Fill(beta_truth, beta_reco_w2best);
                                    if (beta_truth > 1e-6) {
                                        const double rel_beta = (beta_reco_w2best - beta_truth) / beta_truth;
                                        h_beta_RelRes_binned_W2Best_RP->Fill(beta_truth, rel_beta);
                                    }
                                }
                            }
                        }
                    }

                    if (!xpom_reco_da_rp_filled && xpom_denominator_DA > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                        const double xpom_reco_def = (electron_Q2_DA + MX2_reco + t_reco_abs) / xpom_denominator_DA;
                        if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                            h_xpom_reco_DA_RP->Fill(xpom_reco_def);
                            xpom_reco_da_rp_filled = true;
                            // DA RP response/resolution (Set A only)
                            if (is_mc_subsample && xpom_truth_def > 0.0) {
                                h_Response_xpom_DA_RP->Fill(xpom_truth_def, xpom_reco_def);
                                g_xpom_DA_RP->SetPoint(n_g_xpom_da_rp++, xpom_truth_def, xpom_reco_def);
                            }
                            int k_rel = -1;
                            if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0 &&
                                xpom_truth_def > 0.0 && validTruthQ2) {
                                k_rel = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_def, beta_truth);
                                if (k_rel >= 0) {
                                    res_xpom_DA_RP_k.Fill(k_rel, (xpom_reco_def - xpom_truth_def) / xpom_truth_def);
                                }
                            }
                            if (validDAX) {
                                const double beta_reco = electron_x_DA / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_DA_RP->Fill(beta_reco);
                                    if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                        g_beta_DA_RP->SetPoint(n_g_beta_da_rp++, beta_truth, beta_reco);
                                        h_Response_beta_DA_RP->Fill(beta_truth, beta_reco);
                                    }
                                }
                                if (is_mc_subsample && k_rel >= 0 && beta_reco > 0.0 && beta_reco <= 1.0) {
                                    res_beta_DA_RP_k.Fill(k_rel, (beta_reco - beta_truth) / beta_truth);
                                }
                            }
                            if (!xpom_reco_da_all_filled) {
                                h_xpom_reco_DA_all->Fill(xpom_reco_def);
                                xpom_reco_da_all_filled = true;
                                if (validDAX) {
                                    const double beta_reco = electron_x_DA / xpom_reco_def;
                                    if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                        h_beta_reco_DA_all->Fill(beta_reco);
                                    }
                                }
                            }
                        }
                    }

                    if (!xpom_reco_sigma_rp_filled && xpom_denominator_Sigma > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                        const double xpom_reco_def = (electron_Q2_Sigma + MX2_reco + t_reco_abs) / xpom_denominator_Sigma;
                        if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                            h_xpom_reco_Sigma_RP->Fill(xpom_reco_def);
                            xpom_reco_sigma_rp_filled = true;
                            // Sigma RP response/resolution (Set A only)
                            if (is_mc_subsample && xpom_truth_def > 0.0) {
                                h_Response_xpom_Sigma_RP->Fill(xpom_truth_def, xpom_reco_def);
                                g_xpom_Sigma_RP->SetPoint(n_g_xpom_sigma_rp++, xpom_truth_def, xpom_reco_def);
                            }
                            int k_rel = -1;
                            if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0 &&
                                xpom_truth_def > 0.0 && validTruthQ2) {
                                k_rel = GetRelResGlobalBin(electron_Q2_truth, xpom_truth_def, beta_truth);
                                if (k_rel >= 0) {
                                    res_xpom_Sigma_RP_k.Fill(k_rel, (xpom_reco_def - xpom_truth_def) / xpom_truth_def);
                                }
                            }
                            if (hasSigmaMethod && validSigmaX) {
                                const double beta_reco = electron_x_Sigma / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_Sigma_RP->Fill(beta_reco);
                                    if (is_mc_subsample && beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                        g_beta_Sigma_RP->SetPoint(n_g_beta_sigma_rp++, beta_truth, beta_reco);
                                        h_Response_beta_Sigma_RP->Fill(beta_truth, beta_reco);
                                    }
                                }
                                if (is_mc_subsample && k_rel >= 0 && beta_reco > 0.0 && beta_reco <= 1.0) {
                                    res_beta_Sigma_RP_k.Fill(k_rel, (beta_reco - beta_truth) / beta_truth);
                                }
                            }
                            if (!xpom_reco_sigma_all_filled) {
                                h_xpom_reco_Sigma_all->Fill(xpom_reco_def);
                                xpom_reco_sigma_all_filled = true;
                                if (hasSigmaMethod && validSigmaX) {
                                    const double beta_reco = electron_x_Sigma / xpom_reco_def;
                                    if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                        h_beta_reco_Sigma_all->Fill(beta_reco);
                                    }
                                }
                            }
                        }
                    }

                    const double xL_reco_pz = p_rp.Pz() / beams.p_beam.Pz();
                    const double xL_truth_pz = best_match->Pz() / beams.p_beam.Pz();
                    if (std::isfinite(xL_reco_pz)) {
                        h_xL_RP->Fill(xL_reco_pz);
                    }
                    if (std::isfinite(xL_truth_pz) && std::isfinite(xL_reco_pz)) {
                        const double xpom_reco_from_xL = 1.0 - xL_reco_pz;
                        if (xpom_reco_from_xL > 0.0 && std::isfinite(xpom_reco_from_xL)) {
                            h_xpom_RP->Fill(xpom_reco_from_xL);
                        }
                        // xL/xpom RP correlation and resolution (Set A only)
                        if (is_mc_subsample) {
                        g_xL_RP->SetPoint(n_g_xL_RP++, xL_truth_pz, xL_reco_pz);
                        h_xL_corr_RP->Fill(xL_truth_pz, xL_reco_pz);
                        if (xL_truth_pz > 1e-6) {
                            const double rel_xL = (xL_reco_pz - xL_truth_pz) / xL_truth_pz;
                            h_xL_res_RP->Fill(rel_xL);
                            h_xL_RelRes_binned_RP->Fill(xL_truth_pz, rel_xL);
                            if (validTruthQ2) {
                                h_xL_RelRes_vs_xLQ2_RP->Fill(xL_truth_pz, electron_Q2_truth, rel_xL);
                            }
                        }
                        const double xpom_truth_from_xL = 1.0 - xL_truth_pz;
                        if (xpom_truth_from_xL > 1e-9 && std::isfinite(xpom_truth_from_xL) &&
                            xpom_reco_from_xL > 0.0 && std::isfinite(xpom_reco_from_xL)) {
                            h_xpom_corr_RP->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
                            const double rel_xpom = (xpom_reco_from_xL - xpom_truth_from_xL) / xpom_truth_from_xL;
                            h_xpom_res_RP->Fill(rel_xpom);
                            h_xpom_RelRes_binned_RP->Fill(xpom_truth_from_xL, rel_xpom);
                        }
                        } // end is_mc_subsample gate for RP xL/xpom resolution
                    }
                    if (is_mc_subsample) {
                        h_theta_reco_mc_RP->Fill(theta_mrad);
                    } else {
                        h_theta_reco_pdata_RP->Fill(theta_mrad);
                    }
                }
            }
        }

        double xpom_reco_best = -1.0;
        if (has_reco_proton && xpom_denominator_EM > 0.0 && valid_MX2_reco && std::isfinite(t_reco_best)) {
            xpom_reco_best = (electron_Q2_EM + MX2_reco + t_reco_best) / xpom_denominator_EM;
            if (!(std::isfinite(xpom_reco_best) && xpom_reco_best > 0.0)) {
                xpom_reco_best = -1.0;
            }
        }
        double beta_reco_best = -1.0;
        if (validEMX && xpom_reco_best > 0.0) {
            beta_reco_best = electron_x_EM / xpom_reco_best;
            if (!(std::isfinite(beta_reco_best) && beta_reco_best > 0.0 && beta_reco_best <= 1.0)) {
                beta_reco_best = -1.0;
            }
        }

        if (validEMQ2 && beta_reco_best > 0.0) {
            h_beta_Q2_reco->Fill(beta_reco_best, electron_Q2_EM);
        }
        if (validEMQ2 && t_reco_best > 0.0) {
            h_t_Q2_reco->Fill(t_reco_best, electron_Q2_EM);
        }
        if (validEMQ2 && xpom_reco_best > 0.0) {
            h_xpom_Q2_reco->Fill(xpom_reco_best, electron_Q2_EM);
        }
        if (t_reco_best > 0.0 && beta_reco_best > 0.0) {
            h_beta_t_reco->Fill(beta_reco_best, t_reco_best);
        }
        if (t_reco_best > 0.0 && validEMX) {
            h_xbj_t_reco->Fill(electron_x_EM, t_reco_best);
        }
        if (t_reco_best > 0.0 && xpom_reco_best > 0.0) {
            h_xpom_t_reco->Fill(xpom_reco_best, t_reco_best);
        }
        if (beta_reco_best > 0.0 && xpom_reco_best > 0.0) {
            h_xpom_beta_reco->Fill(xpom_reco_best, beta_reco_best);
        }
        if (beta_reco_best > 0.0 && validEMX) {
            h_xbj_beta_reco->Fill(electron_x_EM, beta_reco_best);
        }
        if (xpom_reco_best > 0.0 && validEMX) {
            h_xbj_xpom_reco->Fill(electron_x_EM, xpom_reco_best);
        }
        if (validEMQ2 && xpom_reco_best > 0.0 && beta_reco_best > 0.0) {
            h_phase3D_reco->Fill(electron_Q2_EM, xpom_reco_best, beta_reco_best);
            if (h_phase3D_reco_yaml) {
                h_phase3D_reco_yaml->Fill(electron_Q2_EM, xpom_reco_best, beta_reco_best);
            }

            const int k_meas_relres = GetRelResGlobalBin(electron_Q2_EM, xpom_reco_best, beta_reco_best);
            if (k_meas_relres >= 0) {
                h_phase_bin_meas->Fill(k_meas_relres + 1);
                if (is_mc_subsample) {
                    h_phase_bin_meas_setA->Fill(k_meas_relres + 1);
                } else {
                    h_phase_bin_meas_setB->Fill(k_meas_relres + 1);
                }
                if (k_relres >= 0 && k_relres == k_meas_relres) {
                    h_phase_bin_gen_meas_same->Fill(k_relres + 1);
                    if (is_mc_subsample) {
                        h_phase_bin_gen_meas_same_setA->Fill(k_relres + 1);
                    } else {
                        h_phase_bin_gen_meas_same_setB->Fill(k_relres + 1);
                    }
                }
            }
        }

        if (is_mc_subsample) {
            if (validTruthQ2) h_Q2_truth_mc->Fill(electron_Q2_truth);
            if (validEMQ2) h_Q2_reco_mc->Fill(electron_Q2_EM);
            if (validTruthX) h_x_truth_mc->Fill(electron_x_truth);
            if (validEMX) h_x_reco_mc->Fill(electron_x_EM);
            if (validTruthY) h_y_truth_mc->Fill(electron_y_truth);
            if (validEMY) h_y_reco_mc->Fill(electron_y_EM);
            if (valid_MX2_truth) h_MX2_truth_mc->Fill(MX2_truth);
            if (valid_MX2_reco) h_MX2_reco_mc->Fill(MX2_reco);
            if (t_truth_best > 0.0) h_t_truth_mc->Fill(t_truth_best);
            if (t_reco_best > 0.0) h_t_reco_mc->Fill(t_reco_best);
            if (xpom_truth_best > 0.0) h_xpom_truth_mc->Fill(xpom_truth_best);
            if (xpom_reco_best > 0.0) h_xpom_reco_mc->Fill(xpom_reco_best);
            if (beta_truth_best > 0.0) h_beta_truth_mc->Fill(beta_truth_best);
            if (beta_reco_best > 0.0) h_beta_reco_mc->Fill(beta_reco_best);
            if (validWTruth) h_W2_truth_mc->Fill(W2_truth);
            if (validWBest) h_W2_reco_mc->Fill(W2_best);
        } else {
            if (validEMQ2) h_Q2_reco_pdata->Fill(electron_Q2_EM);
            if (validEMX) h_x_reco_pdata->Fill(electron_x_EM);
            if (validEMY) h_y_reco_pdata->Fill(electron_y_EM);
            if (valid_MX2_reco) h_MX2_reco_pdata->Fill(MX2_reco);
            if (t_reco_best > 0.0) h_t_reco_pdata->Fill(t_reco_best);
            if (xpom_reco_best > 0.0) h_xpom_reco_pdata->Fill(xpom_reco_best);
            if (beta_reco_best > 0.0) h_beta_reco_pdata->Fill(beta_reco_best);
            if (validWBest) h_W2_reco_pdata->Fill(W2_best);
            if (validTruthQ2) h_Q2_truth_pdata->Fill(electron_Q2_truth);
            if (validTruthX) h_x_truth_pdata->Fill(electron_x_truth);
            if (validTruthY) h_y_truth_pdata->Fill(electron_y_truth);
            if (valid_MX2_truth) h_MX2_truth_pdata->Fill(MX2_truth);
            if (validWTruth) h_W2_truth_pdata->Fill(W2_truth);
            if (t_truth_best > 0.0) h_t_truth_pdata_all->Fill(t_truth_best);
            if (xpom_truth_best > 0.0) h_xpom_truth_pdata_all->Fill(xpom_truth_best);
            if (beta_truth_best > 0.0) h_beta_truth_pdata_all->Fill(beta_truth_best);
        }

    }

    // Calculate and scale differential cross sections d(sigma)/dt.
    const double sigma_total = 4.22; // nb for 10x100 GeV configuration
    const double N_gen = static_cast<double>(std::max<Long64_t>(1, nProcessed));
    const double scale_factor = sigma_total / N_gen;

    // Apply Set-A purity/efficiency corrections in each bin:
    // N_corr = N_meas * (purity/efficiency), with
    // purity = N_same / N_reco and efficiency = N_same / N_truth.
    for (int i = 1; i <= h_dsigma_dt_B0->GetNbinsX(); ++i) {
        const double nTruth = h_t_truth_mc_B0->GetBinContent(i);
        const double nReco = h_t_reco_mc_B0->GetBinContent(i);
        const double nSame = h_t_same_mc_B0->GetBinContent(i);
        double corr = 0.0;
        if (nTruth > 0.0 && nReco > 0.0) {
            if (nSame > 0.0) {
                const double eff = nSame / nTruth;
                const double pur = nSame / nReco;
                corr = (eff > 0.0) ? (pur / eff) : 0.0;
            } else {
                // Same-bin occupancy can be empty in sparse bins; fall back to truth/reco.
                corr = nTruth / nReco;
            }
        }
        if (corr > 0.0) {
            h_dsigma_dt_B0->SetBinContent(i, h_dsigma_dt_B0->GetBinContent(i) * corr);
            h_dsigma_dt_B0->SetBinError(i, h_dsigma_dt_B0->GetBinError(i) * corr);
        } else {
            h_dsigma_dt_B0->SetBinContent(i, 0.0);
            h_dsigma_dt_B0->SetBinError(i, 0.0);
        }
    }
    for (int i = 1; i <= h_dsigma_dt_RP->GetNbinsX(); ++i) {
        const double nTruth = h_t_truth_mc_RP->GetBinContent(i);
        const double nReco = h_t_reco_mc_RP->GetBinContent(i);
        const double nSame = h_t_same_mc_RP->GetBinContent(i);
        double corr = 0.0;
        if (nTruth > 0.0 && nReco > 0.0) {
            if (nSame > 0.0) {
                const double eff = nSame / nTruth;
                const double pur = nSame / nReco;
                corr = (eff > 0.0) ? (pur / eff) : 0.0;
            } else {
                corr = nTruth / nReco;
            }
        }
        if (corr > 0.0) {
            h_dsigma_dt_RP->SetBinContent(i, h_dsigma_dt_RP->GetBinContent(i) * corr);
            h_dsigma_dt_RP->SetBinError(i, h_dsigma_dt_RP->GetBinError(i) * corr);
        } else {
            h_dsigma_dt_RP->SetBinContent(i, 0.0);
            h_dsigma_dt_RP->SetBinError(i, 0.0);
        }
    }

    for (int ix = 1; ix <= h_d3sigma_B0->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_d3sigma_B0->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h_d3sigma_B0->GetNbinsZ(); ++iz) {
                const double nTruth = h_d3sigma_truth_mc_B0->GetBinContent(ix, iy, iz);
                const double nReco = h_d3sigma_reco_mc_B0->GetBinContent(ix, iy, iz);
                const double nSame = h_d3sigma_same_mc_B0->GetBinContent(ix, iy, iz);
                double corr = 0.0;
                if (nTruth > 0.0 && nReco > 0.0) {
                    if (nSame > 0.0) {
                        const double eff = nSame / nTruth;
                        const double pur = nSame / nReco;
                        corr = (eff > 0.0) ? (pur / eff) : 0.0;
                    } else {
                        corr = nTruth / nReco;
                    }
                }
                if (corr > 0.0) {
                    h_d3sigma_B0->SetBinContent(ix, iy, iz, h_d3sigma_B0->GetBinContent(ix, iy, iz) * corr);
                    h_d3sigma_B0->SetBinError(ix, iy, iz, h_d3sigma_B0->GetBinError(ix, iy, iz) * corr);
                } else {
                    h_d3sigma_B0->SetBinContent(ix, iy, iz, 0.0);
                    h_d3sigma_B0->SetBinError(ix, iy, iz, 0.0);
                }
            }
        }
    }
    for (int ix = 1; ix <= h_d3sigma_RP->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_d3sigma_RP->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h_d3sigma_RP->GetNbinsZ(); ++iz) {
                const double nTruth = h_d3sigma_truth_mc_RP->GetBinContent(ix, iy, iz);
                const double nReco = h_d3sigma_reco_mc_RP->GetBinContent(ix, iy, iz);
                const double nSame = h_d3sigma_same_mc_RP->GetBinContent(ix, iy, iz);
                double corr = 0.0;
                if (nTruth > 0.0 && nReco > 0.0) {
                    if (nSame > 0.0) {
                        const double eff = nSame / nTruth;
                        const double pur = nSame / nReco;
                        corr = (eff > 0.0) ? (pur / eff) : 0.0;
                    } else {
                        corr = nTruth / nReco;
                    }
                }
                if (corr > 0.0) {
                    h_d3sigma_RP->SetBinContent(ix, iy, iz, h_d3sigma_RP->GetBinContent(ix, iy, iz) * corr);
                    h_d3sigma_RP->SetBinError(ix, iy, iz, h_d3sigma_RP->GetBinError(ix, iy, iz) * corr);
                } else {
                    h_d3sigma_RP->SetBinContent(ix, iy, iz, 0.0);
                    h_d3sigma_RP->SetBinError(ix, iy, iz, 0.0);
                }
            }
        }
    }

    for(int i = 1; i <= h_dsigma_dt_MC->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_MC->GetBinContent(i);
        double bin_error = h_dsigma_dt_MC->GetBinError(i);
        double bin_width = h_dsigma_dt_MC->GetBinWidth(i);

        double dsigma_dt = (bin_content * scale_factor) / bin_width;
        double dsigma_dt_error = (bin_error * scale_factor) / bin_width;

        h_dsigma_dt_MC->SetBinContent(i, dsigma_dt);
        h_dsigma_dt_MC->SetBinError(i, dsigma_dt_error);
    }

    for(int i = 1; i <= h_dsigma_dt_B0->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_B0->GetBinContent(i);
        double bin_error = h_dsigma_dt_B0->GetBinError(i);
        double bin_width = h_dsigma_dt_B0->GetBinWidth(i);

        double dsigma_dt = (bin_content * scale_factor) / bin_width;
        double dsigma_dt_error = (bin_error * scale_factor) / bin_width;

        h_dsigma_dt_B0->SetBinContent(i, dsigma_dt);
        h_dsigma_dt_B0->SetBinError(i, dsigma_dt_error);
    }

    for(int i = 1; i <= h_dsigma_dt_RP->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_RP->GetBinContent(i);
        double bin_error = h_dsigma_dt_RP->GetBinError(i);
        double bin_width = h_dsigma_dt_RP->GetBinWidth(i);

        double dsigma_dt = (bin_content * scale_factor) / bin_width;
        double dsigma_dt_error = (bin_error * scale_factor) / bin_width;

        h_dsigma_dt_RP->SetBinContent(i, dsigma_dt);
        h_dsigma_dt_RP->SetBinError(i, dsigma_dt_error);
    }
    h_dsigma_dt_Sum = (TH1D*)h_dsigma_dt_B0->Clone("dsigma_dt_Sum");
    h_dsigma_dt_Sum->SetDirectory(nullptr);
    h_dsigma_dt_Sum->SetTitle("B0+RP Reco d#sigma/dt (purity/eff-corrected);|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]");
    h_dsigma_dt_Sum->Add(h_dsigma_dt_RP);

    // Scale triple differential cross sections d^3sigma/(dQ^2 d#beta dx_{pom})
    for (int ix = 1; ix <= h_d3sigma_MC->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_d3sigma_MC->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h_d3sigma_MC->GetNbinsZ(); ++iz) {
                const double content = h_d3sigma_MC->GetBinContent(ix, iy, iz);
                const double err = h_d3sigma_MC->GetBinError(ix, iy, iz);
                const double dV = h_d3sigma_MC->GetXaxis()->GetBinWidth(ix) *
                                  h_d3sigma_MC->GetYaxis()->GetBinWidth(iy) *
                                  h_d3sigma_MC->GetZaxis()->GetBinWidth(iz);
                if (dV <= 0.0) continue;
                h_d3sigma_MC->SetBinContent(ix, iy, iz, (content * scale_factor) / dV);
                h_d3sigma_MC->SetBinError(ix, iy, iz, (err * scale_factor) / dV);
            }
        }
    }
    for (int ix = 1; ix <= h_d3sigma_B0->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_d3sigma_B0->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h_d3sigma_B0->GetNbinsZ(); ++iz) {
                const double content = h_d3sigma_B0->GetBinContent(ix, iy, iz);
                const double err = h_d3sigma_B0->GetBinError(ix, iy, iz);
                const double dV = h_d3sigma_B0->GetXaxis()->GetBinWidth(ix) *
                                  h_d3sigma_B0->GetYaxis()->GetBinWidth(iy) *
                                  h_d3sigma_B0->GetZaxis()->GetBinWidth(iz);
                if (dV <= 0.0) continue;
                h_d3sigma_B0->SetBinContent(ix, iy, iz, (content * scale_factor) / dV);
                h_d3sigma_B0->SetBinError(ix, iy, iz, (err * scale_factor) / dV);
            }
        }
    }
    for (int ix = 1; ix <= h_d3sigma_RP->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_d3sigma_RP->GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= h_d3sigma_RP->GetNbinsZ(); ++iz) {
                const double content = h_d3sigma_RP->GetBinContent(ix, iy, iz);
                const double err = h_d3sigma_RP->GetBinError(ix, iy, iz);
                const double dV = h_d3sigma_RP->GetXaxis()->GetBinWidth(ix) *
                                  h_d3sigma_RP->GetYaxis()->GetBinWidth(iy) *
                                  h_d3sigma_RP->GetZaxis()->GetBinWidth(iz);
                if (dV <= 0.0) continue;
                h_d3sigma_RP->SetBinContent(ix, iy, iz, (content * scale_factor) / dV);
                h_d3sigma_RP->SetBinError(ix, iy, iz, (err * scale_factor) / dV);
            }
        }
    }
    h_d3sigma_Sum = (TH3D*)h_d3sigma_B0->Clone("d3sigma_dQ2dbeta_dxpom_Sum");
    h_d3sigma_Sum->SetDirectory(nullptr);
    h_d3sigma_Sum->SetTitle("B0+RP Reco d^{3}#sigma/(dQ^{2}d#betadx_{pom}) (purity/eff-corrected);Q^{2} [GeV^{2}];#beta;x_{pom}");
    h_d3sigma_Sum->Add(h_d3sigma_RP);

    TH2D* h_Response_Q2_rowNorm = (TH2D*)h_Response_Q2->Clone("Response_Q2_rowNorm");
    TH2D* h_Response_Q2_colNorm = (TH2D*)h_Response_Q2->Clone("Response_Q2_colNorm");
    TH2D* h_Response_x_rowNorm  = (TH2D*)h_Response_x->Clone("Response_x_rowNorm");
    TH2D* h_Response_x_colNorm  = (TH2D*)h_Response_x->Clone("Response_x_colNorm");
    TH2D* h_Response_y_rowNorm  = (TH2D*)h_Response_y->Clone("Response_y_rowNorm");
    TH2D* h_Response_y_colNorm  = (TH2D*)h_Response_y->Clone("Response_y_colNorm");
    TH2D* h_Response_W2_EM_rowNorm = (TH2D*)h_Response_W2_EM->Clone("Response_W2_EM_rowNorm");
    TH2D* h_Response_W2_EM_colNorm = (TH2D*)h_Response_W2_EM->Clone("Response_W2_EM_colNorm");
    TH2D* h_Response_W2_DA_rowNorm = (TH2D*)h_Response_W2_DA->Clone("Response_W2_DA_rowNorm");
    TH2D* h_Response_W2_DA_colNorm = (TH2D*)h_Response_W2_DA->Clone("Response_W2_DA_colNorm");
    TH2D* h_Response_W2_Best_rowNorm = (TH2D*)h_Response_W2_Best->Clone("Response_W2_Best_rowNorm");
    TH2D* h_Response_W2_Best_colNorm = (TH2D*)h_Response_W2_Best->Clone("Response_W2_Best_colNorm");
    TH2D* h_Response_W2_Sigma_rowNorm = (TH2D*)h_Response_W2_Sigma->Clone("Response_W2_Sigma_rowNorm");
    TH2D* h_Response_W2_Sigma_colNorm = (TH2D*)h_Response_W2_Sigma->Clone("Response_W2_Sigma_colNorm");

    // Row normalization (normalize each truth bin -> efficiency/purity)
    for(int i = 1; i <= h_Response_Q2_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_Q2_rowNorm->GetNbinsY(); j++) rowSum += h_Response_Q2_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_Q2_rowNorm->GetNbinsY(); j++)
                h_Response_Q2_rowNorm->SetBinContent(i, j, h_Response_Q2_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }
    for(int i = 1; i <= h_Response_x_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_x_rowNorm->GetNbinsY(); j++) rowSum += h_Response_x_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_x_rowNorm->GetNbinsY(); j++)
                h_Response_x_rowNorm->SetBinContent(i, j, h_Response_x_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }
    for(int i = 1; i <= h_Response_y_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_y_rowNorm->GetNbinsY(); j++) rowSum += h_Response_y_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_y_rowNorm->GetNbinsY(); j++)
                h_Response_y_rowNorm->SetBinContent(i, j, h_Response_y_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }
    for(int i = 1; i <= h_Response_W2_EM_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_W2_EM_rowNorm->GetNbinsY(); j++) rowSum += h_Response_W2_EM_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_W2_EM_rowNorm->GetNbinsY(); j++)
                h_Response_W2_EM_rowNorm->SetBinContent(i, j, h_Response_W2_EM_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }
    for(int i = 1; i <= h_Response_W2_DA_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_W2_DA_rowNorm->GetNbinsY(); j++) rowSum += h_Response_W2_DA_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_W2_DA_rowNorm->GetNbinsY(); j++)
                h_Response_W2_DA_rowNorm->SetBinContent(i, j, h_Response_W2_DA_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }
    for(int i = 1; i <= h_Response_W2_Best_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_W2_Best_rowNorm->GetNbinsY(); j++) rowSum += h_Response_W2_Best_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_W2_Best_rowNorm->GetNbinsY(); j++)
                h_Response_W2_Best_rowNorm->SetBinContent(i, j, h_Response_W2_Best_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }
    for(int i = 1; i <= h_Response_W2_Sigma_rowNorm->GetNbinsX(); i++) {
        double rowSum = 0;
        for(int j = 1; j <= h_Response_W2_Sigma_rowNorm->GetNbinsY(); j++) rowSum += h_Response_W2_Sigma_rowNorm->GetBinContent(i, j);
        if(rowSum > 0) {
            for(int j = 1; j <= h_Response_W2_Sigma_rowNorm->GetNbinsY(); j++)
                h_Response_W2_Sigma_rowNorm->SetBinContent(i, j, h_Response_W2_Sigma_rowNorm->GetBinContent(i, j) / rowSum);
        }
    }

    // Column normalization (normalize each reco bin -> stability)
    for(int j = 1; j <= h_Response_Q2_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_Q2_colNorm->GetNbinsX(); i++) colSum += h_Response_Q2_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_Q2_colNorm->GetNbinsX(); i++)
                h_Response_Q2_colNorm->SetBinContent(i, j, h_Response_Q2_colNorm->GetBinContent(i, j) / colSum);
        }
    }
    for(int j = 1; j <= h_Response_x_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_x_colNorm->GetNbinsX(); i++) colSum += h_Response_x_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_x_colNorm->GetNbinsX(); i++)
                h_Response_x_colNorm->SetBinContent(i, j, h_Response_x_colNorm->GetBinContent(i, j) / colSum);
        }
    }
    for(int j = 1; j <= h_Response_y_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_y_colNorm->GetNbinsX(); i++) colSum += h_Response_y_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_y_colNorm->GetNbinsX(); i++)
                h_Response_y_colNorm->SetBinContent(i, j, h_Response_y_colNorm->GetBinContent(i, j) / colSum);
        }
    }
    for(int j = 1; j <= h_Response_W2_EM_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_W2_EM_colNorm->GetNbinsX(); i++) colSum += h_Response_W2_EM_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_W2_EM_colNorm->GetNbinsX(); i++)
                h_Response_W2_EM_colNorm->SetBinContent(i, j, h_Response_W2_EM_colNorm->GetBinContent(i, j) / colSum);
        }
    }
    for(int j = 1; j <= h_Response_W2_DA_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_W2_DA_colNorm->GetNbinsX(); i++) colSum += h_Response_W2_DA_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_W2_DA_colNorm->GetNbinsX(); i++)
                h_Response_W2_DA_colNorm->SetBinContent(i, j, h_Response_W2_DA_colNorm->GetBinContent(i, j) / colSum);
        }
    }
    for(int j = 1; j <= h_Response_W2_Best_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_W2_Best_colNorm->GetNbinsX(); i++) colSum += h_Response_W2_Best_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_W2_Best_colNorm->GetNbinsX(); i++)
                h_Response_W2_Best_colNorm->SetBinContent(i, j, h_Response_W2_Best_colNorm->GetBinContent(i, j) / colSum);
        }
    }
    for(int j = 1; j <= h_Response_W2_Sigma_colNorm->GetNbinsY(); j++) {
        double colSum = 0;
        for(int i = 1; i <= h_Response_W2_Sigma_colNorm->GetNbinsX(); i++) colSum += h_Response_W2_Sigma_colNorm->GetBinContent(i, j);
        if(colSum > 0) {
            for(int i = 1; i <= h_Response_W2_Sigma_colNorm->GetNbinsX(); i++)
                h_Response_W2_Sigma_colNorm->SetBinContent(i, j, h_Response_W2_Sigma_colNorm->GetBinContent(i, j) / colSum);
        }
    }

    // Build relative resolution vs global bin index histograms
    TH1D* h_Q2_RelRes_vs_k_EM = BuildRelResVsKHist(
        "Q2_RelRes_vs_k_EM",
        "Q^{2} relative resolution vs k (EM);k;RMS((Q^{2}_{reco}-Q^{2}_{truth})/Q^{2}_{truth})",
        res_Q2_EM_k);
    TH1D* h_Q2_RelRes_vs_k_DA = BuildRelResVsKHist(
        "Q2_RelRes_vs_k_DA",
        "Q^{2} relative resolution vs k (DA);k;RMS((Q^{2}_{reco}-Q^{2}_{truth})/Q^{2}_{truth})",
        res_Q2_DA_k);
    TH1D* h_Q2_RelRes_vs_k_Sigma = BuildRelResVsKHist(
        "Q2_RelRes_vs_k_Sigma",
        "Q^{2} relative resolution vs k (#Sigma);k;RMS((Q^{2}_{reco}-Q^{2}_{truth})/Q^{2}_{truth})",
        res_Q2_Sigma_k);

    TH1D* h_xpom_RelRes_vs_k_EM_B0 = BuildRelResVsKHist(
        "xpom_RelRes_vs_k_EM_B0",
        "x_{pom} relative resolution vs k (EM, B0);k;RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_EM_B0_k);
    TH1D* h_xpom_RelRes_vs_k_DA_B0 = BuildRelResVsKHist(
        "xpom_RelRes_vs_k_DA_B0",
        "x_{pom} relative resolution vs k (DA, B0);k;RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_DA_B0_k);
    TH1D* h_xpom_RelRes_vs_k_Sigma_B0 = BuildRelResVsKHist(
        "xpom_RelRes_vs_k_Sigma_B0",
        "x_{pom} relative resolution vs k (#Sigma, B0);k;RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_Sigma_B0_k);
    TH1D* h_xpom_RelRes_vs_k_EM_RP = BuildRelResVsKHist(
        "xpom_RelRes_vs_k_EM_RP",
        "x_{pom} relative resolution vs k (EM, RP);k;RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_EM_RP_k);
    TH1D* h_xpom_RelRes_vs_k_DA_RP = BuildRelResVsKHist(
        "xpom_RelRes_vs_k_DA_RP",
        "x_{pom} relative resolution vs k (DA, RP);k;RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_DA_RP_k);
    TH1D* h_xpom_RelRes_vs_k_Sigma_RP = BuildRelResVsKHist(
        "xpom_RelRes_vs_k_Sigma_RP",
        "x_{pom} relative resolution vs k (#Sigma, RP);k;RMS((x_{pom,reco}-x_{pom,truth})/x_{pom,truth})",
        res_xpom_Sigma_RP_k);

    TH1D* h_beta_RelRes_vs_k_EM_B0 = BuildRelResVsKHist(
        "beta_RelRes_vs_k_EM_B0",
        "#beta relative resolution vs k (EM, B0);k;RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_EM_B0_k);
    TH1D* h_beta_RelRes_vs_k_DA_B0 = BuildRelResVsKHist(
        "beta_RelRes_vs_k_DA_B0",
        "#beta relative resolution vs k (DA, B0);k;RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_DA_B0_k);
    TH1D* h_beta_RelRes_vs_k_Sigma_B0 = BuildRelResVsKHist(
        "beta_RelRes_vs_k_Sigma_B0",
        "#beta relative resolution vs k (#Sigma, B0);k;RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_Sigma_B0_k);
    TH1D* h_beta_RelRes_vs_k_EM_RP = BuildRelResVsKHist(
        "beta_RelRes_vs_k_EM_RP",
        "#beta relative resolution vs k (EM, RP);k;RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_EM_RP_k);
    TH1D* h_beta_RelRes_vs_k_DA_RP = BuildRelResVsKHist(
        "beta_RelRes_vs_k_DA_RP",
        "#beta relative resolution vs k (DA, RP);k;RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_DA_RP_k);
    TH1D* h_beta_RelRes_vs_k_Sigma_RP = BuildRelResVsKHist(
        "beta_RelRes_vs_k_Sigma_RP",
        "#beta relative resolution vs k (#Sigma, RP);k;RMS((#beta_{reco}-#beta_{truth})/#beta_{truth})",
        res_beta_Sigma_RP_k);

    // Write all histograms and TTree to the output file
    outputFile->cd();
    TParameter<Long64_t> nEventsParam("nEventsProcessed", nProcessed);
    nEventsParam.Write();
    TParameter<double> nGenParam("N_gen", N_gen);
    nGenParam.Write();
    TParameter<double> sigmaTotalParam("sigma_total_nb", sigma_total);
    sigmaTotalParam.Write();
    TParameter<double>("DISCut_q2_min", cutCfg.q2_min).Write();
    TParameter<double>("DISCut_q2_max", cutCfg.q2_max).Write();
    TParameter<double>("DISCut_y_min", cutCfg.y_min).Write();
    TParameter<double>("DISCut_y_max", cutCfg.y_max).Write();
    tree->Write();

    // Write Q2/xy histograms
    h_RelRes_Q2_EM->Write();
    h_RelRes_Q2_DA->Write();
    h_RelRes_Q2_Sigma->Write();
    h_RelRes_Q2_binned_EM->Write();
    h_RelRes_Q2_binned_DA->Write();
    h_RelRes_Q2_binned_Sigma->Write();
    h_Corr_Q2_EM->Write();
    h_Corr_Q2_DA->Write();
    h_Corr_Q2_Sigma->Write();
    h_Q2_RelRes_vs_k_EM->Write();
    h_Q2_RelRes_vs_k_DA->Write();
    h_Q2_RelRes_vs_k_Sigma->Write();
    h_Q2_truth->Write();
    h_Q2_EM->Write();
    h_Q2_DA->Write();
    h_Q2_Sigma->Write();
    
    h_x_truth->Write();  h_x_EM->Write();  h_x_DA->Write();  h_x_Sigma->Write();
    h_y_truth->Write();  h_y_EM->Write();  h_y_DA->Write();  h_y_Sigma->Write();
    h_RelRes_x_EM->Write();
    h_RelRes_x_DA->Write();
    h_RelRes_x_Sigma->Write();
    h_RelRes_y_EM->Write();
    h_RelRes_y_DA->Write();
    h_RelRes_y_Sigma->Write();
    h_RelRes_x_binned_EM->Write();
    h_RelRes_x_binned_DA->Write();
    h_RelRes_x_binned_Sigma->Write();
    h_RelRes_y_binned_EM->Write();
    h_RelRes_y_binned_DA->Write();
    h_RelRes_y_binned_Sigma->Write();
    h_Corr_x_EM->Write();
    h_Corr_x_DA->Write();
    h_Corr_x_Sigma->Write();
    h_Corr_y_EM->Write();
    h_Corr_y_DA->Write();
    h_Corr_y_Sigma->Write();

    // Backward-compatible aliases for legacy ESigma naming
    auto h_Q2_ESigma_alias = (TH1D*)h_Q2_Sigma->Clone("h_Q2_ESigma");
    auto h_x_ESigma_alias = (TH1D*)h_x_Sigma->Clone("x_ESigma");
    auto h_y_ESigma_alias = (TH1D*)h_y_Sigma->Clone("y_ESigma");
    auto h_RelRes_Q2_ESigma_alias = (TH1D*)h_RelRes_Q2_Sigma->Clone("Q2_RelRes_ESigma");
    auto h_Q2_Res_EM_alias = (TH1D*)h_RelRes_Q2_EM->Clone("Q2_Res_EM");
    auto h_RelRes_Q2_binned_ESigma_alias = (TH2D*)h_RelRes_Q2_binned_Sigma->Clone("Q2_RelRes_binned_ESigma");
    auto h_RelRes_x_ESigma_alias = (TH1D*)h_RelRes_x_Sigma->Clone("x_RelRes_ESigma");
    auto h_RelRes_x_binned_ESigma_alias = (TH2D*)h_RelRes_x_binned_Sigma->Clone("x_RelRes_binned_ESigma");
    auto h_RelRes_y_ESigma_alias = (TH1D*)h_RelRes_y_Sigma->Clone("y_RelRes_ESigma");
    auto h_RelRes_y_binned_ESigma_alias = (TH2D*)h_RelRes_y_binned_Sigma->Clone("y_RelRes_binned_ESigma");
    auto h_Corr_Q2_ESigma_alias = (TH2D*)h_Corr_Q2_Sigma->Clone("Corr_Q2_ESigma");
    auto h_Corr_x_ESigma_alias = (TH2D*)h_Corr_x_Sigma->Clone("x_Corr_ESigma");
    auto h_Corr_y_ESigma_alias = (TH2D*)h_Corr_y_Sigma->Clone("y_Corr_ESigma");
    auto h_MX2_alias = (TH1D*)h_MX2_reco->Clone("h_MX2");
    auto h_MX2_truth_alias = (TH1D*)h_MX2_truth->Clone("h_MX2_truth");
    h_Q2_ESigma_alias->Write();
    h_x_ESigma_alias->Write();
    h_y_ESigma_alias->Write();
    h_RelRes_Q2_ESigma_alias->Write();
    h_Q2_Res_EM_alias->Write();
    h_RelRes_Q2_binned_ESigma_alias->Write();
    h_RelRes_x_ESigma_alias->Write();
    h_RelRes_x_binned_ESigma_alias->Write();
    h_RelRes_y_ESigma_alias->Write();
    h_RelRes_y_binned_ESigma_alias->Write();
    h_Corr_Q2_ESigma_alias->Write();
    h_Corr_x_ESigma_alias->Write();
    h_Corr_y_ESigma_alias->Write();
    h_MX2_alias->Write();
    h_MX2_truth_alias->Write();
    h_xQ2_truth->Write();
    h_xQ2_reco->Write();

    // Write Q2 resolution profile histograms (binned in x, Q2)
    h_Q2_RelRes_vs_xy_EM->Write();
    h_Q2_RelRes_vs_xy_DA->Write();
    h_Q2_RelRes_vs_xy_Sigma->Write();

    // Write x resolution profile histograms (binned in x, Q2)
    h_x_RelRes_vs_xQ2_EM->Write();
    h_x_RelRes_vs_xQ2_DA->Write();
    h_x_RelRes_vs_xQ2_Sigma->Write();

    // Write y resolution profile histograms (binned in x, Q2)
    h_y_RelRes_vs_xQ2_EM->Write();
    h_y_RelRes_vs_xQ2_DA->Write();
    h_y_RelRes_vs_xQ2_Sigma->Write();
    h_EPz_truth->Write();
    h_EPz->Write();
    h_EPz_2D->Write();
    h_eta_max->Write();
    h_eta_max_truth->Write();

    // Write NEW inclusive DIS histograms
    h_W2_EM->Write(); h_W2_DA->Write(); h_W2_Best->Write(); h_W2_Sigma->Write(); h_W2_truth->Write();
    h_RelRes_W2_EM->Write(); h_RelRes_W2_DA->Write(); h_RelRes_W2_Best->Write(); h_RelRes_W2_Sigma->Write();
    h_RelRes_W2_binned_EM->Write(); h_RelRes_W2_binned_DA->Write(); h_RelRes_W2_binned_Best->Write(); h_RelRes_W2_binned_Sigma->Write();
    h_Ep_e->Write(); h_phi_e->Write(); h_pT_e->Write();
    h_Ep_e_truth->Write(); h_phi_e_truth->Write(); h_pT_e_truth->Write();
    h_Response_Q2->Write(); h_Response_x->Write(); h_Response_y->Write();
    h_Response_W2_EM->Write(); h_Response_W2_DA->Write(); h_Response_W2_Best->Write(); h_Response_W2_Sigma->Write();
    h_Response_Q2_rowNorm->Write(); h_Response_Q2_colNorm->Write();
    h_Response_x_rowNorm->Write(); h_Response_x_colNorm->Write();
    h_Response_y_rowNorm->Write(); h_Response_y_colNorm->Write();
    h_Response_W2_EM_rowNorm->Write(); h_Response_W2_EM_colNorm->Write();
    h_Response_W2_DA_rowNorm->Write(); h_Response_W2_DA_colNorm->Write();
    h_Response_W2_Best_rowNorm->Write(); h_Response_W2_Best_colNorm->Write();
    h_Response_W2_Sigma_rowNorm->Write(); h_Response_W2_Sigma_colNorm->Write();

    g_Q2_EM->Write(); g_Q2_DA->Write(); g_Q2_Sigma->Write();
    g_x_EM->Write(); g_x_DA->Write(); g_x_Sigma->Write();
    g_y_EM->Write(); g_y_DA->Write(); g_y_Sigma->Write();
    g_W2_EM->Write(); g_W2_DA->Write(); g_W2_Best->Write(); g_W2_Sigma->Write();
    g_Ep_e->Write(); g_phi_e->Write(); g_pT_e->Write();

    TParameter<double>("W2Best_W0", kW2BestW0).Write();
    TParameter<double>("W2Best_Gamma", kW2BestGamma).Write();

    // Write diffractive t histograms/graphs
    h_t_MC->Write();
    h_t_B0->Write();
    h_t_RP_histo->Write();
    h_dsigma_dt_MC->Write();
    h_dsigma_dt_B0->Write();
    h_dsigma_dt_RP->Write();
    if (h_dsigma_dt_Sum) h_dsigma_dt_Sum->Write();
    h_theta_MC->Write();
    h_theta_all_TS->Write();
    h_theta_B0->Write();
    h_theta_RP->Write();
    h_t_res_B0->Write();
    h_t_res_RP->Write();
    h_t_RelRes_binned_B0->Write();
    h_t_RelRes_binned_RP->Write();
    h_t_corr_B0->Write();
    h_t_corr_RP->Write();
    g_t_B0->Write();
    g_t_RP->Write();
    h_xL_MC->Write();
    h_xL_B0->Write();
    h_xL_RP->Write();
    h_xL_res_B0->Write();
    h_xL_res_RP->Write();
    h_xL_corr_B0->Write();
    h_xL_corr_RP->Write();
    h_xL_RelRes_binned_B0->Write();
    h_xL_RelRes_binned_RP->Write();
    g_xL_B0->Write();
    g_xL_RP->Write();
    h_xpom_truth_all->Write();
    h_xpom_truth_B0->Write();
    h_xpom_truth_RP->Write();
    h_xpom_reco_EM_all->Write();
    h_xpom_reco_DA_all->Write();
    h_xpom_reco_Sigma_all->Write();
    h_xpom_reco_EM_B0->Write();
    h_xpom_reco_EM_RP->Write();
    h_xpom_reco_DA_B0->Write();
    h_xpom_reco_DA_RP->Write();
    h_xpom_reco_Sigma_B0->Write();
    h_xpom_reco_Sigma_RP->Write();
    h_xpom_MC->Write();
    h_xpom_B0->Write();
    h_xpom_RP->Write();
    h_xpom_def_MC->Write();
    h_xpom_def_B0->Write();
    h_xpom_def_RP->Write();
    h_xpom_comp_MC->Write();
    h_xpom_comp_B0->Write();
    h_xpom_comp_RP->Write();
    h_xpom_corr_B0->Write();
    h_xpom_corr_RP->Write();
    h_xpom_res_B0->Write();
    h_xpom_res_RP->Write();
    h_xpom_RelRes_binned_B0->Write();
    h_xpom_RelRes_binned_RP->Write();
    h_xpom_RelRes_binned_W2Best_B0->Write();
    h_xpom_RelRes_binned_W2Best_RP->Write();
    h_Response_xpom_EM_B0->Write();
    h_Response_xpom_EM_RP->Write();
    h_Response_xpom_DA_B0->Write();
    h_Response_xpom_DA_RP->Write();
    h_Response_xpom_Sigma_B0->Write();
    h_Response_xpom_Sigma_RP->Write();
    h_Response_xpom_W2Best_B0->Write();
    h_Response_xpom_W2Best_RP->Write();
    h_Response_beta_EM_B0->Write();
    h_Response_beta_EM_RP->Write();
    h_Response_beta_DA_B0->Write();
    h_Response_beta_DA_RP->Write();
    h_Response_beta_Sigma_B0->Write();
    h_Response_beta_Sigma_RP->Write();
    h_Response_beta_W2Best_B0->Write();
    h_Response_beta_W2Best_RP->Write();
    h_xpom_RelRes_vs_k_EM_B0->Write();
    h_xpom_RelRes_vs_k_DA_B0->Write();
    h_xpom_RelRes_vs_k_Sigma_B0->Write();
    h_xpom_RelRes_vs_k_EM_RP->Write();
    h_xpom_RelRes_vs_k_DA_RP->Write();
    h_xpom_RelRes_vs_k_Sigma_RP->Write();
    h_beta_RelRes_vs_k_EM_B0->Write();
    h_beta_RelRes_vs_k_DA_B0->Write();
    h_beta_RelRes_vs_k_Sigma_B0->Write();
    h_beta_RelRes_vs_k_EM_RP->Write();
    h_beta_RelRes_vs_k_DA_RP->Write();
    h_beta_RelRes_vs_k_Sigma_RP->Write();
    h_beta_truth_all->Write();
    h_beta_truth_B0->Write();
    h_beta_truth_RP->Write();
    h_beta_reco_EM_all->Write();
    h_beta_reco_DA_all->Write();
    h_beta_reco_Sigma_all->Write();
    h_beta_reco_EM_B0->Write();
    h_beta_reco_EM_RP->Write();
    h_beta_reco_DA_B0->Write();
    h_beta_reco_DA_RP->Write();
    h_beta_reco_Sigma_B0->Write();
    h_beta_reco_Sigma_RP->Write();
    h_beta_MC->Write();
    h_beta_B0->Write();
    h_beta_RP->Write();
    h_beta_res_B0->Write();
    h_beta_res_RP->Write();
    h_beta_corr_B0->Write();
    h_beta_corr_RP->Write();
    h_beta_RelRes_binned_B0->Write();
    h_beta_RelRes_binned_RP->Write();
    h_beta_RelRes_binned_W2Best_B0->Write();
    h_beta_RelRes_binned_W2Best_RP->Write();
    h_beta_vs_Q2->Write();
    h_beta_vs_xpom->Write();
    h_beta_vs_t->Write();
    h_t_truth_mc->Write();
    h_t_reco_mc->Write();
    h_t_reco_pdata->Write();
    h_t_truth_pdata_all->Write();
    h_beta_truth_mc->Write();
    h_beta_reco_mc->Write();
    h_beta_reco_pdata->Write();
    h_beta_truth_pdata_all->Write();
    h_xpom_truth_mc->Write();
    h_xpom_reco_mc->Write();
    h_xpom_reco_pdata->Write();
    h_xpom_truth_pdata_all->Write();
    h_Q2_truth_mc->Write();
    h_Q2_reco_mc->Write();
    h_Q2_reco_pdata->Write();
    h_Q2_truth_pdata->Write();
    h_MX2_truth_mc->Write();
    h_MX2_reco_mc->Write();
    h_MX2_reco_pdata->Write();
    h_MX2_truth_pdata->Write();
    h_x_truth_mc->Write();
    h_x_reco_mc->Write();
    h_x_reco_pdata->Write();
    h_x_truth_pdata->Write();
    h_y_truth_mc->Write();
    h_y_reco_mc->Write();
    h_y_reco_pdata->Write();
    h_y_truth_pdata->Write();
    h_W2_truth_mc->Write();
    h_W2_reco_mc->Write();
    h_W2_reco_pdata->Write();
    h_W2_truth_pdata->Write();
    h_t_truth_mc_B0->Write();
    h_t_reco_mc_B0->Write();
    h_t_same_mc_B0->Write();
    h_t_reco_pdata_B0->Write();
    h_t_truth_pdata_B0->Write();
    h_t_truth_mc_RP->Write();
    h_t_reco_mc_RP->Write();
    h_t_same_mc_RP->Write();
    h_t_reco_pdata_RP->Write();
    h_t_truth_pdata_RP->Write();
    h_beta_truth_mc_B0->Write();
    h_beta_reco_mc_B0->Write();
    h_beta_reco_pdata_B0->Write();
    h_beta_truth_pdata_B0->Write();
    h_beta_truth_mc_RP->Write();
    h_beta_reco_mc_RP->Write();
    h_beta_reco_pdata_RP->Write();
    h_beta_truth_pdata_RP->Write();
    h_xpom_truth_mc_B0->Write();
    h_xpom_reco_mc_B0->Write();
    h_xpom_reco_pdata_B0->Write();
    h_xpom_truth_pdata_B0->Write();
    h_xpom_truth_mc_RP->Write();
    h_xpom_reco_mc_RP->Write();
    h_xpom_reco_pdata_RP->Write();
    h_xpom_truth_pdata_RP->Write();
    h_theta_truth_mc_B0->Write();
    h_theta_reco_mc_B0->Write();
    h_theta_reco_pdata_B0->Write();
    h_theta_truth_pdata_B0->Write();
    h_theta_truth_mc_RP->Write();
    h_theta_reco_mc_RP->Write();
    h_theta_reco_pdata_RP->Write();
    h_theta_truth_pdata_RP->Write();
    h_theta_truth_pdata_all->Write();
    h_t_RelRes_vs_xpomQ2_B0->Write();
    h_t_RelRes_vs_xpomQ2_RP->Write();
    h_xpom_RelRes_vs_xpomQ2_B0->Write();
    h_xpom_RelRes_vs_xpomQ2_RP->Write();
    h_beta_RelRes_vs_betaQ2_B0->Write();
    h_beta_RelRes_vs_betaQ2_RP->Write();
    h_xL_RelRes_vs_xLQ2_B0->Write();
    h_xL_RelRes_vs_xLQ2_RP->Write();
    h_d3sigma_MC->Write();
    h_d3sigma_B0->Write();
    h_d3sigma_RP->Write();
    h_d3sigma_same_mc_B0->Write();
    h_d3sigma_same_mc_RP->Write();
    if (h_d3sigma_Sum) h_d3sigma_Sum->Write();
    h_beta_Q2_truth->Write();
    h_t_Q2_truth->Write();
    h_xpom_Q2_truth->Write();
    h_beta_t_truth->Write();
    h_xbj_t_truth->Write();
    h_xpom_t_truth->Write();
    h_xpom_beta_truth->Write();
    h_xbj_beta_truth->Write();
    h_xbj_xpom_truth->Write();
    h_beta_Q2_reco->Write();
    h_t_Q2_reco->Write();
    h_xpom_Q2_reco->Write();
    h_beta_t_reco->Write();
    h_xbj_t_reco->Write();
    h_xpom_t_reco->Write();
    h_xpom_beta_reco->Write();
    h_xbj_beta_reco->Write();
    h_xbj_xpom_reco->Write();
    h_phase3D_reco->Write();
    if (h_phase3D_reco_yaml) {
        h_phase3D_reco_yaml->Write();
    }
    h_phase_bin_gen->Write();
    h_phase_bin_meas->Write();
    h_phase_bin_gen_meas_same->Write();
    h_phase_bin_gen_setA->Write();
    h_phase_bin_meas_setA->Write();
    h_phase_bin_gen_meas_same_setA->Write();
    h_phase_bin_gen_setB->Write();
    h_phase_bin_meas_setB->Write();
    h_phase_bin_gen_meas_same_setB->Write();
    g_xpom_EM_B0->Write();
    g_xpom_EM_RP->Write();
    g_xpom_DA_B0->Write();
    g_xpom_DA_RP->Write();
    g_xpom_Sigma_B0->Write();
    g_xpom_Sigma_RP->Write();
    g_beta_EM_B0->Write();
    g_beta_EM_RP->Write();
    g_beta_DA_B0->Write();
    g_beta_DA_RP->Write();
    g_beta_Sigma_B0->Write();
    g_beta_Sigma_RP->Write();
    g_MX2->Write();

    // Write diffractive M_X^2 histograms/correlations
    h_MX2_truth->Write();
    h_MX2_reco->Write();
    h_MX2_corr->Write();
    h_MX2_RelRes->Write();
    h_MX2_RelRes_binned->Write();
    h_MX2_RelRes_vs_MX2Q2->Write();
    h_MX2_t_truth->Write();
    h_MX2_t_B0->Write();
    h_MX2_t_RP->Write();

    bool wroteOccupancyYAML = WriteBinOccupancyYAML("/data/bin_occupancy.yaml", yaml_bins, occupancy_truth_k);
    std::string occupancyPath = "/data/bin_occupancy.yaml";
    if (!wroteOccupancyYAML) {
        wroteOccupancyYAML = WriteBinOccupancyYAML("data/bin_occupancy.yaml", yaml_bins, occupancy_truth_k);
        occupancyPath = "data/bin_occupancy.yaml";
    }
    if (!wroteOccupancyYAML) {
        wroteOccupancyYAML = WriteBinOccupancyYAML("bin_occupancy.yaml", yaml_bins, occupancy_truth_k);
        occupancyPath = "bin_occupancy.yaml";
    }
    if (wroteOccupancyYAML) {
        std::cout << "Wrote bin occupancy YAML: " << occupancyPath << std::endl;
    } else {
        std::cerr << "WARNING: Failed to write bin occupancy YAML." << std::endl;
    }

    outputFile->Close();
    delete events;
    delete outputFile;

    std::cout << "\nInclusive kinematics summary:" << std::endl;
    std::cout << "  DIS cutflow config: Q2>(" << cutCfg.q2_min << "), y in (" << cutCfg.y_min
              << ", " << cutCfg.y_max << "), E-pz in (" << cutCfg.epz_min << ", " << cutCfg.epz_max
              << "), apply_epz=" << (cutCfg.apply_epz_cut ? "true" : "false")
              << ", proton_acceptance=" << (cutCfg.apply_proton_tag_cuts ? "true" : "false") << std::endl;
    std::cout << "  Passed DIS cuts (truth/EM/DA/Sigma): "
              << nPassTruthDIS << " / " << nPassEMDIS << " / " << nPassDADIS;
    if (sigmaAvailable) std::cout << " / " << nPassSigmaDIS;
    std::cout << std::endl;
    std::cout << "  Valid truth events (Q2,x,y): " << nTruthValid << " / " << nentries << std::endl;
    std::cout << "  Valid EM events   (Q2,x,y): " << nEMValid << " / " << nentries << std::endl;
    std::cout << "  Valid DA events   (Q2,x,y): " << nDAValid << " / " << nentries << std::endl;
    if (sigmaAvailable) {
        std::cout << "  Valid Sigma events(Q2,x,y): " << nSigmaValid << " / " << nentries << std::endl;
    }
    std::cout << "  Truth kinematics derived from MC: " << nTruthFromMC << std::endl;
    std::cout << "  Out-of-range truth counts: Q2=" << nOutQ2Truth
              << " x=" << nOutXTruth << " y=" << nOutYTruth << std::endl;
    std::cout << "  Out-of-range EM counts:    Q2=" << nOutQ2EM
              << " x=" << nOutXEM << " y=" << nOutYEM << std::endl;
    std::cout << "  Out-of-range DA counts:    Q2=" << nOutQ2DA
              << " x=" << nOutXDA << " y=" << nOutYDA << std::endl;
    if (sigmaAvailable) {
        std::cout << "  Out-of-range Sigma counts: Q2=" << nOutQ2Sigma
                  << " x=" << nOutXSigma << " y=" << nOutYSigma << std::endl;
    }
    
    std::cout << "\nAnalysis complete! Output saved to DDIS_Combined_output.root" << std::endl;
    return 0;
}
