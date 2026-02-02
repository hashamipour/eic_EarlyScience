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
#include <algorithm>
#include <random>
#include <cctype>
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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fileList.txt> [N]" << std::endl;
        return 1;
    }
    TString fileList = argv[1];
    int sampleN = -1;
    if (argc >= 3) {
        char* end = nullptr;
        long parsed = std::strtol(argv[2], &end, 10);
        if (!end || *end != '\0' || parsed <= 0) {
            std::cerr << "Error: N must be a positive integer." << std::endl;
            return 1;
        }
        sampleN = static_cast<int>(parsed);
    }

    std::cout<< " __ __ __ __ __ __ __ __ __ __" <<std::endl;
    std::cout<< "|                             |"<<std::endl;
    std::cout<< "|   Inclusive DIS Skim        |"<<std::endl;
    std::cout<< "|   Q2/x/y/W2 + e- vars       |"<<std::endl;
    std::cout<< "|__ __ __ __ __ __ __ __ __ __|"<<std::endl;
    std::cout<< "\nInput filelist: " << fileList <<std::endl;

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
    for (const auto& fileName : selectedFiles) {
        TString tmp = fileName;
        // Skip existence check for XRootD - TFile::Open handles remote files
        if(tmp.BeginsWith("root://")) {
            events->Add(tmp);
            nFiles++;
            continue;
        }
        // Local file check
        if (!std::filesystem::exists(tmp.Data())) {
            std::cerr << "\nWarning: File does not exist: " << fileName << std::endl;
            continue;
        }
        events->Add(tmp);
        nFiles++;
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
    std::vector<Double_t> bin_edges_Q2 = GetRoundedLogBins(3.4, 150.0, n_bins);
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

    // Relative resolution histograms
    TH1D* h_RelRes_x_EM     = new TH1D("x_RelRes_EM",     "electron method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_x_DA     = new TH1D("x_RelRes_DA",     "DA method;#frac{x(Reco)-x(MC)}{x(MC)}",       101, -0.15, 0.15);
    TH1D* h_RelRes_x_Sigma = new TH1D("x_RelRes_Sigma", "Sigma method;#frac{x(Reco)-x(MC)}{x(MC)}",   101, -0.15, 0.15);

    TH1D* h_RelRes_y_EM     = new TH1D("y_RelRes_EM",     "electron method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_y_DA     = new TH1D("y_RelRes_DA",     "DA method;#frac{y(Reco)-y(MC)}{y(MC)}",       101, -0.15, 0.15);
    TH1D* h_RelRes_y_Sigma = new TH1D("y_RelRes_Sigma", "Sigma method;#frac{y(Reco)-y(MC)}{y(MC)}",   101, -0.15, 0.15);

    // 2D binned relres vs truth
    int n_binned = 51;
    TH2D* h_RelRes_x_binned_EM     = new TH2D("x_RelRes_binned_EM",     "x: truth vs rel. res (EM);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_DA     = new TH2D("x_RelRes_binned_DA",     "x: truth vs rel. res (DA);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",       
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_Sigma = new TH2D("x_RelRes_binned_Sigma", "x: truth vs rel. res (Sigma);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",   
        x_bins.size()-1, x_bins.data(), n_binned, -0.15, 0.15);

    TH2D* h_RelRes_y_binned_EM     = new TH2D("y_RelRes_binned_EM",     "y: truth vs rel. res (EM);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_DA     = new TH2D("y_RelRes_binned_DA",     "y: truth vs rel. res (DA);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",       n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_Sigma = new TH2D("y_RelRes_binned_Sigma", "y: truth vs rel. res (Sigma);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",   n_y_bins, 0.0, 1.0, n_binned, -0.15, 0.15);

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

    TH1D* h_RelRes_Q2_EM = new TH1D("Q2_RelRes_EM","electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    TH1D* h_RelRes_Q2_DA = new TH1D("Q2_RelRes_DA","DA method;#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    TH1D* h_RelRes_Q2_Sigma = new TH1D("Q2_RelRes_Sigma","Sigma method;#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    
    TH2D* h_RelRes_Q2_binned_EM = new TH2D("Q2_RelRes_binned_EM",";Q^{2} [GeV^{2}];#frac{Q^{2}(EM)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_DA = new TH2D("Q2_RelRes_binned_DA",";Q^{2} [GeV^{2}];#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_Sigma = new TH2D("Q2_RelRes_binned_Sigma",";Q^{2} [GeV^{2}];#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);

    
    TH1D* h_Q2_truth    = new TH1D("h_Q2_truth","Q^2;# of events",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_EM       = new TH1D("h_Q2_EM",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_DA       = new TH1D("h_Q2_DA",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_Sigma   = new TH1D("h_Q2_Sigma",";Q^{2}",n_bins, bin_edges_Q2.data());

    //---------------------------------------------------------
    // NEW INCLUSIVE DIS HISTOGRAMS
    //---------------------------------------------------------

    // W² histograms
    const int n_w2_bins = 60;
    const double w2_min = 10.0;
    const double w2_max = 1.0e4;
    std::vector<Double_t> w2_bins = GetLogBins(w2_min, w2_max, n_w2_bins);
    TH1D* h_W2_EM    = new TH1D("W2_EM",    "Electron method;W^{2} [GeV^{2}]", w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_DA    = new TH1D("W2_DA",    "DA method;W^{2} [GeV^{2}]",       w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_Sigma = new TH1D("W2_Sigma", "Sigma method;W^{2} [GeV^{2}]",    w2_bins.size()-1, w2_bins.data());
    TH1D* h_W2_truth = new TH1D("W2_truth", "Truth;W^{2} [GeV^{2}]",           w2_bins.size()-1, w2_bins.data());

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
    TGraph* g_W2_Sigma = new TGraph(); g_W2_Sigma->SetName("g_W2_Sigma");

    TGraph* g_Ep_e = new TGraph(); g_Ep_e->SetName("g_Ep_e");
    TGraph* g_phi_e = new TGraph(); g_phi_e->SetName("g_phi_e");
    TGraph* g_pT_e = new TGraph(); g_pT_e->SetName("g_pT_e");

    int n_g_Q2_EM = 0, n_g_Q2_DA = 0, n_g_Q2_Sigma = 0;
    int n_g_x_EM = 0, n_g_x_DA = 0, n_g_x_Sigma = 0;
    int n_g_y_EM = 0, n_g_y_DA = 0, n_g_y_Sigma = 0;
    int n_g_W2_EM = 0, n_g_W2_DA = 0, n_g_W2_Sigma = 0;
    int n_g_Ep_e = 0, n_g_phi_e = 0, n_g_pT_e = 0;

    //---------------------------------------------------------
    // DIFFRACTIVE: Mandelstam t
    //---------------------------------------------------------
    std::vector<Double_t> t_bins_low = GetLogBins(1e-3, 0.5, 20);
    std::vector<Double_t> t_bins_high = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.6, 1.8, 2.0};
    std::vector<Double_t> t_bins = t_bins_low;
    for(size_t i = 1; i < t_bins_high.size(); i++){
        t_bins.push_back(t_bins_high[i]);
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

    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_RP = new TH1D("theta_RP", "RP Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);

    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_B0 = new TH1D("xL_B0", "B0 Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_RP = new TH1D("xL_RP", "RP Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);

    // x_pom histograms from definition: x_pom = (Q^2 + M_X^2 - t)/(Q^2 + W^2 - m_p^2)
    std::vector<Double_t> xpom_bins = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1};
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

    TH1D* h_beta_truth_all = new TH1D("beta_truth_all", "Truth #beta (All);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_truth_B0 = new TH1D("beta_truth_B0", "Truth #beta (B0);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_truth_RP = new TH1D("beta_truth_RP", "Truth #beta (RP);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_EM_all = new TH1D("beta_reco_EM_all", "Reco #beta EM (All);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_EM_B0 = new TH1D("beta_reco_EM_B0", "Reco #beta EM (B0);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_EM_RP = new TH1D("beta_reco_EM_RP", "Reco #beta EM (RP);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_DA_all = new TH1D("beta_reco_DA_all", "Reco #beta DA (All);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_DA_B0 = new TH1D("beta_reco_DA_B0", "Reco #beta DA (B0);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_DA_RP = new TH1D("beta_reco_DA_RP", "Reco #beta DA (RP);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_Sigma_all = new TH1D("beta_reco_Sigma_all", "Reco #beta #Sigma (All);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_Sigma_B0 = new TH1D("beta_reco_Sigma_B0", "Reco #beta #Sigma (B0);#beta;Counts", 50, 0.0, 1.0);
    TH1D* h_beta_reco_Sigma_RP = new TH1D("beta_reco_Sigma_RP", "Reco #beta #Sigma (RP);#beta;Counts", 50, 0.0, 1.0);

    // Phase-space density plots (truth, unbinned-style)
    const int n_beta_bins_density = 120;
    const int n_xpom_bins_density = 120;
    std::vector<Double_t> xpom_bins_density = GetLogBins(1.0e-4, 0.3, n_xpom_bins_density);
    TH2D* h_beta_Q2_truth = new TH2D("beta_Q2_truth", "Event Density;#beta;Q^{2} [GeV^{2}]",
                                    n_beta_bins_density, 0.0, 1.0,
                                    q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_t_Q2_truth = new TH2D("t_Q2_truth", "Event Density;|t| [GeV^{2}];Q^{2} [GeV^{2}]",
                                 t_bins.size()-1, t_bins.data(),
                                 q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_xpom_Q2_truth = new TH2D("xpom_Q2_truth", "Event Density;x_{pom};Q^{2} [GeV^{2}]",
                                    xpom_bins_density.size()-1, xpom_bins_density.data(),
                                    q2_bins_density.size()-1, q2_bins_density.data());
    TH2D* h_beta_t_truth = new TH2D("beta_t_truth", "Event Density;#beta;|t| [GeV^{2}]",
                                   n_beta_bins_density, 0.0, 1.0,
                                   t_bins.size()-1, t_bins.data());
    TH2D* h_xbj_t_truth = new TH2D("xbj_t_truth", "Event Density;x_{Bj};|t| [GeV^{2}]",
                                  x_bins_density.size()-1, x_bins_density.data(),
                                  t_bins.size()-1, t_bins.data());
    TH2D* h_xpom_t_truth = new TH2D("xpom_t_truth", "Event Density;x_{pom};|t| [GeV^{2}]",
                                   xpom_bins_density.size()-1, xpom_bins_density.data(),
                                   t_bins.size()-1, t_bins.data());
    TH2D* h_xpom_beta_truth = new TH2D("xpom_beta_truth", "Event Density;x_{pom};#beta",
                                      xpom_bins_density.size()-1, xpom_bins_density.data(),
                                      n_beta_bins_density, 0.0, 1.0);
    TH2D* h_xbj_beta_truth = new TH2D("xbj_beta_truth", "Event Density;x_{Bj};#beta",
                                     x_bins_density.size()-1, x_bins_density.data(),
                                     n_beta_bins_density, 0.0, 1.0);
    TH2D* h_xbj_xpom_truth = new TH2D("xbj_xpom_truth", "Event Density;x_{Bj};x_{pom}",
                                     x_bins_density.size()-1, x_bins_density.data(),
                                     xpom_bins_density.size()-1, xpom_bins_density.data());

    TH1D* h_t_res_B0 = new TH1D("t_res_B0", "B0 t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 100, -2.0, 2.0);
    TH1D* h_t_res_RP = new TH1D("t_res_RP", "RP t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 100, -2.0, 2.0);

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
    std::cout << "Finding beam particles..." << std::endl;
    
    BeamInfo beams;
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
    for (Long64_t i = 0; i < nentries; i++) {
        // Update the counter every 1000 events
        if (i % 1000 == 0) {
        printf("\rProcessing event %lld of %lld; %.2f percent done.", i, nentries, 100.0*i/nentries);
        fflush(stdout);
        }
        if (!tree_reader.Next()) break;
        nProcessed++;

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

        const bool validTruthQ2 = isValidQ2(electron_Q2_truth);
        const bool validTruthX = isValidX(electron_x_truth);
        const bool validTruthY = isValidY(electron_y_truth);
        if (!validTruthQ2) nOutQ2Truth++;
        if (!validTruthX) nOutXTruth++;
        if (!validTruthY) nOutYTruth++;
        if (validTruthQ2 && validTruthX && validTruthY) nTruthValid++;

        const bool validEMQ2 = isValidQ2(electron_Q2_EM);
        const bool validEMX = isValidX(electron_x_EM);
        const bool validEMY = isValidY(electron_y_EM);
        if (!validEMQ2) nOutQ2EM++;
        if (!validEMX) nOutXEM++;
        if (!validEMY) nOutYEM++;
        if (validEMQ2 && validEMX && validEMY) nEMValid++;

        const bool validDAQ2 = isValidQ2(electron_Q2_DA);
        const bool validDAX = isValidX(electron_x_DA);
        const bool validDAY = isValidY(electron_y_DA);
        if (!validDAQ2) nOutQ2DA++;
        if (!validDAX) nOutXDA++;
        if (!validDAY) nOutYDA++;
        if (validDAQ2 && validDAX && validDAY) nDAValid++;

        const bool validSigmaQ2 = isValidQ2(electron_Q2_Sigma);
        const bool validSigmaX = isValidX(electron_x_Sigma);
        const bool validSigmaY = isValidY(electron_y_Sigma);
        if (hasSigmaMethod) {
            if (!validSigmaQ2) nOutQ2Sigma++;
            if (!validSigmaX) nOutXSigma++;
            if (!validSigmaY) nOutYSigma++;
            if (validSigmaQ2 && validSigmaX && validSigmaY) nSigmaValid++;
        }

        if (validTruthX) h_x_truth->Fill(electron_x_truth);
        if (validEMX) h_x_EM->Fill(electron_x_EM);
        if (validDAX) h_x_DA->Fill(electron_x_DA);
        if (hasSigmaMethod && validSigmaX) h_x_Sigma->Fill(electron_x_Sigma);

        if (validTruthY) h_y_truth->Fill(electron_y_truth);
        if (validEMY) h_y_EM->Fill(electron_y_EM);
        if (validDAY) h_y_DA->Fill(electron_y_DA);
        if (hasSigmaMethod && validSigmaY) h_y_Sigma->Fill(electron_y_Sigma);

        // Relative resolutions vs truth
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

        // Correlations
        if (validTruthX && validEMX) g_x_EM->SetPoint(n_g_x_EM++, electron_x_truth, electron_x_EM);
        if (validTruthX && validDAX) g_x_DA->SetPoint(n_g_x_DA++, electron_x_truth, electron_x_DA);
        if (hasSigmaMethod && validTruthX && validSigmaX) g_x_Sigma->SetPoint(n_g_x_Sigma++, electron_x_truth, electron_x_Sigma);

        if (validTruthY && validEMY) g_y_EM->SetPoint(n_g_y_EM++, electron_y_truth, electron_y_EM);
        if (validTruthY && validDAY) g_y_DA->SetPoint(n_g_y_DA++, electron_y_truth, electron_y_DA);
        if (hasSigmaMethod && validTruthY && validSigmaY) g_y_Sigma->SetPoint(n_g_y_Sigma++, electron_y_truth, electron_y_Sigma);

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

        if (validTruthQ2 && validEMQ2) {
            h_RelRes_Q2_EM->Fill((electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth);
            h_RelRes_Q2_binned_EM->Fill(electron_Q2_truth, (electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth);
        }
        if (validTruthQ2 && validDAQ2) {
            h_RelRes_Q2_DA->Fill((electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth);
            h_RelRes_Q2_binned_DA->Fill(electron_Q2_truth, (electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth);
        }
        if (hasSigmaMethod && validTruthQ2 && validSigmaQ2) {
            h_RelRes_Q2_Sigma->Fill((electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth);
            h_RelRes_Q2_binned_Sigma->Fill(electron_Q2_truth, (electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth);
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

        //=================================================================
        // NEW INCLUSIVE DIS FILLS
        //=================================================================

        // Fill W² histograms
        const bool validWTruth = std::isfinite(electron_W_truth) && electron_W_truth > 0.0f;
        const bool validWEM = std::isfinite(electron_W_EM) && electron_W_EM > 0.0f;
        const bool validWDA = std::isfinite(electron_W_DA) && electron_W_DA > 0.0f;
        const bool validWSigma = std::isfinite(electron_W_Sigma) && electron_W_Sigma > 0.0f;

        const double W2_truth = validWTruth ? electron_W_truth * electron_W_truth : -1.0;
        const double W2_EM = validWEM ? electron_W_EM * electron_W_EM : -1.0;
        const double W2_DA = validWDA ? electron_W_DA * electron_W_DA : -1.0;
        const double W2_Sigma = validWSigma ? electron_W_Sigma * electron_W_Sigma : -1.0;

        if (validWEM) h_W2_EM->Fill(W2_EM);
        if (validWDA) h_W2_DA->Fill(W2_DA);
        if (hasSigmaMethod && validWSigma) h_W2_Sigma->Fill(W2_Sigma);
        if (validWTruth) h_W2_truth->Fill(W2_truth);

        if (validWTruth && validWEM) g_W2_EM->SetPoint(n_g_W2_EM++, W2_truth, W2_EM);
        if (validWTruth && validWDA) g_W2_DA->SetPoint(n_g_W2_DA++, W2_truth, W2_DA);
        if (hasSigmaMethod && validWTruth && validWSigma) g_W2_Sigma->SetPoint(n_g_W2_Sigma++, W2_truth, W2_Sigma);

        const double m_p_sq = MASS_PROTON * MASS_PROTON;
        const double xpom_denominator_EM    = (validEMQ2 && validWEM)       ? (electron_Q2_EM    + W2_EM    - m_p_sq) : -1.0;
        const double xpom_denominator_truth = (validTruthQ2 && validWTruth) ? (electron_Q2_truth + W2_truth - m_p_sq) : -1.0;
        const double xpom_denominator_DA    = (validDAQ2 && validWDA)       ? (electron_Q2_DA    + W2_DA    - m_p_sq) : -1.0;
        const double xpom_denominator_Sigma = (hasSigmaMethod && validSigmaQ2 && validWSigma) ? (electron_Q2_Sigma + W2_Sigma - m_p_sq): -1.0;

        // Response matrices (Electron method)
        if (validTruthQ2 && validEMQ2) h_Response_Q2->Fill(electron_Q2_truth, electron_Q2_EM);
        if (validTruthX && validEMX) h_Response_x->Fill(electron_x_truth, electron_x_EM);
        if (validTruthY && validEMY) h_Response_y->Fill(electron_y_truth, electron_y_EM);

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

            h_Ep_e->Fill(E);
            h_phi_e->Fill(reco_phi);
            h_pT_e->Fill(reco_pT);
            hasRecoElectron = true;
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

            h_Ep_e_truth->Fill(mc_E);
            h_phi_e_truth->Fill(truth_phi);
            h_pT_e_truth->Fill(truth_pT);
            hasTruthElectron = true;
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
        if (valid_MX2_reco && valid_MX2_truth) h_MX2_corr->Fill(MX2_truth, MX2_reco);
        if (valid_MX2_reco && valid_MX2_truth) g_MX2->SetPoint(n_g_MX2++, MX2_truth, MX2_reco);

        //=================================================================
        // DIFFRACTIVE: Mandelstam t analysis
        //=================================================================
        bool xpom_truth_all_filled = false;
        bool xpom_truth_b0_filled = false;
        bool xpom_truth_rp_filled = false;
        bool xpom_reco_em_b0_filled = false;
        bool xpom_reco_da_b0_filled = false;
        bool xpom_reco_sigma_b0_filled = false;
        bool xpom_reco_em_rp_filled = false;
        bool xpom_reco_da_rp_filled = false;
        bool xpom_reco_sigma_rp_filled = false;
        bool xpom_reco_em_all_filled = false;
        bool xpom_reco_da_all_filled = false;
        bool xpom_reco_sigma_all_filled = false;
        bool has_truth_lead = false;
        double lead_truth_pz = -1.0;
        P3MVector lead_truth(0.0, 0.0, 0.0, 0.0);
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
            if (std::isfinite(t_val)) {
                h_t_MC->Fill(t_val);
                h_dsigma_dt_MC->Fill(t_val);
                if (valid_MX2_truth) {
                    h_MX2_t_truth->Fill(MX2_truth, t_val);
                }
            }
            h_theta_MC->Fill(p.Theta() * 1000.0);

            const double xL_pz = p.Pz() / beams.p_beam.Pz();
            if (std::isfinite(xL_pz)) {
                h_xL_MC->Fill(xL_pz);
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

        if (hasTSProtons && tsassoc_rec_id && tsassoc_sim_id && tsre_px_array && tsre_py_array && tsre_pz_array) {
            for (unsigned int j = 0; j < tsassoc_rec_id->GetSize(); j++) {
                auto mc_idx = (*tsassoc_sim_id)[j];
                if (mc_idx >= (unsigned)mc_pdg_array.GetSize()) continue;
                if (mc_genStatus_array[mc_idx] != 1 || mc_pdg_array[mc_idx] != 2212) continue;

                P3MVector p_reco((*tsre_px_array)[j], (*tsre_py_array)[j], (*tsre_pz_array)[j], mc_mass_array[mc_idx]);
                undoAfterburn(p_reco);

                if (p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;
                h_theta_B0->Fill(p_reco.Theta() * 1000.0);

                P3MVector p_truth(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
                undoAfterburn(p_truth);

                const double t_reco_abs = TMath::Abs(CalcT(beams.p_beam, p_reco));
                const double t_truth_abs = TMath::Abs(CalcT(beams.p_beam, p_truth));
                if (std::isfinite(t_reco_abs)) {
                    h_t_B0->Fill(t_reco_abs);
                    h_dsigma_dt_B0->Fill(t_reco_abs);
                }
                if (std::isfinite(t_truth_abs) && std::isfinite(t_reco_abs)) {
                    h_t_corr_B0->Fill(t_truth_abs, t_reco_abs);
                    g_t_B0->SetPoint(n_g_t_B0++, t_truth_abs, t_reco_abs);
                    if (t_truth_abs > 1e-6) {
                        const double rel = (t_reco_abs - t_truth_abs) / t_truth_abs;
                        h_t_res_B0->Fill(rel);
                        h_t_RelRes_binned_B0->Fill(t_truth_abs, rel);
                    }
                }
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

                if (!xpom_truth_b0_filled && xpom_truth_def > 0.0) {
                    h_xpom_truth_B0->Fill(xpom_truth_def);
                    if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                        h_beta_truth_B0->Fill(beta_truth);
                    }
                    xpom_truth_b0_filled = true;
                }
                if (!xpom_truth_all_filled && xpom_truth_def > 0.0) {
                    h_xpom_truth_all->Fill(xpom_truth_def);
                    if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                        h_beta_truth_all->Fill(beta_truth);
                    }
                    xpom_truth_all_filled = true;
                }

                if (!xpom_reco_em_b0_filled && xpom_denominator_EM > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                    const double xpom_reco_def = (electron_Q2_EM + MX2_reco + t_reco_abs) / xpom_denominator_EM;
                    if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                        h_xpom_reco_EM_B0->Fill(xpom_reco_def);
                        xpom_reco_em_b0_filled = true;
                        if (xpom_truth_def > 0.0) {
                            h_Response_xpom_EM_B0->Fill(xpom_truth_def, xpom_reco_def);
                            g_xpom_EM_B0->SetPoint(n_g_xpom_em_b0++, xpom_truth_def, xpom_reco_def);
                        }
                        if (validEMX) {
                            const double beta_reco = electron_x_EM / xpom_reco_def;
                            if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                h_beta_reco_EM_B0->Fill(beta_reco);
                                if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    g_beta_EM_B0->SetPoint(n_g_beta_em_b0++, beta_truth, beta_reco);
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

                if (!xpom_reco_da_b0_filled && xpom_denominator_DA > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                    const double xpom_reco_def = (electron_Q2_DA + MX2_reco + t_reco_abs) / xpom_denominator_DA;
                    if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                        h_xpom_reco_DA_B0->Fill(xpom_reco_def);
                        xpom_reco_da_b0_filled = true;
                        if (xpom_truth_def > 0.0) {
                            h_Response_xpom_DA_B0->Fill(xpom_truth_def, xpom_reco_def);
                            g_xpom_DA_B0->SetPoint(n_g_xpom_da_b0++, xpom_truth_def, xpom_reco_def);
                        }
                        if (validDAX) {
                            const double beta_reco = electron_x_DA / xpom_reco_def;
                            if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                h_beta_reco_DA_B0->Fill(beta_reco);
                                if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    g_beta_DA_B0->SetPoint(n_g_beta_da_b0++, beta_truth, beta_reco);
                                }
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
                        if (xpom_truth_def > 0.0) {
                            h_Response_xpom_Sigma_B0->Fill(xpom_truth_def, xpom_reco_def);
                            g_xpom_Sigma_B0->SetPoint(n_g_xpom_sigma_b0++, xpom_truth_def, xpom_reco_def);
                        }
                        if (hasSigmaMethod && validSigmaX) {
                            const double beta_reco = electron_x_Sigma / xpom_reco_def;
                            if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                h_beta_reco_Sigma_B0->Fill(beta_reco);
                                if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                    g_beta_Sigma_B0->SetPoint(n_g_beta_sigma_b0++, beta_truth, beta_reco);
                                }
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
                    g_xL_B0->SetPoint(n_g_xL_B0++, xL_truth_pz, xL_reco_pz);
                }
            }
        }

        if (hasRPProtons && rp_px_array && rp_py_array && rp_pz_array && rp_mass_array && rp_pdg_array) {
            for (int j = 0; j < rp_px_array->GetSize(); j++) {
                if ((*rp_pdg_array)[j] != 2212) continue;

                P3MVector p_rp((*rp_px_array)[j], (*rp_py_array)[j], (*rp_pz_array)[j], (*rp_mass_array)[j]);
                const double theta_mrad = p_rp.Theta() * 1000.0;
                if (theta_mrad <= 5.0) {
                    h_theta_RP->Fill(theta_mrad);
                }

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
                    }
                    if (std::isfinite(t_truth_abs) && std::isfinite(t_reco_abs)) {
                        h_t_corr_RP->Fill(t_truth_abs, t_reco_abs);
                        g_t_RP->SetPoint(n_g_t_RP++, t_truth_abs, t_reco_abs);
                        if (t_truth_abs > 1e-6) {
                            const double rel = (t_reco_abs - t_truth_abs) / t_truth_abs;
                            h_t_res_RP->Fill(rel);
                            h_t_RelRes_binned_RP->Fill(t_truth_abs, rel);
                        }
                    }
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

                    if (!xpom_truth_rp_filled && xpom_truth_def > 0.0) {
                        h_xpom_truth_RP->Fill(xpom_truth_def);
                        if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                            h_beta_truth_RP->Fill(beta_truth);
                        }
                        xpom_truth_rp_filled = true;
                    }
                    if (!xpom_truth_all_filled && xpom_truth_def > 0.0) {
                        h_xpom_truth_all->Fill(xpom_truth_def);
                        if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                            h_beta_truth_all->Fill(beta_truth);
                        }
                        xpom_truth_all_filled = true;
                    }

                    if (!xpom_reco_em_rp_filled && xpom_denominator_EM > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                        const double xpom_reco_def = (electron_Q2_EM + MX2_reco + t_reco_abs) / xpom_denominator_EM;
                        if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                            h_xpom_reco_EM_RP->Fill(xpom_reco_def);
                            xpom_reco_em_rp_filled = true;
                            if (xpom_truth_def > 0.0) {
                                h_Response_xpom_EM_RP->Fill(xpom_truth_def, xpom_reco_def);
                                g_xpom_EM_RP->SetPoint(n_g_xpom_em_rp++, xpom_truth_def, xpom_reco_def);
                            }
                            if (validEMX) {
                                const double beta_reco = electron_x_EM / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_EM_RP->Fill(beta_reco);
                                    if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                        g_beta_EM_RP->SetPoint(n_g_beta_em_rp++, beta_truth, beta_reco);
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

                    if (!xpom_reco_da_rp_filled && xpom_denominator_DA > 0.0 && valid_MX2_reco && std::isfinite(t_reco_abs)) {
                        const double xpom_reco_def = (electron_Q2_DA + MX2_reco + t_reco_abs) / xpom_denominator_DA;
                        if (std::isfinite(xpom_reco_def) && xpom_reco_def > 0.0) {
                            h_xpom_reco_DA_RP->Fill(xpom_reco_def);
                            xpom_reco_da_rp_filled = true;
                            if (xpom_truth_def > 0.0) {
                                h_Response_xpom_DA_RP->Fill(xpom_truth_def, xpom_reco_def);
                                g_xpom_DA_RP->SetPoint(n_g_xpom_da_rp++, xpom_truth_def, xpom_reco_def);
                            }
                            if (validDAX) {
                                const double beta_reco = electron_x_DA / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_DA_RP->Fill(beta_reco);
                                    if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                        g_beta_DA_RP->SetPoint(n_g_beta_da_rp++, beta_truth, beta_reco);
                                    }
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
                            if (xpom_truth_def > 0.0) {
                                h_Response_xpom_Sigma_RP->Fill(xpom_truth_def, xpom_reco_def);
                                g_xpom_Sigma_RP->SetPoint(n_g_xpom_sigma_rp++, xpom_truth_def, xpom_reco_def);
                            }
                            if (hasSigmaMethod && validSigmaX) {
                                const double beta_reco = electron_x_Sigma / xpom_reco_def;
                                if (beta_reco > 0.0 && beta_reco <= 1.0) {
                                    h_beta_reco_Sigma_RP->Fill(beta_reco);
                                    if (beta_truth_valid && beta_truth > 0.0 && beta_truth <= 1.0) {
                                        g_beta_Sigma_RP->SetPoint(n_g_beta_sigma_rp++, beta_truth, beta_reco);
                                    }
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
                        g_xL_RP->SetPoint(n_g_xL_RP++, xL_truth_pz, xL_reco_pz);
                    }
                }
            }
        }

    }

    // Calculate and scale differential cross sections d(sigma)/dt
    const double sigma_total = 4.22; // nb for 10x100 GeV configuration
    const double N_gen = 100000.0;   // Number of generated events
    const double scale_factor = sigma_total / N_gen;

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

    TH2D* h_Response_Q2_rowNorm = (TH2D*)h_Response_Q2->Clone("Response_Q2_rowNorm");
    TH2D* h_Response_Q2_colNorm = (TH2D*)h_Response_Q2->Clone("Response_Q2_colNorm");
    TH2D* h_Response_x_rowNorm  = (TH2D*)h_Response_x->Clone("Response_x_rowNorm");
    TH2D* h_Response_x_colNorm  = (TH2D*)h_Response_x->Clone("Response_x_colNorm");
    TH2D* h_Response_y_rowNorm  = (TH2D*)h_Response_y->Clone("Response_y_rowNorm");
    TH2D* h_Response_y_colNorm  = (TH2D*)h_Response_y->Clone("Response_y_colNorm");

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

    // Write all histograms and TTree to the output file
    outputFile->cd();
    TParameter<Long64_t> nEventsParam("nEventsProcessed", nProcessed);
    nEventsParam.Write();
    tree->Write();

    // Write Q2/xy histograms
    h_RelRes_Q2_EM->Write();
    h_RelRes_Q2_DA->Write();
    h_RelRes_Q2_Sigma->Write();
    h_RelRes_Q2_binned_EM->Write();
    h_RelRes_Q2_binned_DA->Write();
    h_RelRes_Q2_binned_Sigma->Write();
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
    h_xQ2_truth->Write();

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

    // Write NEW inclusive DIS histograms
    h_W2_EM->Write(); h_W2_DA->Write(); h_W2_Sigma->Write(); h_W2_truth->Write();
    h_Ep_e->Write(); h_phi_e->Write(); h_pT_e->Write();
    h_Ep_e_truth->Write(); h_phi_e_truth->Write(); h_pT_e_truth->Write();
    h_Response_Q2->Write(); h_Response_x->Write(); h_Response_y->Write();
    h_Response_Q2_rowNorm->Write(); h_Response_Q2_colNorm->Write();
    h_Response_x_rowNorm->Write(); h_Response_x_colNorm->Write();
    h_Response_y_rowNorm->Write(); h_Response_y_colNorm->Write();

    g_Q2_EM->Write(); g_Q2_DA->Write(); g_Q2_Sigma->Write();
    g_x_EM->Write(); g_x_DA->Write(); g_x_Sigma->Write();
    g_y_EM->Write(); g_y_DA->Write(); g_y_Sigma->Write();
    g_W2_EM->Write(); g_W2_DA->Write(); g_W2_Sigma->Write();
    g_Ep_e->Write(); g_phi_e->Write(); g_pT_e->Write();

    // Write diffractive t histograms/graphs
    h_t_MC->Write();
    h_t_B0->Write();
    h_t_RP_histo->Write();
    h_dsigma_dt_MC->Write();
    h_dsigma_dt_B0->Write();
    h_dsigma_dt_RP->Write();
    h_theta_MC->Write();
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
    h_Response_xpom_EM_B0->Write();
    h_Response_xpom_EM_RP->Write();
    h_Response_xpom_DA_B0->Write();
    h_Response_xpom_DA_RP->Write();
    h_Response_xpom_Sigma_B0->Write();
    h_Response_xpom_Sigma_RP->Write();
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
    h_beta_Q2_truth->Write();
    h_t_Q2_truth->Write();
    h_xpom_Q2_truth->Write();
    h_beta_t_truth->Write();
    h_xbj_t_truth->Write();
    h_xpom_t_truth->Write();
    h_xpom_beta_truth->Write();
    h_xbj_beta_truth->Write();
    h_xbj_xpom_truth->Write();
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
    h_MX2_t_truth->Write();
    h_MX2_t_B0->Write();
    h_MX2_t_RP->Write();

    outputFile->Close();
    delete events;
    delete outputFile;

    std::cout << "\nInclusive kinematics summary:" << std::endl;
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
