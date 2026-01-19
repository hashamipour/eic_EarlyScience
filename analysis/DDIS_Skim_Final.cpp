// Unified DDIS skimmer: Q2/xy + t analysis + x_pom/beta comparisons
// g++ DDIS_Skim_Final.cpp -o DDIS_Skim_Final $(root-config --cflags --glibs)
// ./DDIS_Skim_Final filelist.txt [output.root]

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <TGraph.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>

#include "Utility.hpp"
#include "RecoMethods.hpp"

// ROOT Math objects
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"

using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

// Global afterburner correction parameters
Float_t fXAngle{-0.025};
RotationX rotAboutX;
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

static void undoAfterburnAndCalc(P3MVector& p, P3MVector& k) {
    P3MVector p_beam(fXAngle * p.E(), 0., p.E(), p.M());
    P3MVector e_beam(0., 0., -k.E(), k.M());

    P3MVector CoM_boost = p_beam + e_beam;
    vBoostToCoM.SetXYZ(-CoM_boost.X() / CoM_boost.E(),
                       -CoM_boost.Y() / CoM_boost.E(),
                       -CoM_boost.Z() / CoM_boost.E());

    p_beam = boost(p_beam, vBoostToCoM);
    e_beam = boost(e_beam, vBoostToCoM);

    Float_t fRotY = -1.0 * TMath::ATan2(p_beam.X(), p_beam.Z());
    Float_t fRotX =  1.0 * TMath::ATan2(p_beam.Y(), p_beam.Z());

    rotAboutY = RotationY(fRotY);
    rotAboutX = RotationX(fRotX);

    p_beam = rotAboutY(p_beam);
    p_beam = rotAboutX(p_beam);
    e_beam = rotAboutY(e_beam);
    e_beam = rotAboutX(e_beam);

    P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
    vBoostToHoF.SetXYZ(HoF_boost.X() / HoF_boost.E(),
                       HoF_boost.Y() / HoF_boost.E(),
                       HoF_boost.Z() / HoF_boost.E());

    p_beam = boost(p_beam, vBoostToHoF);
    e_beam = boost(e_beam, vBoostToHoF);

    p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
    k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}

static void undoAfterburn(P3MVector& a) {
    a = boost(a, vBoostToCoM);
    a = rotAboutY(a);
    a = rotAboutX(a);
    a = boost(a, vBoostToHoF);
}

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
                mc_idx = assoc_sim_id[j];
                break;
            }
        }
        if (mc_idx < 0) continue;

        double E_reco = re_energy[i];
        double pz_reco = re_pz[i];
        sumEPz_reco += (E_reco - pz_reco);

        double px_mc = mc_px[mc_idx];
        double py_mc = mc_py[mc_idx];
        double pz_mc = mc_pz[mc_idx];
        double m_mc = mc_mass[mc_idx];
        double E_mc = TMath::Sqrt(px_mc * px_mc + py_mc * py_mc + pz_mc * pz_mc + m_mc * m_mc);
        sumEPz_truth += (E_mc - pz_mc);
    }
}

static double XpomFromDefinition(double MX2, double Q2, double t_abs, double W) {
    double m_p_sq = fMass_proton * fMass_proton;
    double W2 = W * W;
    double denominator = W2 + Q2 - m_p_sq;
    if (denominator <= 0.0) return -999.0;
    return (MX2 + Q2 + t_abs) / denominator;
}

static bool ValidUnitInterval(double val) {
    return (val > 0.0 && val <= 1.0);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fileList.txt> [output.root]" << std::endl;
        return 1;
    }

    TString fileList = argv[1];
    TString outputName = (argc > 2) ? argv[2] : "DDIS_Final_output.root";

    std::ifstream fileListStream(fileList.Data());
    if (!fileListStream.is_open()) {
        std::cerr << "Error: Could not open file list " << fileList << std::endl;
        return 1;
    }

    TChain* events = new TChain("events");
    std::string fileName;
    int nFiles = 0;
    while (getline(fileListStream, fileName)) {
        TString tmp = fileName;
        if (!std::filesystem::exists(tmp.Data())) {
            std::cerr << "Warning: File does not exist: " << fileName << std::endl;
            continue;
        }
        auto inputRootFile = TFile::Open(tmp);
        events->Add((TString)fileName);
        inputRootFile->Close();
        nFiles++;
    }
    std::cout << "No. of files: " << nFiles << "; no. of events: " << events->GetEntries() << std::endl;

    TFile* outputFile = new TFile(outputName, "RECREATE");

    TTree* tree = new TTree("Q2_tree", "Q2 Kinematics Data");
    float out_x_truth, out_x_EM, out_x_DA, out_x_Sigma;
    float out_y_truth, out_y_EM, out_y_DA, out_y_Sigma;
    tree->Branch("x_truth",    &out_x_truth,   "x_truth/F");
    tree->Branch("x_EM",       &out_x_EM,      "x_EM/F");
    tree->Branch("x_DA",       &out_x_DA,      "x_DA/F");
    tree->Branch("x_Sigma",    &out_x_Sigma,   "x_Sigma/F");
    tree->Branch("y_truth",    &out_y_truth,   "y_truth/F");
    tree->Branch("y_EM",       &out_y_EM,      "y_EM/F");
    tree->Branch("y_DA",       &out_y_DA,      "y_DA/F");
    tree->Branch("y_Sigma",    &out_y_Sigma,   "y_Sigma/F");

    // Q2 / x / y binning
    int n_bins = 10;
    std::vector<Double_t> bin_edges_Q2 = GetRoundedLogBins(3.4, 150.0, n_bins);
    n_bins = bin_edges_Q2.size() - 1;
    std::vector<Double_t> x_bins = GetLogBins(1.0e-4, 1.0, 25);
    const int n_y_bins = 20;

    // Q2 histograms
    TH1D* h_Q2_truth = new TH1D("h_Q2_truth", "Q^2;# of events", n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_EM = new TH1D("h_Q2_EM", ";Q^{2}", n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_DA = new TH1D("h_Q2_DA", ";Q^{2}", n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_Sigma = new TH1D("h_Q2_Sigma", ";Q^{2}", n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_EICRecon = new TH1D("Q2_EICRecon", "EICRecon Q^{2};Q^{2};Counts", n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_calc = new TH1D("Q2_calc", "Q^{2} from s*x*y;Q^{2};Counts", n_bins, bin_edges_Q2.data());

    TH1D* h_RelRes_Q2_EM = new TH1D("Q2_RelRes_EM", "electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_Q2_DA = new TH1D("Q2_RelRes_DA", "DA method;#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_Q2_Sigma = new TH1D("Q2_RelRes_Sigma", "Sigma method;#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}", 101, -0.15, 0.15);

    TH2D* h_RelRes_Q2_binned_EM = new TH2D("Q2_RelRes_binned_EM", ";Q^{2} [GeV^{2}];#frac{Q^{2}(EM)-Q^{2}(MC)}{Q^{2}(MC)}",
                                          n_bins, bin_edges_Q2.data(), 51, -0.15, 0.15);
    TH2D* h_RelRes_Q2_binned_DA = new TH2D("Q2_RelRes_binned_DA", ";Q^{2} [GeV^{2}];#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",
                                          n_bins, bin_edges_Q2.data(), 51, -0.15, 0.15);
    TH2D* h_RelRes_Q2_binned_Sigma = new TH2D("Q2_RelRes_binned_Sigma", ";Q^{2} [GeV^{2}];#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",
                                             n_bins, bin_edges_Q2.data(), 51, -0.15, 0.15);

    TH2D* h_Corr_Q2_EM = new TH2D("Corr_Q2_EM", ";Q^{2}_{MC};Q^{2}_{EM}", n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_DA = new TH2D("Corr_Q2_DA", ";Q^{2}_{MC};Q^{2}_{DA}", n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_Sigma = new TH2D("Corr_Q2_Sigma", ";Q^{2}_{MC};Q^{2}_{Sigma}", n_bins, bin_edges_Q2.data(), n_bins, bin_edges_Q2.data());

    // x histograms
    TH1D* h_x_EM = new TH1D("x_EM", "electron method;x_{Bj}", x_bins.size() - 1, x_bins.data());
    TH1D* h_x_DA = new TH1D("x_DA", "DA method;x_{Bj}", x_bins.size() - 1, x_bins.data());
    TH1D* h_x_Sigma = new TH1D("x_Sigma", "Sigma method;x_{Bj}", x_bins.size() - 1, x_bins.data());
    TH1D* h_x_truth = new TH1D("x_truth", "truth;x_{Bj}", x_bins.size() - 1, x_bins.data());

    TH1D* h_RelRes_x_EM = new TH1D("x_RelRes_EM", "electron method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_x_DA = new TH1D("x_RelRes_DA", "DA method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_x_Sigma = new TH1D("x_RelRes_Sigma", "Sigma method;#frac{x(Reco)-x(MC)}{x(MC)}", 101, -0.15, 0.15);

    TH2D* h_RelRes_x_binned_EM = new TH2D("x_RelRes_binned_EM", "x: truth vs rel. res (EM);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",
                                          x_bins.size() - 1, x_bins.data(), 51, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_DA = new TH2D("x_RelRes_binned_DA", "x: truth vs rel. res (DA);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",
                                          x_bins.size() - 1, x_bins.data(), 51, -0.15, 0.15);
    TH2D* h_RelRes_x_binned_Sigma = new TH2D("x_RelRes_binned_Sigma", "x: truth vs rel. res (Sigma);x_{truth};#frac{x(Reco)-x(MC)}{x(MC)}",
                                             x_bins.size() - 1, x_bins.data(), 51, -0.15, 0.15);

    TH2D* h_Corr_x_EM = new TH2D("x_Corr_EM", "x correlation (EM);x_{truth};x_{EM}",
                                 x_bins.size() - 1, x_bins.data(), x_bins.size() - 1, x_bins.data());
    TH2D* h_Corr_x_DA = new TH2D("x_Corr_DA", "x correlation (DA);x_{truth};x_{DA}",
                                 x_bins.size() - 1, x_bins.data(), x_bins.size() - 1, x_bins.data());
    TH2D* h_Corr_x_Sigma = new TH2D("x_Corr_Sigma", "x correlation (Sigma);x_{truth};x_{Sigma}",
                                    x_bins.size() - 1, x_bins.data(), x_bins.size() - 1, x_bins.data());

    // y histograms
    TH1D* h_y_EM = new TH1D("y_EM", "electron method;y", n_y_bins, 0.0, 1.0);
    TH1D* h_y_DA = new TH1D("y_DA", "DA method;y", n_y_bins, 0.0, 1.0);
    TH1D* h_y_Sigma = new TH1D("y_Sigma", "Sigma method;y", n_y_bins, 0.0, 1.0);
    TH1D* h_y_truth = new TH1D("y_truth", "truth;y", n_y_bins, 0.0, 1.0);

    TH1D* h_RelRes_y_EM = new TH1D("y_RelRes_EM", "electron method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_y_DA = new TH1D("y_RelRes_DA", "DA method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);
    TH1D* h_RelRes_y_Sigma = new TH1D("y_RelRes_Sigma", "Sigma method;#frac{y(Reco)-y(MC)}{y(MC)}", 101, -0.15, 0.15);

    TH2D* h_RelRes_y_binned_EM = new TH2D("y_RelRes_binned_EM", "y: truth vs rel. res (EM);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",
                                          n_y_bins, 0.0, 1.0, 51, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_DA = new TH2D("y_RelRes_binned_DA", "y: truth vs rel. res (DA);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",
                                          n_y_bins, 0.0, 1.0, 51, -0.15, 0.15);
    TH2D* h_RelRes_y_binned_Sigma = new TH2D("y_RelRes_binned_Sigma", "y: truth vs rel. res (Sigma);y_{truth};#frac{y(Reco)-y(MC)}{y(MC)}",
                                             n_y_bins, 0.0, 1.0, 51, -0.15, 0.15);

    TH2D* h_Corr_y_EM = new TH2D("y_Corr_EM", "y correlation (EM);y_{truth};y_{EM}", n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_DA = new TH2D("y_Corr_DA", "y correlation (DA);y_{truth};y_{DA}", n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_Sigma = new TH2D("y_Corr_Sigma", "y correlation (Sigma);y_{truth};y_{Sigma}", n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);

    // Q2/x/y resolution profiles for circle maps
    TProfile2D* h_Q2_RelRes_vs_xy_EM = new TProfile2D("Q2_RelRes_vs_xy_EM",
        "Q^{2} Rel. Res. binned in (x_{Bj}, y) (EM);x_{Bj};y",
        x_bins.size() - 1, x_bins.data(), n_y_bins, 0.0, 1.0, "s");
    TProfile2D* h_Q2_RelRes_vs_xy_DA = new TProfile2D("Q2_RelRes_vs_xy_DA",
        "Q^{2} Rel. Res. binned in (x_{Bj}, y) (DA);x_{Bj};y",
        x_bins.size() - 1, x_bins.data(), n_y_bins, 0.0, 1.0, "s");
    TProfile2D* h_Q2_RelRes_vs_xy_Sigma = new TProfile2D("Q2_RelRes_vs_xy_Sigma",
        "Q^{2} Rel. Res. binned in (x_{Bj}, y) (Sigma);x_{Bj};y",
        x_bins.size() - 1, x_bins.data(), n_y_bins, 0.0, 1.0, "s");

    TProfile2D* h_x_RelRes_vs_xQ2_EM = new TProfile2D("x_RelRes_vs_xQ2_EM",
        "x Rel. Res. binned in (x_{Bj}, Q^{2}) (EM);x_{Bj};Q^{2} [GeV^{2}]",
        x_bins.size() - 1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_x_RelRes_vs_xQ2_DA = new TProfile2D("x_RelRes_vs_xQ2_DA",
        "x Rel. Res. binned in (x_{Bj}, Q^{2}) (DA);x_{Bj};Q^{2} [GeV^{2}]",
        x_bins.size() - 1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_x_RelRes_vs_xQ2_Sigma = new TProfile2D("x_RelRes_vs_xQ2_Sigma",
        "x Rel. Res. binned in (x_{Bj}, Q^{2}) (Sigma);x_{Bj};Q^{2} [GeV^{2}]",
        x_bins.size() - 1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");

    TProfile2D* h_y_RelRes_vs_xQ2_EM = new TProfile2D("y_RelRes_vs_xQ2_EM",
        "y Rel. Res. binned in (x_{Bj}, Q^{2}) (EM);x_{Bj};Q^{2} [GeV^{2}]",
        x_bins.size() - 1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_y_RelRes_vs_xQ2_DA = new TProfile2D("y_RelRes_vs_xQ2_DA",
        "y Rel. Res. binned in (x_{Bj}, Q^{2}) (DA);x_{Bj};Q^{2} [GeV^{2}]",
        x_bins.size() - 1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_y_RelRes_vs_xQ2_Sigma = new TProfile2D("y_RelRes_vs_xQ2_Sigma",
        "y Rel. Res. binned in (x_{Bj}, Q^{2}) (Sigma);x_{Bj};Q^{2} [GeV^{2}]",
        x_bins.size() - 1, x_bins.data(), n_bins, bin_edges_Q2.data(), "s");

    // Hadronic final state variables
    TH1D* h_EPz_truth = new TH1D("h_EPz_truth", "MC Truth Sum(E-p_{z}) - Matched Particles Only;#Sigma(E-p_{z}) [GeV];Counts", 50, 0, 25);
    TH1D* h_EPz = new TH1D("h_EPz", "Reco Sum(E-p_{z}) - Matched Particles Only;#Sigma(E-p_{z}) [GeV];Counts", 50, 0, 25);
    TH2D* h_EPz_2D = new TH2D("h_EPz_2D", "E-p_{z} Truth vs Reco;#Sigma(E-p_{z})_{truth} [GeV];#Sigma(E-p_{z})_{reco} [GeV]",
                              50, 0, 25, 50, 0, 25);

    TH1D* h_eta_max = new TH1D("h_eta_max", "Maximum Pseudorapidity per Event (Reco);#eta_{max};Counts", 50, -4.0, 6.0);
    TH1D* h_eta_max_truth = new TH1D("h_eta_max_truth", "Maximum Pseudorapidity per Event (Truth);#eta_{max};Counts", 50, -4.0, 6.0);
    TH1D* h_MX2 = new TH1D("h_MX2", "Hadronic Invariant Mass Squared (Reco);M_{X}^{2} [GeV^{2}];Counts", 100, 0.0, 200.0);
    TH1D* h_MX2_truth = new TH1D("h_MX2_truth", "Hadronic Invariant Mass Squared (Truth);M_{X}^{2} [GeV^{2}];Counts", 100, 0.0, 200.0);
    TH1D* h_MX2_MC = new TH1D("MX2_MC", "Truth M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts", 100, 0.0, 200.0);
    TH1D* h_MX2_eX = new TH1D("MX2_eX", "Reco M_{X}^{2};M_{X}^{2} [GeV^{2}];Counts", 100, 0.0, 200.0);

    TH2D* h_MX2_corr = new TH2D("MX2_corr", "M_{X}^{2} Correlation;Truth M_{X}^{2} [GeV^{2}];Reco M_{X}^{2} [GeV^{2}]",
                                100, 0.0, 200.0, 100, 0.0, 200.0);
    TH2D* h_MX2_RelRes_binned = new TH2D("MX2_RelRes_binned",
        "M_{X}^{2} resolution;M_{X,truth}^{2} [GeV^{2}];(M_{X,reco}^{2}-M_{X,truth}^{2})/M_{X,truth}^{2}",
        100, 0.0, 200.0, 101, -0.5, 0.5);
    TProfile2D* h_MX2_RelRes_vs_MX2Q2 = new TProfile2D("MX2_RelRes_vs_MX2Q2",
        "M_{X}^{2} resolution vs (M_{X}^{2}, Q^{2});M_{X}^{2} [GeV^{2}];Q^{2} [GeV^{2}]",
        50, 0.0, 200.0, n_bins, bin_edges_Q2.data(), "s");

    // t-analysis binning
    std::vector<Double_t> t_bins_low = GetLogBins(1e-3, 0.5, 20);
    std::vector<Double_t> t_bins_high = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.6};
    std::vector<Double_t> t_bins = t_bins_low;
    for (size_t i = 1; i < t_bins_high.size(); i++) {
        t_bins.push_back(t_bins_high[i]);
    }

    TH1D* h_t_MC = new TH1D("t_MC", "Truth Mandelstam t;|t| [GeV^{2}];Counts", t_bins.size() - 1, t_bins.data());
    TH1D* h_t_B0 = new TH1D("t_B0", "B0 Reco Mandelstam t;|t| [GeV^{2}];Counts", t_bins.size() - 1, t_bins.data());
    TH1D* h_t_RP_histo = new TH1D("t_RP_histo", "RP Reco Mandelstam t;|t| [GeV^{2}];Counts", t_bins.size() - 1, t_bins.data());

    TH1D* h_t_B0_cutFirstBin = new TH1D("t_B0_cutFirstBin", "B0 Reco Mandelstam t (first bin cut);|t| [GeV^{2}];Counts",
                                        t_bins.size() - 1, t_bins.data());

    TH1D* h_dsigma_dt_MC = new TH1D("dsigma_dt_MC", "Truth d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                    t_bins.size() - 1, t_bins.data());
    TH1D* h_dsigma_dt_B0 = new TH1D("dsigma_dt_B0", "B0 Reco d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                    t_bins.size() - 1, t_bins.data());
    TH1D* h_dsigma_dt_RP = new TH1D("dsigma_dt_RP", "RP Reco d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                    t_bins.size() - 1, t_bins.data());

    TH1D* h_t_res_B0 = new TH1D("t_res_B0", "B0 t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 100, -2.0, 2.0);
    TH1D* h_t_res_RP = new TH1D("t_res_RP", "RP t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 100, -2.0, 2.0);

    TH2D* h_t_corr_B0 = new TH2D("t_corr_B0", "B0 |t| Correlation;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                 t_bins.size() - 1, t_bins.data(), t_bins.size() - 1, t_bins.data());
    TH2D* h_t_corr_RP = new TH2D("t_corr_RP", "RP |t| Correlation;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                 t_bins.size() - 1, t_bins.data(), t_bins.size() - 1, t_bins.data());
    TH2D* h_t_corr_B0_cutFirstBin = new TH2D("t_corr_B0_cutFirstBin", "B0 |t| Correlation (first bin cut);Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                            t_bins.size() - 1, t_bins.data(), t_bins.size() - 1, t_bins.data());

    TH2D* h_t_RelRes_binned_B0 = new TH2D("t_RelRes_binned_B0", "t relres vs truth (B0);|t|_{truth};(|t|_{reco}-|t|_{truth})/|t|_{truth}",
                                          t_bins.size() - 1, t_bins.data(), 101, -0.5, 0.5);
    TH2D* h_t_RelRes_binned_RP = new TH2D("t_RelRes_binned_RP", "t relres vs truth (RP);|t|_{truth};(|t|_{reco}-|t|_{truth})/|t|_{truth}",
                                          t_bins.size() - 1, t_bins.data(), 101, -0.5, 0.5);

    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_all_TS = new TH1D("theta_all_TS", "All Truth-Seeded Proton Angles;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_RP = new TH1D("theta_RP", "RP Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);

    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_B0 = new TH1D("xL_B0", "B0 Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_RP = new TH1D("xL_RP", "RP Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_B0_cutFirstBin = new TH1D("xL_B0_cutFirstBin", "B0 Reco x_{L} (first bin cut);x_{L};Counts", 30, 0.75, 1.05);

    TH1D* h_xL_res_B0 = new TH1D("xL_res_B0", "B0 x_L Resolution;(x_{L,reco}-x_{L,truth})/x_{L,truth};Counts", 100, -2.0, 2.0);
    TH1D* h_xL_res_RP = new TH1D("xL_res_RP", "RP x_L Resolution;(x_{L,reco}-x_{L,truth})/x_{L,truth};Counts", 100, -2.0, 2.0);

    TH2D* h_xL_corr_B0 = new TH2D("xL_corr_B0", "B0 x_L Correlation;Truth x_L;Reco x_L", 100, 0.0, 2.0, 100, 0.0, 2.0);
    TH2D* h_xL_corr_RP = new TH2D("xL_corr_RP", "RP x_L Correlation;Truth x_L;Reco x_L", 100, 0.0, 0.5, 100, 0.0, 0.5);
    TH2D* h_xL_corr_B0_cutFirstBin = new TH2D("xL_corr_B0_cutFirstBin",
                                              "B0 x_L Correlation (first bin cut);Truth x_L;Reco x_L",
                                              100, 0.0, 2.0, 100, 0.0, 2.0);

    TH2D* h_xL_RelRes_binned_B0 = new TH2D("xL_RelRes_binned_B0", "xL relres vs truth (B0);x_{L,truth};(x_{L,reco}-x_{L,truth})/x_{L,truth}",
                                           30, 0.75, 1.05, 101, -0.5, 0.5);
    TH2D* h_xL_RelRes_binned_RP = new TH2D("xL_RelRes_binned_RP", "xL relres vs truth (RP);x_{L,truth};(x_{L,reco}-x_{L,truth})/x_{L,truth}",
                                           30, 0.75, 1.05, 101, -0.5, 0.5);

    const int n_xpom_bins = 20;
    double xpom_bins[n_xpom_bins + 1];
    double xpom_min = 1e-4;
    double xpom_max = 0.4;
    double log_min = TMath::Log10(xpom_min);
    double log_max = TMath::Log10(xpom_max);
    for (int i = 0; i <= n_xpom_bins; i++) {
        xpom_bins[i] = TMath::Power(10, log_min + i * (log_max - log_min) / n_xpom_bins);
    }

    TH1D* h_xpom_MC = new TH1D("xpom_MC", "Truth x_{pom} (from x_L);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_B0 = new TH1D("xpom_B0", "B0 Reco x_{pom} (from x_L);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_RP = new TH1D("xpom_RP", "RP Reco x_{pom} (from x_L);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_B0_cutFirstBin = new TH1D("xpom_B0_cutFirstBin", "B0 Reco x_{pom} (first bin cut);x_{pom};Counts", n_xpom_bins, xpom_bins);

    TH1D* h_xpom_def_MC = new TH1D("xpom_def_MC", "Truth x_{pom} (from definition);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_def_B0 = new TH1D("xpom_def_B0", "B0 Reco x_{pom} (from definition);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_def_RP = new TH1D("xpom_def_RP", "RP Reco x_{pom} (from definition);x_{pom};Counts", n_xpom_bins, xpom_bins);

    TH2D* h_xpom_comp_MC = new TH2D("xpom_comp_MC", "Truth: x_{pom} from x_L vs from definition;x_{pom} (1-x_L);x_{pom} (definition)",
                                    n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_comp_B0 = new TH2D("xpom_comp_B0", "B0 Reco: x_{pom} from x_L vs from definition;x_{pom} (1-x_L);x_{pom} (definition)",
                                    n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_comp_RP = new TH2D("xpom_comp_RP", "RP Reco: x_{pom} from x_L vs from definition;x_{pom} (1-x_L);x_{pom} (definition)",
                                    n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);

    TH1D* h_xpom_res_B0 = new TH1D("xpom_res_B0", "B0 x_pom Resolution (from x_L);(x_{pom,reco}-x_{pom,truth})/x_{pom,truth};Counts", 100, -2.0, 2.0);
    TH1D* h_xpom_res_RP = new TH1D("xpom_res_RP", "RP x_pom Resolution (from x_L);(x_{pom,reco}-x_{pom,truth})/x_{pom,truth};Counts", 100, -2.0, 2.0);

    TH2D* h_xpom_corr_B0 = new TH2D("xpom_corr_B0", "B0 x_pom Correlation (from x_L);Truth x_{pom};Reco x_{pom}",
                                    n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_corr_RP = new TH2D("xpom_corr_RP", "RP x_pom Correlation (from x_L);Truth x_{pom};Reco x_{pom}",
                                    n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_corr_B0_cutFirstBin = new TH2D("xpom_corr_B0_cutFirstBin",
                                                "B0 x_pom Correlation (first bin cut);Truth x_{pom};Reco x_{pom}",
                                                n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);

    TH2D* h_xpom_RelRes_binned_B0 = new TH2D("xpom_RelRes_binned_B0", "xpom relres vs truth (B0);x_{pom,truth};(x_{pom,reco}-x_{pom,truth})/x_{pom,truth}",
                                             n_xpom_bins, xpom_bins, 101, -0.5, 0.5);
    TH2D* h_xpom_RelRes_binned_RP = new TH2D("xpom_RelRes_binned_RP", "xpom relres vs truth (RP);x_{pom,truth};(x_{pom,reco}-x_{pom,truth})/x_{pom,truth}",
                                             n_xpom_bins, xpom_bins, 101, -0.5, 0.5);

    // Beta histograms (from x_pom definition)
    TH1D* h_beta_MC = new TH1D("beta_MC", "Truth #beta (from x_{pom} def);#beta;Counts", 10, 0.0, 1.0);
    TH1D* h_beta_B0 = new TH1D("beta_B0", "B0 Reco #beta (from x_{pom} def);#beta;Counts", 10, 0.0, 1.0);
    TH1D* h_beta_RP = new TH1D("beta_RP", "RP Reco #beta (from x_{pom} def);#beta;Counts", 10, 0.0, 1.0);

    TH1D* h_beta_res_B0 = new TH1D("beta_res_B0", "B0 #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts", 100, -2.0, 2.0);
    TH1D* h_beta_res_RP = new TH1D("beta_res_RP", "RP #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts", 100, -2.0, 2.0);

    TH2D* h_beta_corr_B0 = new TH2D("beta_corr_B0", "B0 #beta Correlation;Truth #beta;Reco #beta", 10, 0.0, 1.0, 10, 0.0, 1.0);
    TH2D* h_beta_corr_RP = new TH2D("beta_corr_RP", "RP #beta Correlation;Truth #beta;Reco #beta", 10, 0.0, 1.0, 10, 0.0, 1.0);

    TH2D* h_beta_RelRes_binned_B0 = new TH2D("beta_RelRes_binned_B0", "beta relres vs truth (B0);#beta_{truth};(#beta_{reco}-#beta_{truth})/#beta_{truth}",
                                             10, 0.0, 1.0, 101, -0.5, 0.5);
    TH2D* h_beta_RelRes_binned_RP = new TH2D("beta_RelRes_binned_RP", "beta relres vs truth (RP);#beta_{truth};(#beta_{reco}-#beta_{truth})/#beta_{truth}",
                                             10, 0.0, 1.0, 101, -0.5, 0.5);

    TH2D* h_beta_vs_Q2 = new TH2D("beta_vs_Q2", "#beta vs Q^{2};Q^{2} [GeV^{2}];#beta", 50, 0.1, 100.0, 50, 0.0, 1.0);
    TH2D* h_beta_vs_xpom = new TH2D("beta_vs_xpom", "#beta vs x_{pom};x_{pom};#beta", 50, 1e-4, 0.4, 50, 0.0, 1.0);
    TH2D* h_beta_vs_t = new TH2D("beta_vs_t", "#beta vs |t|;|t| [GeV^{2}];#beta", 50, 0.0, 2.0, 50, 0.0, 1.0);

    TProfile2D* h_t_RelRes_vs_xpomQ2_B0 = new TProfile2D("t_RelRes_vs_xpomQ2_B0",
        "t relres vs (x_{pom}, Q^{2}) B0;x_{pom};Q^{2} [GeV^{2}]",
        n_xpom_bins, xpom_bins, n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_t_RelRes_vs_xpomQ2_RP = new TProfile2D("t_RelRes_vs_xpomQ2_RP",
        "t relres vs (x_{pom}, Q^{2}) RP;x_{pom};Q^{2} [GeV^{2}]",
        n_xpom_bins, xpom_bins, n_bins, bin_edges_Q2.data(), "s");

    TProfile2D* h_xpom_RelRes_vs_xpomQ2_B0 = new TProfile2D("xpom_RelRes_vs_xpomQ2_B0",
        "x_{pom} relres vs (x_{pom}, Q^{2}) B0;x_{pom};Q^{2} [GeV^{2}]",
        n_xpom_bins, xpom_bins, n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_xpom_RelRes_vs_xpomQ2_RP = new TProfile2D("xpom_RelRes_vs_xpomQ2_RP",
        "x_{pom} relres vs (x_{pom}, Q^{2}) RP;x_{pom};Q^{2} [GeV^{2}]",
        n_xpom_bins, xpom_bins, n_bins, bin_edges_Q2.data(), "s");

    TProfile2D* h_beta_RelRes_vs_betaQ2_B0 = new TProfile2D("beta_RelRes_vs_betaQ2_B0",
        "#beta relres vs (beta, Q^{2}) B0;#beta;Q^{2} [GeV^{2}]",
        10, 0.0, 1.0, n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_beta_RelRes_vs_betaQ2_RP = new TProfile2D("beta_RelRes_vs_betaQ2_RP",
        "#beta relres vs (beta, Q^{2}) RP;#beta;Q^{2} [GeV^{2}]",
        10, 0.0, 1.0, n_bins, bin_edges_Q2.data(), "s");

    TProfile2D* h_xL_RelRes_vs_xLQ2_B0 = new TProfile2D("xL_RelRes_vs_xLQ2_B0",
        "x_{L} relres vs (x_{L}, Q^{2}) B0;x_{L};Q^{2} [GeV^{2}]",
        30, 0.75, 1.05, n_bins, bin_edges_Q2.data(), "s");
    TProfile2D* h_xL_RelRes_vs_xLQ2_RP = new TProfile2D("xL_RelRes_vs_xLQ2_RP",
        "x_{L} relres vs (x_{L}, Q^{2}) RP;x_{L};Q^{2} [GeV^{2}]",
        30, 0.75, 1.05, n_bins, bin_edges_Q2.data(), "s");

    // Triple differential cross section histograms
    const int n_beta_bins = 5;
    double beta_bins[n_beta_bins + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    const int n_Q2_3d_bins = 10;
    std::vector<Double_t> Q2_3d_bins = GetLogBins(1.0, 100.0, n_Q2_3d_bins);
    const int n_xpom_3d_bins = 8;
    double xpom_3d_bins[n_xpom_3d_bins + 1];
    double xpom_3d_min = 1e-3;
    double xpom_3d_max = 0.1;
    double log_3d_min = TMath::Log10(xpom_3d_min);
    double log_3d_max = TMath::Log10(xpom_3d_max);
    for (int i = 0; i <= n_xpom_3d_bins; i++) {
        xpom_3d_bins[i] = TMath::Power(10, log_3d_min + i * (log_3d_max - log_3d_min) / n_xpom_3d_bins);
    }

    TH3D* h_d3sigma_MC = new TH3D("d3sigma_dQ2dbeta_dxpom_MC",
                                  "Truth d^{3}#sigma/(dQ^{2}d#betadx_{pom});Q^{2} [GeV^{2}];#beta;x_{pom}",
                                  n_Q2_3d_bins, Q2_3d_bins.data(),
                                  n_beta_bins, beta_bins,
                                  n_xpom_3d_bins, xpom_3d_bins);
    TH3D* h_d3sigma_B0 = new TH3D("d3sigma_dQ2dbeta_dxpom_B0",
                                  "B0 Reco d^{3}#sigma/(dQ^{2}d#betadx_{pom});Q^{2} [GeV^{2}];#beta;x_{pom}",
                                  n_Q2_3d_bins, Q2_3d_bins.data(),
                                  n_beta_bins, beta_bins,
                                  n_xpom_3d_bins, xpom_3d_bins);
    TH3D* h_d3sigma_RP = new TH3D("d3sigma_dQ2dbeta_dxpom_RP",
                                  "RP Reco d^{3}#sigma/(dQ^{2}d#betadx_{pom});Q^{2} [GeV^{2}];#beta;x_{pom}",
                                  n_Q2_3d_bins, Q2_3d_bins.data(),
                                  n_beta_bins, beta_bins,
                                  n_xpom_3d_bins, xpom_3d_bins);

    // Tree reader
    TTreeReader tree_reader(events);
    TTreeReaderArray<double> mc_px_array(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py_array(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz_array(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> mc_mass_array(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> mc_genStatus_array(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> mc_pdg_array(tree_reader, "MCParticles.PDG");

    TTreeReaderArray<unsigned int> assoc_rec_id(tree_reader, "ReconstructedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> assoc_sim_id(tree_reader, "ReconstructedParticleAssociations.simID");
    TTreeReaderArray<float> re_px_array(tree_reader, "ReconstructedParticles.momentum.x");
    TTreeReaderArray<float> re_py_array(tree_reader, "ReconstructedParticles.momentum.y");
    TTreeReaderArray<float> re_pz_array(tree_reader, "ReconstructedParticles.momentum.z");
    TTreeReaderArray<float> re_energy_array(tree_reader, "ReconstructedParticles.energy");
    TTreeReaderArray<int> re_pdg_array(tree_reader, "ReconstructedParticles.PDG");
    TTreeReaderArray<int> electron_scat_index(tree_reader, "ScatteredElectronsTruth_objIdx.index");

    TTreeReaderArray<unsigned int> tsassoc_rec_id(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> tsassoc_sim_id(tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID");
    TTreeReaderArray<float> tsre_px_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> tsre_py_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> tsre_pz_array(tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");

    TTreeReaderArray<float> rp_px_array(tree_reader, "ForwardRomanPotRecParticles.momentum.x");
    TTreeReaderArray<float> rp_py_array(tree_reader, "ForwardRomanPotRecParticles.momentum.y");
    TTreeReaderArray<float> rp_pz_array(tree_reader, "ForwardRomanPotRecParticles.momentum.z");
    TTreeReaderArray<float> rp_mass_array(tree_reader, "ForwardRomanPotRecParticles.mass");
    TTreeReaderArray<int> rp_pdg_array(tree_reader, "ForwardRomanPotRecParticles.PDG");

    Float_t electron_Q2_EM, electron_Q2_DA, electron_Q2_Sigma, electron_Q2_truth;
    Float_t electron_x_EM, electron_x_DA, electron_x_Sigma, electron_x_truth;
    Float_t electron_y_EM, electron_y_DA, electron_y_Sigma, electron_y_truth;
    Float_t electron_W_EM, electron_W_DA, electron_W_Sigma, electron_W_truth;

    events->SetBranchAddress("InclusiveKinematicsElectron.Q2", &electron_Q2_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.Q2", &electron_Q2_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.Q2", &electron_Q2_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.Q2", &electron_Q2_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.x", &electron_x_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.x", &electron_x_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.x", &electron_x_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.x", &electron_x_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.y", &electron_y_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.y", &electron_y_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.y", &electron_y_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.y", &electron_y_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.W", &electron_W_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.W", &electron_W_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.W", &electron_W_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.W", &electron_W_truth);

    // Find beam particles
    std::cout << "Finding beam particles..." << std::endl;
    BeamInfo beams;
    P3MVector beame4_acc(0, 0, 0, 0), beamp4_acc(0, 0, 0, 0);
    while (tree_reader.Next()) {
        for (int i = 0; i < mc_px_array.GetSize(); i++) {
            if (mc_genStatus_array[i] != 4) continue;
            if (mc_pdg_array[i] == 2212) {
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], beams.fMass_proton);
                beamp4_acc += p;
            } else if (mc_pdg_array[i] == 11) {
                P3MVector p(mc_px_array[i], mc_py_array[i], mc_pz_array[i], beams.fMass_electron);
                beame4_acc += p;
            }
        }
    }

    auto nEntries = std::max<Long64_t>(1, events->GetEntries());
    beams.e_beam.SetCoordinates(beame4_acc.X() / nEntries, beame4_acc.Y() / nEntries, beame4_acc.Z() / nEntries, beams.fMass_electron);
    beams.p_beam.SetCoordinates(beamp4_acc.X() / nEntries, beamp4_acc.Y() / nEntries, beamp4_acc.Z() / nEntries, beams.fMass_proton);
    std::cout << "Found beam energies " << beams.e_beam.E() << "x" << beams.p_beam.E() << " GeV" << std::endl;
    undoAfterburnAndCalc(beams.p_beam, beams.e_beam);
    const double s_cm = (beams.p_beam + beams.e_beam).M2();
    tree_reader.Restart();

    std::vector<double> t_truth_vec_B0;
    std::vector<double> t_reco_vec_B0;
    std::vector<double> t_truth_vec_RP;
    std::vector<double> t_reco_vec_RP;

    // Event loop
    Long64_t nentries = events->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        if (i % 1000 == 0) {
            printf("\rProcessing event %lld of %lld; %.2f percent done.", i, nentries, 100.0 * i / nentries);
            fflush(stdout);
        }
        tree_reader.Next();
        events->GetEntry(i);

        // Fill Q2/x/y histograms
        h_x_EM->Fill(electron_x_EM);
        h_x_DA->Fill(electron_x_DA);
        h_x_Sigma->Fill(electron_x_Sigma);
        h_x_truth->Fill(electron_x_truth);

        h_y_EM->Fill(electron_y_EM);
        h_y_DA->Fill(electron_y_DA);
        h_y_Sigma->Fill(electron_y_Sigma);
        h_y_truth->Fill(electron_y_truth);

        if (electron_x_truth != 0.0f) {
            double res_x_EM = (electron_x_EM - electron_x_truth) / electron_x_truth;
            double res_x_DA = (electron_x_DA - electron_x_truth) / electron_x_truth;
            double res_x_Sigma = (electron_x_Sigma - electron_x_truth) / electron_x_truth;
            h_RelRes_x_EM->Fill(res_x_EM);
            h_RelRes_x_DA->Fill(res_x_DA);
            h_RelRes_x_Sigma->Fill(res_x_Sigma);
            h_RelRes_x_binned_EM->Fill(electron_x_truth, res_x_EM);
            h_RelRes_x_binned_DA->Fill(electron_x_truth, res_x_DA);
            h_RelRes_x_binned_Sigma->Fill(electron_x_truth, res_x_Sigma);
            if (electron_Q2_truth > 0.0f) {
                h_x_RelRes_vs_xQ2_EM->Fill(electron_x_truth, electron_Q2_truth, res_x_EM);
                h_x_RelRes_vs_xQ2_DA->Fill(electron_x_truth, electron_Q2_truth, res_x_DA);
                h_x_RelRes_vs_xQ2_Sigma->Fill(electron_x_truth, electron_Q2_truth, res_x_Sigma);
            }
        }
        if (electron_y_truth != 0.0f) {
            double res_y_EM = (electron_y_EM - electron_y_truth) / electron_y_truth;
            double res_y_DA = (electron_y_DA - electron_y_truth) / electron_y_truth;
            double res_y_Sigma = (electron_y_Sigma - electron_y_truth) / electron_y_truth;
            h_RelRes_y_EM->Fill(res_y_EM);
            h_RelRes_y_DA->Fill(res_y_DA);
            h_RelRes_y_Sigma->Fill(res_y_Sigma);
            h_RelRes_y_binned_EM->Fill(electron_y_truth, res_y_EM);
            h_RelRes_y_binned_DA->Fill(electron_y_truth, res_y_DA);
            h_RelRes_y_binned_Sigma->Fill(electron_y_truth, res_y_Sigma);
            if (electron_Q2_truth > 0.0f && electron_x_truth > 0.0f) {
                h_y_RelRes_vs_xQ2_EM->Fill(electron_x_truth, electron_Q2_truth, res_y_EM);
                h_y_RelRes_vs_xQ2_DA->Fill(electron_x_truth, electron_Q2_truth, res_y_DA);
                h_y_RelRes_vs_xQ2_Sigma->Fill(electron_x_truth, electron_Q2_truth, res_y_Sigma);
            }
        }

        h_Corr_x_EM->Fill(electron_x_truth, electron_x_EM);
        h_Corr_x_DA->Fill(electron_x_truth, electron_x_DA);
        h_Corr_x_Sigma->Fill(electron_x_truth, electron_x_Sigma);

        h_Corr_y_EM->Fill(electron_y_truth, electron_y_EM);
        h_Corr_y_DA->Fill(electron_y_truth, electron_y_DA);
        h_Corr_y_Sigma->Fill(electron_y_truth, electron_y_Sigma);

        h_Q2_EM->Fill(electron_Q2_EM);
        h_Q2_DA->Fill(electron_Q2_DA);
        h_Q2_Sigma->Fill(electron_Q2_Sigma);
        h_Q2_EICRecon->Fill(electron_Q2_EM);
        if (electron_x_EM > 0.0f && electron_y_EM > 0.0f) {
            double Q2_calc = s_cm * electron_x_EM * electron_y_EM;
            h_Q2_calc->Fill(Q2_calc);
        }

        out_x_truth = electron_x_truth;
        out_x_EM = electron_x_EM;
        out_x_DA = electron_x_DA;
        out_x_Sigma = electron_x_Sigma;
        out_y_truth = electron_y_truth;
        out_y_EM = electron_y_EM;
        out_y_DA = electron_y_DA;
        out_y_Sigma = electron_y_Sigma;
        tree->Fill();

        h_Q2_truth->Fill(electron_Q2_truth);

        if (electron_Q2_truth != 0.0f) {
            double res_Q2_EM = (electron_Q2_EM - electron_Q2_truth) / electron_Q2_truth;
            double res_Q2_DA = (electron_Q2_DA - electron_Q2_truth) / electron_Q2_truth;
            double res_Q2_Sigma = (electron_Q2_Sigma - electron_Q2_truth) / electron_Q2_truth;
            h_RelRes_Q2_EM->Fill(res_Q2_EM);
            h_RelRes_Q2_DA->Fill(res_Q2_DA);
            h_RelRes_Q2_Sigma->Fill(res_Q2_Sigma);
            h_RelRes_Q2_binned_EM->Fill(electron_Q2_truth, res_Q2_EM);
            h_RelRes_Q2_binned_DA->Fill(electron_Q2_truth, res_Q2_DA);
            h_RelRes_Q2_binned_Sigma->Fill(electron_Q2_truth, res_Q2_Sigma);
            if (electron_x_truth > 0.0f && electron_y_truth > 0.0f) {
                h_Q2_RelRes_vs_xy_EM->Fill(electron_x_truth, electron_y_truth, res_Q2_EM);
                h_Q2_RelRes_vs_xy_DA->Fill(electron_x_truth, electron_y_truth, res_Q2_DA);
                h_Q2_RelRes_vs_xy_Sigma->Fill(electron_x_truth, electron_y_truth, res_Q2_Sigma);
            }
        }

        h_Corr_Q2_EM->Fill(electron_Q2_truth, electron_Q2_EM);
        h_Corr_Q2_DA->Fill(electron_Q2_truth, electron_Q2_DA);
        h_Corr_Q2_Sigma->Fill(electron_Q2_truth, electron_Q2_Sigma);

        // E-pz for matched particles
        double sumEPz_truth_matched = 0.0;
        double sumEPz_reco_matched = 0.0;
        CalculateSumEPz_Matched(re_px_array, re_py_array, re_pz_array, re_energy_array,
                                mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
                                assoc_rec_id, assoc_sim_id, sumEPz_truth_matched, sumEPz_reco_matched);
        h_EPz_truth->Fill(sumEPz_truth_matched);
        h_EPz->Fill(sumEPz_reco_matched);
        h_EPz_2D->Fill(sumEPz_truth_matched, sumEPz_reco_matched);

        // Eta max and MX2
        double eta_max_reco = -999.0;
        double eta_max_truth = -999.0;
        for (unsigned int j = 0; j < re_energy_array.GetSize(); j++) {
            P3MVector particle_reco(re_px_array[j], re_py_array[j], re_pz_array[j], 0.0);
            double eta_reco = particle_reco.Eta();
            if (eta_reco > eta_max_reco) eta_max_reco = eta_reco;

            int mc_idx = -1;
            for (unsigned int k = 0; k < assoc_rec_id.GetSize(); k++) {
                if (assoc_rec_id[k] == j) {
                    mc_idx = assoc_sim_id[k];
                    break;
                }
            }
            if (mc_idx >= 0 && mc_idx < mc_px_array.GetSize()) {
                P3MVector particle_mc(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
                double eta_mc = particle_mc.Eta();
                if (eta_mc > eta_max_truth) eta_max_truth = eta_mc;
            }
        }
        if (eta_max_reco > -999.0) h_eta_max->Fill(eta_max_reco);
        if (eta_max_truth > -999.0) h_eta_max_truth->Fill(eta_max_truth);

        int scat_e_idx = -1;
        if (electron_scat_index.GetSize() > 0) {
            scat_e_idx = electron_scat_index[0];
        }

        P3EVector total_hadrons_reco(0.0, 0.0, 0.0, 0.0);
        for (unsigned int j = 0; j < re_energy_array.GetSize(); j++) {
            if ((int)j == scat_e_idx) continue;
            P3EVector particle(re_px_array[j], re_py_array[j], re_pz_array[j], re_energy_array[j]);
            total_hadrons_reco += particle;
        }
        double MX2_reco = total_hadrons_reco.M2();

        P3EVector total_hadrons_truth(0.0, 0.0, 0.0, 0.0);
        for (unsigned int j = 0; j < re_energy_array.GetSize(); j++) {
            if ((int)j == scat_e_idx) continue;
            int mc_idx = -1;
            for (unsigned int k = 0; k < assoc_rec_id.GetSize(); k++) {
                if (assoc_rec_id[k] == j) {
                    mc_idx = assoc_sim_id[k];
                    break;
                }
            }
            if (mc_idx >= 0 && mc_idx < mc_px_array.GetSize()) {
                double px_mc = mc_px_array[mc_idx];
                double py_mc = mc_py_array[mc_idx];
                double pz_mc = mc_pz_array[mc_idx];
                double m_mc = mc_mass_array[mc_idx];
                double E_mc = TMath::Sqrt(px_mc * px_mc + py_mc * py_mc + pz_mc * pz_mc + m_mc * m_mc);
                P3EVector particle_mc(px_mc, py_mc, pz_mc, E_mc);
                total_hadrons_truth += particle_mc;
            }
        }
        double MX2_truth = total_hadrons_truth.M2();

        h_MX2->Fill(MX2_reco);
        h_MX2_truth->Fill(MX2_truth);
        h_MX2_MC->Fill(MX2_truth);
        h_MX2_eX->Fill(MX2_reco);
        h_MX2_corr->Fill(MX2_truth, MX2_reco);
        if (MX2_truth > 0.0) {
            double relres_MX2 = (MX2_reco - MX2_truth) / MX2_truth;
            h_MX2_RelRes_binned->Fill(MX2_truth, relres_MX2);
            if (electron_Q2_truth > 0.0) {
                h_MX2_RelRes_vs_MX2Q2->Fill(MX2_truth, electron_Q2_truth, relres_MX2);
            }
        }

        // Build truth proton list for RP matching
        std::vector<P3MVector> truth_protons;
        truth_protons.reserve(mc_px_array.GetSize());
        for (int j = 0; j < mc_px_array.GetSize(); j++) {
            if (mc_genStatus_array[j] == 1 && mc_pdg_array[j] == 2212) {
                P3MVector p(mc_px_array[j], mc_py_array[j], mc_pz_array[j], mc_mass_array[j]);
                undoAfterburn(p);
                truth_protons.push_back(p);

                double xL_truth_pz = p.Pz() / beams.p_beam.Pz();
                double xpom_from_xL = 1.0 - xL_truth_pz;
                double t_val = TMath::Abs(CalcT(beams.p_beam, p));
                double xpom_from_def = XpomFromDefinition(MX2_truth, electron_Q2_truth, t_val, electron_W_truth);

                h_theta_MC->Fill(p.Theta() * 1000.0);
                h_t_MC->Fill(t_val);
                h_dsigma_dt_MC->Fill(t_val);
                h_xL_MC->Fill(xL_truth_pz);
                h_xpom_MC->Fill(xpom_from_xL);

                if (ValidUnitInterval(xpom_from_def)) {
                    h_xpom_def_MC->Fill(xpom_from_def);
                    h_xpom_comp_MC->Fill(xpom_from_xL, xpom_from_def);
                    if (ValidUnitInterval(electron_x_truth)) {
                        double beta_truth = electron_x_truth / xpom_from_def;
                        if (ValidUnitInterval(beta_truth)) {
                            h_beta_MC->Fill(beta_truth);
                            if (electron_Q2_truth > 0.0) {
                                h_d3sigma_MC->Fill(electron_Q2_truth, beta_truth, xpom_from_def);
                            }
                        }
                    }
                }
            }
        }

        // Process B0 protons (truth-seeded)
        for (unsigned int j = 0; j < tsassoc_rec_id.GetSize(); j++) {
            auto mc_idx = tsassoc_sim_id[j];
            if (mc_idx >= (unsigned)mc_pdg_array.GetSize() || mc_genStatus_array[mc_idx] != 1 || mc_pdg_array[mc_idx] != 2212)
                continue;

            P3MVector p_reco(tsre_px_array[j], tsre_py_array[j], tsre_pz_array[j], mc_mass_array[mc_idx]);
            undoAfterburn(p_reco);
            h_theta_all_TS->Fill(p_reco.Theta() * 1000.0);

            if (p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;
            h_theta_B0->Fill(p_reco.Theta() * 1000.0);

            P3MVector p_truth(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
            undoAfterburn(p_truth);

            double t_reco = CalcT(beams.p_beam, p_reco);
            double t_truth = CalcT(beams.p_beam, p_truth);
            double t_reco_abs = TMath::Abs(t_reco);
            double t_truth_abs = TMath::Abs(t_truth);

            double xL_reco_pz = p_reco.Pz() / beams.p_beam.Pz();
            double xL_truth_pz = p_truth.Pz() / beams.p_beam.Pz();
            double xpom_reco_from_xL = 1.0 - xL_reco_pz;
            double xpom_truth_from_xL = 1.0 - xL_truth_pz;

            double xpom_reco_from_def = XpomFromDefinition(MX2_reco, electron_Q2_EM, t_reco_abs, electron_W_EM);
            double xpom_truth_from_def = XpomFromDefinition(MX2_truth, electron_Q2_truth, t_truth_abs, electron_W_truth);

            h_t_B0->Fill(t_reco_abs);
            h_dsigma_dt_B0->Fill(t_reco_abs);
            h_t_corr_B0->Fill(t_truth_abs, t_reco_abs);
            t_truth_vec_B0.push_back(t_truth_abs);
            t_reco_vec_B0.push_back(t_reco_abs);
            if (t_truth_abs > 1e-6) {
                double relres_t = (t_reco_abs - t_truth_abs) / t_truth_abs;
                h_t_res_B0->Fill(relres_t);
                h_t_RelRes_binned_B0->Fill(t_truth_abs, relres_t);
                if (ValidUnitInterval(xpom_truth_from_xL) && electron_Q2_truth > 0.0) {
                    h_t_RelRes_vs_xpomQ2_B0->Fill(xpom_truth_from_xL, electron_Q2_truth, relres_t);
                }
            }

            h_xL_B0->Fill(xL_reco_pz);
            h_xL_corr_B0->Fill(xL_truth_pz, xL_reco_pz);
            if (xL_truth_pz > 1e-6) {
                double relres_xL = (xL_reco_pz - xL_truth_pz) / xL_truth_pz;
                h_xL_res_B0->Fill(relres_xL);
                h_xL_RelRes_binned_B0->Fill(xL_truth_pz, relres_xL);
                if (electron_Q2_truth > 0.0) {
                    h_xL_RelRes_vs_xLQ2_B0->Fill(xL_truth_pz, electron_Q2_truth, relres_xL);
                }
            }

            h_xpom_B0->Fill(xpom_reco_from_xL);
            h_xpom_corr_B0->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
            if (xpom_truth_from_xL > 1e-6) {
                double relres_xpom = (xpom_reco_from_xL - xpom_truth_from_xL) / xpom_truth_from_xL;
                h_xpom_res_B0->Fill(relres_xpom);
                h_xpom_RelRes_binned_B0->Fill(xpom_truth_from_xL, relres_xpom);
                if (electron_Q2_truth > 0.0) {
                    h_xpom_RelRes_vs_xpomQ2_B0->Fill(xpom_truth_from_xL, electron_Q2_truth, relres_xpom);
                }
            }

            if (ValidUnitInterval(xpom_reco_from_def)) {
                h_xpom_def_B0->Fill(xpom_reco_from_def);
                h_xpom_comp_B0->Fill(xpom_reco_from_xL, xpom_reco_from_def);

                if (ValidUnitInterval(electron_x_EM)) {
                    double beta_reco = electron_x_EM / xpom_reco_from_def;
                    if (ValidUnitInterval(beta_reco)) {
                        h_beta_B0->Fill(beta_reco);
                        if (electron_Q2_EM > 0.0) {
                            h_d3sigma_B0->Fill(electron_Q2_EM, beta_reco, xpom_reco_from_def);
                        }
                        h_beta_vs_Q2->Fill(electron_Q2_EM, beta_reco);
                        h_beta_vs_xpom->Fill(xpom_reco_from_def, beta_reco);
                        h_beta_vs_t->Fill(t_reco_abs, beta_reco);

                        if (ValidUnitInterval(xpom_truth_from_def) && ValidUnitInterval(electron_x_truth)) {
                            double beta_truth = electron_x_truth / xpom_truth_from_def;
                            if (ValidUnitInterval(beta_truth)) {
                                double relres_beta = (beta_reco - beta_truth) / beta_truth;
                                h_beta_corr_B0->Fill(beta_truth, beta_reco);
                                h_beta_res_B0->Fill(relres_beta);
                                h_beta_RelRes_binned_B0->Fill(beta_truth, relres_beta);
                                if (electron_Q2_truth > 0.0) {
                                    h_beta_RelRes_vs_betaQ2_B0->Fill(beta_truth, electron_Q2_truth, relres_beta);
                                }
                            }
                        }
                    }
                }
            }

            if (t_reco_abs > t_bins[1]) {
                h_t_B0_cutFirstBin->Fill(t_reco_abs);
                h_t_corr_B0_cutFirstBin->Fill(t_truth_abs, t_reco_abs);
                h_xL_B0_cutFirstBin->Fill(xL_reco_pz);
                h_xL_corr_B0_cutFirstBin->Fill(xL_truth_pz, xL_reco_pz);
                h_xpom_B0_cutFirstBin->Fill(xpom_reco_from_xL);
                h_xpom_corr_B0_cutFirstBin->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
            }
        }

        // Process RP protons
        for (int j = 0; j < rp_px_array.GetSize(); j++) {
            if (rp_pdg_array[j] != 2212) continue;

            P3MVector p_rp(rp_px_array[j], rp_py_array[j], rp_pz_array[j], rp_mass_array[j]);
            if (p_rp.Theta() >= 0.005) continue;
            h_theta_RP->Fill(p_rp.Theta() * 1000.0);

            double min_dtheta = 1e9;
            int best_match = -1;
            for (size_t k = 0; k < truth_protons.size(); k++) {
                if (truth_protons[k].Theta() >= 0.006) continue;
                double dtheta = TMath::Abs(p_rp.Theta() - truth_protons[k].Theta());
                if (dtheta < min_dtheta) {
                    min_dtheta = dtheta;
                    best_match = static_cast<int>(k);
                }
            }
            if (best_match < 0 || min_dtheta >= 0.001) continue;

            const P3MVector& p_truth = truth_protons[best_match];
            double t_reco = CalcT(beams.p_beam, p_rp);
            double t_truth = CalcT(beams.p_beam, p_truth);
            double t_reco_abs = TMath::Abs(t_reco);
            double t_truth_abs = TMath::Abs(t_truth);

            double xL_reco_pz = p_rp.Pz() / beams.p_beam.Pz();
            double xL_truth_pz = p_truth.Pz() / beams.p_beam.Pz();
            double xpom_reco_from_xL = 1.0 - xL_reco_pz;
            double xpom_truth_from_xL = 1.0 - xL_truth_pz;

            double xpom_reco_from_def = XpomFromDefinition(MX2_reco, electron_Q2_EM, t_reco_abs, electron_W_EM);
            double xpom_truth_from_def = XpomFromDefinition(MX2_truth, electron_Q2_truth, t_truth_abs, electron_W_truth);

            h_t_RP_histo->Fill(t_reco_abs);
            h_dsigma_dt_RP->Fill(t_reco_abs);
            h_t_corr_RP->Fill(t_truth_abs, t_reco_abs);
            t_truth_vec_RP.push_back(t_truth_abs);
            t_reco_vec_RP.push_back(t_reco_abs);
            if (t_truth_abs > 1e-6) {
                double relres_t = (t_reco_abs - t_truth_abs) / t_truth_abs;
                h_t_res_RP->Fill(relres_t);
                h_t_RelRes_binned_RP->Fill(t_truth_abs, relres_t);
                if (ValidUnitInterval(xpom_truth_from_xL) && electron_Q2_truth > 0.0) {
                    h_t_RelRes_vs_xpomQ2_RP->Fill(xpom_truth_from_xL, electron_Q2_truth, relres_t);
                }
            }

            h_xL_RP->Fill(xL_reco_pz);
            h_xL_corr_RP->Fill(xL_truth_pz, xL_reco_pz);
            if (xL_truth_pz > 1e-6) {
                double relres_xL = (xL_reco_pz - xL_truth_pz) / xL_truth_pz;
                h_xL_res_RP->Fill(relres_xL);
                h_xL_RelRes_binned_RP->Fill(xL_truth_pz, relres_xL);
                if (electron_Q2_truth > 0.0) {
                    h_xL_RelRes_vs_xLQ2_RP->Fill(xL_truth_pz, electron_Q2_truth, relres_xL);
                }
            }

            h_xpom_RP->Fill(xpom_reco_from_xL);
            h_xpom_corr_RP->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
            if (xpom_truth_from_xL > 1e-6) {
                double relres_xpom = (xpom_reco_from_xL - xpom_truth_from_xL) / xpom_truth_from_xL;
                h_xpom_res_RP->Fill(relres_xpom);
                h_xpom_RelRes_binned_RP->Fill(xpom_truth_from_xL, relres_xpom);
                if (electron_Q2_truth > 0.0) {
                    h_xpom_RelRes_vs_xpomQ2_RP->Fill(xpom_truth_from_xL, electron_Q2_truth, relres_xpom);
                }
            }

            if (ValidUnitInterval(xpom_reco_from_def)) {
                h_xpom_def_RP->Fill(xpom_reco_from_def);
                h_xpom_comp_RP->Fill(xpom_reco_from_xL, xpom_reco_from_def);

                if (ValidUnitInterval(electron_x_EM)) {
                    double beta_reco = electron_x_EM / xpom_reco_from_def;
                    if (ValidUnitInterval(beta_reco)) {
                        h_beta_RP->Fill(beta_reco);
                        if (electron_Q2_EM > 0.0) {
                            h_d3sigma_RP->Fill(electron_Q2_EM, beta_reco, xpom_reco_from_def);
                        }
                        if (ValidUnitInterval(xpom_truth_from_def) && ValidUnitInterval(electron_x_truth)) {
                            double beta_truth = electron_x_truth / xpom_truth_from_def;
                            if (ValidUnitInterval(beta_truth)) {
                                double relres_beta = (beta_reco - beta_truth) / beta_truth;
                                h_beta_corr_RP->Fill(beta_truth, beta_reco);
                                h_beta_res_RP->Fill(relres_beta);
                                h_beta_RelRes_binned_RP->Fill(beta_truth, relres_beta);
                                if (electron_Q2_truth > 0.0) {
                                    h_beta_RelRes_vs_betaQ2_RP->Fill(beta_truth, electron_Q2_truth, relres_beta);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "\nDone looping over events." << std::endl;

    // Scale differential cross sections
    const double sigma_total = 4.22;
    const double N_gen = 100000.0;
    const double scale_factor = sigma_total / N_gen;

    for (int i = 1; i <= h_dsigma_dt_MC->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_MC->GetBinContent(i);
        double bin_error = h_dsigma_dt_MC->GetBinError(i);
        double bin_width = h_dsigma_dt_MC->GetBinWidth(i);
        h_dsigma_dt_MC->SetBinContent(i, (bin_content * scale_factor) / bin_width);
        h_dsigma_dt_MC->SetBinError(i, (bin_error * scale_factor) / bin_width);
    }
    for (int i = 1; i <= h_dsigma_dt_B0->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_B0->GetBinContent(i);
        double bin_error = h_dsigma_dt_B0->GetBinError(i);
        double bin_width = h_dsigma_dt_B0->GetBinWidth(i);
        h_dsigma_dt_B0->SetBinContent(i, (bin_content * scale_factor) / bin_width);
        h_dsigma_dt_B0->SetBinError(i, (bin_error * scale_factor) / bin_width);
    }
    for (int i = 1; i <= h_dsigma_dt_RP->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_RP->GetBinContent(i);
        double bin_error = h_dsigma_dt_RP->GetBinError(i);
        double bin_width = h_dsigma_dt_RP->GetBinWidth(i);
        h_dsigma_dt_RP->SetBinContent(i, (bin_content * scale_factor) / bin_width);
        h_dsigma_dt_RP->SetBinError(i, (bin_error * scale_factor) / bin_width);
    }

    // Scale triple differential cross sections
    for (int i = 1; i <= h_d3sigma_MC->GetNbinsX(); i++) {
        for (int j = 1; j <= h_d3sigma_MC->GetNbinsY(); j++) {
            for (int k = 1; k <= h_d3sigma_MC->GetNbinsZ(); k++) {
                double bin_content = h_d3sigma_MC->GetBinContent(i, j, k);
                double bin_error = h_d3sigma_MC->GetBinError(i, j, k);
                double width_Q2 = h_d3sigma_MC->GetXaxis()->GetBinWidth(i);
                double width_beta = h_d3sigma_MC->GetYaxis()->GetBinWidth(j);
                double width_xpom = h_d3sigma_MC->GetZaxis()->GetBinWidth(k);
                double bin_volume = width_Q2 * width_beta * width_xpom;
                h_d3sigma_MC->SetBinContent(i, j, k, (bin_content * scale_factor) / bin_volume);
                h_d3sigma_MC->SetBinError(i, j, k, (bin_error * scale_factor) / bin_volume);
            }
        }
    }
    for (int i = 1; i <= h_d3sigma_B0->GetNbinsX(); i++) {
        for (int j = 1; j <= h_d3sigma_B0->GetNbinsY(); j++) {
            for (int k = 1; k <= h_d3sigma_B0->GetNbinsZ(); k++) {
                double bin_content = h_d3sigma_B0->GetBinContent(i, j, k);
                double bin_error = h_d3sigma_B0->GetBinError(i, j, k);
                double width_Q2 = h_d3sigma_B0->GetXaxis()->GetBinWidth(i);
                double width_beta = h_d3sigma_B0->GetYaxis()->GetBinWidth(j);
                double width_xpom = h_d3sigma_B0->GetZaxis()->GetBinWidth(k);
                double bin_volume = width_Q2 * width_beta * width_xpom;
                h_d3sigma_B0->SetBinContent(i, j, k, (bin_content * scale_factor) / bin_volume);
                h_d3sigma_B0->SetBinError(i, j, k, (bin_error * scale_factor) / bin_volume);
            }
        }
    }
    for (int i = 1; i <= h_d3sigma_RP->GetNbinsX(); i++) {
        for (int j = 1; j <= h_d3sigma_RP->GetNbinsY(); j++) {
            for (int k = 1; k <= h_d3sigma_RP->GetNbinsZ(); k++) {
                double bin_content = h_d3sigma_RP->GetBinContent(i, j, k);
                double bin_error = h_d3sigma_RP->GetBinError(i, j, k);
                double width_Q2 = h_d3sigma_RP->GetXaxis()->GetBinWidth(i);
                double width_beta = h_d3sigma_RP->GetYaxis()->GetBinWidth(j);
                double width_xpom = h_d3sigma_RP->GetZaxis()->GetBinWidth(k);
                double bin_volume = width_Q2 * width_beta * width_xpom;
                h_d3sigma_RP->SetBinContent(i, j, k, (bin_content * scale_factor) / bin_volume);
                h_d3sigma_RP->SetBinError(i, j, k, (bin_error * scale_factor) / bin_volume);
            }
        }
    }

    outputFile->cd();
    tree->Write();

    // Q2/x/y histograms
    h_RelRes_Q2_EM->Write();
    h_RelRes_Q2_DA->Write();
    h_RelRes_Q2_Sigma->Write();
    h_RelRes_Q2_binned_EM->Write();
    h_RelRes_Q2_binned_DA->Write();
    h_RelRes_Q2_binned_Sigma->Write();
    h_Corr_Q2_EM->Write();
    h_Corr_Q2_DA->Write();
    h_Corr_Q2_Sigma->Write();
    h_Q2_truth->Write();
    h_Q2_EM->Write();
    h_Q2_DA->Write();
    h_Q2_Sigma->Write();
    h_Q2_EICRecon->Write();
    h_Q2_calc->Write();

    h_x_truth->Write();
    h_x_EM->Write();
    h_x_DA->Write();
    h_x_Sigma->Write();
    h_y_truth->Write();
    h_y_EM->Write();
    h_y_DA->Write();
    h_y_Sigma->Write();

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

    h_Q2_RelRes_vs_xy_EM->Write();
    h_Q2_RelRes_vs_xy_DA->Write();
    h_Q2_RelRes_vs_xy_Sigma->Write();
    h_x_RelRes_vs_xQ2_EM->Write();
    h_x_RelRes_vs_xQ2_DA->Write();
    h_x_RelRes_vs_xQ2_Sigma->Write();
    h_y_RelRes_vs_xQ2_EM->Write();
    h_y_RelRes_vs_xQ2_DA->Write();
    h_y_RelRes_vs_xQ2_Sigma->Write();

    h_EPz_truth->Write();
    h_EPz->Write();
    h_EPz_2D->Write();
    h_eta_max->Write();
    h_eta_max_truth->Write();
    h_MX2->Write();
    h_MX2_truth->Write();
    h_MX2_MC->Write();
    h_MX2_eX->Write();
    h_MX2_corr->Write();
    h_MX2_RelRes_binned->Write();
    h_MX2_RelRes_vs_MX2Q2->Write();

    // t-analysis histograms
    h_t_MC->Write();
    h_t_B0->Write();
    h_t_RP_histo->Write();
    h_t_B0_cutFirstBin->Write();
    h_dsigma_dt_MC->Write();
    h_dsigma_dt_B0->Write();
    h_dsigma_dt_RP->Write();
    h_t_res_B0->Write();
    h_t_res_RP->Write();
    h_t_corr_B0->Write();
    h_t_corr_RP->Write();
    h_t_corr_B0_cutFirstBin->Write();
    h_t_RelRes_binned_B0->Write();
    h_t_RelRes_binned_RP->Write();

    h_theta_MC->Write();
    h_theta_all_TS->Write();
    h_theta_B0->Write();
    h_theta_RP->Write();

    h_xL_MC->Write();
    h_xL_B0->Write();
    h_xL_RP->Write();
    h_xL_B0_cutFirstBin->Write();
    h_xL_res_B0->Write();
    h_xL_res_RP->Write();
    h_xL_corr_B0->Write();
    h_xL_corr_RP->Write();
    h_xL_corr_B0_cutFirstBin->Write();
    h_xL_RelRes_binned_B0->Write();
    h_xL_RelRes_binned_RP->Write();

    h_xpom_MC->Write();
    h_xpom_B0->Write();
    h_xpom_RP->Write();
    h_xpom_B0_cutFirstBin->Write();
    h_xpom_def_MC->Write();
    h_xpom_def_B0->Write();
    h_xpom_def_RP->Write();
    h_xpom_comp_MC->Write();
    h_xpom_comp_B0->Write();
    h_xpom_comp_RP->Write();
    h_xpom_res_B0->Write();
    h_xpom_res_RP->Write();
    h_xpom_corr_B0->Write();
    h_xpom_corr_RP->Write();
    h_xpom_corr_B0_cutFirstBin->Write();
    h_xpom_RelRes_binned_B0->Write();
    h_xpom_RelRes_binned_RP->Write();

    h_beta_MC->Write();
    h_beta_B0->Write();
    h_beta_RP->Write();
    h_beta_res_B0->Write();
    h_beta_res_RP->Write();
    h_beta_corr_B0->Write();
    h_beta_corr_RP->Write();
    h_beta_RelRes_binned_B0->Write();
    h_beta_RelRes_binned_RP->Write();
    h_beta_vs_Q2->Write();
    h_beta_vs_xpom->Write();
    h_beta_vs_t->Write();

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

    // Combined correlation graphs
    if (!t_truth_vec_B0.empty()) {
        TGraph* g_B0 = new TGraph(t_truth_vec_B0.size(), t_truth_vec_B0.data(), t_reco_vec_B0.data());
        g_B0->SetName("t_corr_B0_graph");
        g_B0->Write();
        delete g_B0;
    }
    if (!t_truth_vec_RP.empty()) {
        TGraph* g_RP = new TGraph(t_truth_vec_RP.size(), t_truth_vec_RP.data(), t_reco_vec_RP.data());
        g_RP->SetName("t_corr_RP_graph");
        g_RP->Write();
        delete g_RP;
    }

    outputFile->Close();
    delete events;
    delete outputFile;

    std::cout << "\nAnalysis complete! Output saved to " << outputName << std::endl;
    return 0;
}
