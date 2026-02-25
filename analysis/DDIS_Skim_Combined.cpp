// Combined DDIS Skimmer: Q2/xy analysis + Mandelstam t analysis + x_pom comparison
// g++ DDIS_Skim_Combined.cpp -o DDIS_Skim_Combined $(root-config --cflags --glibs)
// ./DDIS_Skim_Combined filelist.txt

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TObject.h>
#include <TString.h>
#include <TProfile2D.h>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
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

// Function: E-pz calculation for matched particles
void CalculateSumEPz_Matched(
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
    // Initialize output sums
    sumEPz_truth = 0.0;
    sumEPz_reco = 0.0;

    // Loop over all reconstructed particles
    for(unsigned int i = 0; i < re_energy.GetSize(); i++){
        // Search for this reco particle in associations
        int mc_idx = -1;
        for(unsigned int j = 0; j < assoc_rec_id.GetSize(); j++){
            if(assoc_rec_id[j] == i) {
                mc_idx = assoc_sim_id[j];
                break;
            }
        }

        // If no association found, skip this reco particle
        if(mc_idx < 0) continue;

        // Calculate reco E - pz
        double E_reco = re_energy[i];
        double pz_reco = re_pz[i];
        sumEPz_reco += (E_reco - pz_reco);

        // Calculate MC E - pz for the matched particle
        double px_mc = mc_px[mc_idx];
        double py_mc = mc_py[mc_idx];
        double pz_mc = mc_pz[mc_idx];
        double m_mc = mc_mass[mc_idx];
        double E_mc = TMath::Sqrt(px_mc*px_mc + py_mc*py_mc + pz_mc*pz_mc + m_mc*m_mc);
        sumEPz_truth += (E_mc - pz_mc);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fileList.txt>" << std::endl;
        return 1;
    }
    TString fileList = argv[1];

    std::cout<< " __ __ __ __ __ __ __ __ __ __" <<std::endl;
    std::cout<< "|                             |"<<std::endl;
    std::cout<< "|  Combined ePIC DDIS Skim    |"<<std::endl;
    std::cout<< "|  Q2/xy + t + x_pom comp.    |"<<std::endl;
    std::cout<< "|__ __ __ __ __ __ __ __ __ __|"<<std::endl;
    std::cout<< "\nInput filelist: " << fileList <<std::endl;

    std::ifstream fileListStream;
    fileListStream.open(fileList);
    std::string fileName;

    //---------------------------------------------------------
    // CREATE TCHAIN AND OUTPUT ROOT FILE
    //---------------------------------------------------------
    TChain* events = new TChain("events");
    Int_t nFiles{0};
    std::vector<std::string> addedFiles;
    while(getline(fileListStream, fileName)){
        TString tmp = fileName;
        if (!std::filesystem::exists(tmp.Data())) {
        std::cerr << "\nError: File does not exist: " << fileName << std::endl;
        continue;
        }
        auto inputRootFile = TFile::Open(tmp);
        events->Add((TString)fileName);
        inputRootFile->Close();
        nFiles++;
        addedFiles.push_back(fileName);
    }
    std::cout<<"\nNo. of files: "<<nFiles<<"; no. of events: "<<events->GetEntries()<<std::endl;
    
    // Create output file
    TFile* outputFile = new TFile("DDIS_Combined_output.root", "RECREATE");

    // Create TTree for event-level data
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

    // Correlation plots (truth vs reco)
    TH2D* h_Corr_x_EM     = new TH2D("x_Corr_EM",     "x correlation (EM);x_{truth};x_{EM}",           x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_DA     = new TH2D("x_Corr_DA",     "x correlation (DA);x_{truth};x_{DA}",           x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());
    TH2D* h_Corr_x_Sigma = new TH2D("x_Corr_Sigma", "x correlation (Sigma);x_{truth};x_{Sigma}",   x_bins.size()-1, x_bins.data(), x_bins.size()-1, x_bins.data());

    TH2D* h_Corr_y_EM     = new TH2D("y_Corr_EM",     "y correlation (EM);y_{truth};y_{EM}",           n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_DA     = new TH2D("y_Corr_DA",     "y correlation (DA);y_{truth};y_{DA}",           n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);
    TH2D* h_Corr_y_Sigma = new TH2D("y_Corr_Sigma", "y correlation (Sigma);y_{truth};y_{Sigma}",   n_y_bins, 0.0, 1.0, n_y_bins, 0.0, 1.0);

    // TProfile2D for Q2 resolution - BINNED in (x_Bj, y) but DISPLAYED in (x_Bj, Q2)
    // Since Q2 = s*x*y, binning in (x,y) creates a parallelogram in (x,Q2) space
    TProfile2D* h_Q2_RelRes_vs_xy_EM = new TProfile2D("Q2_RelRes_vs_xy_EM",
        "Q^{2} Rel. Res. binned in (x_{Bj}, y) (EM);x_{Bj};y",
        x_bins.size()-1, x_bins.data(), n_y_bins, 0.0, 1.0, "s");
    TProfile2D* h_Q2_RelRes_vs_xy_DA = new TProfile2D("Q2_RelRes_vs_xy_DA",
        "Q^{2} Rel. Res. binned in (x_{Bj}, y) (DA);x_{Bj};y",
        x_bins.size()-1, x_bins.data(), n_y_bins, 0.0, 1.0, "s");
    TProfile2D* h_Q2_RelRes_vs_xy_Sigma = new TProfile2D("Q2_RelRes_vs_xy_Sigma",
        "Q^{2} Rel. Res. binned in (x_{Bj}, y) (Sigma);x_{Bj};y",
        x_bins.size()-1, x_bins.data(), n_y_bins, 0.0, 1.0, "s");

    TH1D* h_RelRes_Q2_EM = new TH1D("Q2_RelRes_EM","electron method;#frac{Q^{2}(Reco)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    TH1D* h_RelRes_Q2_DA = new TH1D("Q2_RelRes_DA","DA method;#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    TH1D* h_RelRes_Q2_Sigma = new TH1D("Q2_RelRes_Sigma","Sigma method;#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",101,-0.15,0.15);
    
    TH2D* h_RelRes_Q2_binned_EM = new TH2D("Q2_RelRes_binned_EM",";Q^{2} [GeV^{2}];#frac{Q^{2}(EM)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_DA = new TH2D("Q2_RelRes_binned_DA",";Q^{2} [GeV^{2}];#frac{Q^{2}(DA)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);
    TH2D* h_RelRes_Q2_binned_Sigma = new TH2D("Q2_RelRes_binned_Sigma",";Q^{2} [GeV^{2}];#frac{Q^{2}(Sigma)-Q^{2}(MC)}{Q^{2}(MC)}",
        n_bins, bin_edges_Q2.data(), n_binned,-0.15,0.15);

    
    TH2D* h_Corr_Q2_EM = new TH2D("Corr_Q2_EM", ";Q^{2}_{MC};Q^{2}_{EM}",
                                  n_bins, bin_edges_Q2.data(),
                                  n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_DA = new TH2D("Corr_Q2_DA", ";Q^{2}_{MC};Q^{2}_{DA}",
                                    n_bins, bin_edges_Q2.data(),
                                    n_bins, bin_edges_Q2.data());
    TH2D* h_Corr_Q2_Sigma = new TH2D("Corr_Q2_Sigma", ";Q^{2}_{MC};Q^{2}_{Sigma}",
                                    n_bins, bin_edges_Q2.data(),
                                    n_bins, bin_edges_Q2.data());

    TH1D* h_Q2_truth    = new TH1D("h_Q2_truth","Q^2;# of events",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_EM       = new TH1D("h_Q2_EM",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_DA       = new TH1D("h_Q2_DA",";Q^{2}",n_bins, bin_edges_Q2.data());
    TH1D* h_Q2_Sigma   = new TH1D("h_Q2_Sigma",";Q^{2}",n_bins, bin_edges_Q2.data());

    // E-pz histograms for MATCHED particles
    TH1D* h_EPz_truth = new TH1D("h_EPz_truth", "MC Truth Sum(E-p_{z}) - Matched Particles Only;#Sigma(E-p_{z}) [GeV];Counts", 50, 0, 25);
    TH1D* h_EPz = new TH1D("h_EPz", "Reco Sum(E-p_{z}) - Matched Particles Only;#Sigma(E-p_{z}) [GeV];Counts", 50, 0, 25);

    // Eta_max histograms (reco and truth)
    TH1D* h_eta_max = new TH1D("h_eta_max", "Maximum Pseudorapidity per Event (Reco);#eta_{max};Counts", 50, -4.0, 6.0);
    TH1D* h_eta_max_truth = new TH1D("h_eta_max_truth", "Maximum Pseudorapidity per Event (Truth);#eta_{max};Counts", 50, -4.0, 6.0);

    // M_X^2 histograms (hadronic invariant mass squared, excluding scattered electron)
    TH1D* h_MX2 = new TH1D("h_MX2", "Hadronic Invariant Mass Squared (Reco);M_{X}^{2} [GeV^{2}];Counts", 100, 0.0, 200.0);
    TH1D* h_MX2_truth = new TH1D("h_MX2_truth", "Hadronic Invariant Mass Squared (Truth);M_{X}^{2} [GeV^{2}];Counts", 100, 0.0, 200.0);

    //---------------------------------------------------------
    // DECLARE OUTPUT HISTOGRAMS (t part)
    //---------------------------------------------------------
    
    // Setup binning for t
    std::vector<Double_t> t_bins_low = GetLogBins(1e-3, 0.5, 20);
    std::vector<Double_t> t_bins_high = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.6};
    std::vector<Double_t> t_bins = t_bins_low;
    for(size_t i = 1; i < t_bins_high.size(); i++){
        t_bins.push_back(t_bins_high[i]);
    }
    
    // Mandelstam t histograms
    TH1D* h_t_MC = new TH1D("t_MC", "Truth Mandelstam t;|t| [GeV^{2}];Counts",
                            t_bins.size()-1, t_bins.data());
    TH1D* h_t_B0 = new TH1D("t_B0", "B0 Reco Mandelstam t;|t| [GeV^{2}];Counts",
                            t_bins.size()-1, t_bins.data());
    TH1D* h_t_RP_histo = new TH1D("t_RP_histo", "RP Reco Mandelstam t;|t| [GeV^{2}];Counts",
                                   t_bins.size()-1, t_bins.data());

    // Differential cross section d(sigma)/dt histograms
    TH1D* h_dsigma_dt_MC = new TH1D("dsigma_dt_MC", "Truth d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                     t_bins.size()-1, t_bins.data());
    TH1D* h_dsigma_dt_B0 = new TH1D("dsigma_dt_B0", "B0 Reco d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                     t_bins.size()-1, t_bins.data());
    TH1D* h_dsigma_dt_RP = new TH1D("dsigma_dt_RP", "RP Reco d#sigma/dt;|t| [GeV^{2}];d#sigma/dt [nb/GeV^{2}]",
                                     t_bins.size()-1, t_bins.data());

    // Angular histograms
    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);

    // Create logarithmic binning for theta_all_TS (0.01 to 10000 mrad)
    // const int n_theta_bins = 100;
    // double theta_bins[n_theta_bins + 1];
    // double theta_min = 0.01;
    // double theta_max = 10000.0;
    // double log_theta_min = TMath::Log10(theta_min);
    // double log_theta_max = TMath::Log10(theta_max);
    // for(int i = 0; i <= n_theta_bins; i++){
    //     theta_bins[i] = TMath::Power(10, log_theta_min + i * (log_theta_max - log_theta_min) / n_theta_bins);
    // }
    TH1D* h_theta_all_TS = new TH1D("theta_all_TS", "All Truth-Seeded Proton Angles;#theta [mrad];Counts", 100, 0.0, 25.0);

    TH1D* h_theta_B0 = new TH1D("theta_B0", "B0 Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    TH1D* h_theta_RP = new TH1D("theta_RP", "RP Proton Scattering Angle;#theta [mrad];Counts", 100, 0.0, 25.0);
    
    // x_L histograms
    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_B0 = new TH1D("xL_B0", "B0 Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);
    TH1D* h_xL_RP = new TH1D("xL_RP", "RP Reco x_{L};x_{L};Counts", 30, 0.75, 1.05);

    // Create logarithmic binning for x_pom (1e-4 to 0.4)
    const int n_xpom_bins = 20;
    double xpom_bins[n_xpom_bins + 1];
    double xpom_min = 1e-4;
    double xpom_max = 0.4;
    double log_min = TMath::Log10(xpom_min);
    double log_max = TMath::Log10(xpom_max);
    for(int i = 0; i <= n_xpom_bins; i++){
        xpom_bins[i] = TMath::Power(10, log_min + i * (log_max - log_min) / n_xpom_bins);
    }
    
    // x_pom histograms (from x_L: x_pom = 1 - x_L)
    TH1D* h_xpom_MC = new TH1D("xpom_MC", "Truth x_{pom} (from x_L);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_B0 = new TH1D("xpom_B0", "B0 Reco x_{pom} (from x_L);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_RP = new TH1D("xpom_RP", "RP Reco x_{pom} (from x_L);x_{pom};Counts", n_xpom_bins, xpom_bins);
    
    // x_pom histograms from definition: x_pom = (M_X^2 + Q^2 + |t|)/(W^2 + Q^2 - m_p^2)
    TH1D* h_xpom_def_MC = new TH1D("xpom_def_MC", "Truth x_{pom} (from definition);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_def_B0 = new TH1D("xpom_def_B0", "B0 Reco x_{pom} (from definition);x_{pom};Counts", n_xpom_bins, xpom_bins);
    TH1D* h_xpom_def_RP = new TH1D("xpom_def_RP", "RP Reco x_{pom} (from definition);x_{pom};Counts", n_xpom_bins, xpom_bins);
    
    // Comparison histograms: x_pom from definition vs from x_L
    TH2D* h_xpom_comp_MC = new TH2D("xpom_comp_MC", "Truth: x_{pom} from x_L vs from definition;x_{pom} (1-x_L);x_{pom} (definition)",
                                     n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_comp_B0 = new TH2D("xpom_comp_B0", "B0 Reco: x_{pom} from x_L vs from definition;x_{pom} (1-x_L);x_{pom} (definition)",
                                     n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_comp_RP = new TH2D("xpom_comp_RP", "RP Reco: x_{pom} from x_L vs from definition;x_{pom} (1-x_L);x_{pom} (definition)",
                                     n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    
    // Resolution histograms for t
    TH1D* h_t_res_B0 = new TH1D("t_res_B0", "B0 t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 100, -2.0, 2.0);
    TH1D* h_t_res_RP = new TH1D("t_res_RP", "RP t Resolution;(|t|_{reco}-|t|_{truth})/|t|_{truth};Counts", 100, -2.0, 2.0);
    
    // Resolution histograms for x_L
    TH1D* h_xL_res_B0 = new TH1D("xL_res_B0", "B0 x_L Resolution;(x_{L,reco}-x_{L,truth})/x_{L,truth};Counts", 100, -2.0, 2.0);
    TH1D* h_xL_res_RP = new TH1D("xL_res_RP", "RP x_L Resolution;(x_{L,reco}-x_{L,truth})/x_{L,truth};Counts", 100, -2.0, 2.0);
    
    // Resolution histograms for x_pom
    TH1D* h_xpom_res_B0 = new TH1D("xpom_res_B0", "B0 x_pom Resolution (from x_L);(x_{pom,reco}-x_{pom,truth})/x_{pom,truth};Counts", 100, -2.0, 2.0);
    TH1D* h_xpom_res_RP = new TH1D("xpom_res_RP", "RP x_pom Resolution (from x_L);(x_{pom,reco}-x_{pom,truth})/x_{pom,truth};Counts", 100, -2.0, 2.0);
    
    // Correlation histograms for t
    TH2D* h_t_corr_B0 = new TH2D("t_corr_B0", "B0 |t| Correlation;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                  t_bins.size()-1, t_bins.data(), t_bins.size()-1, t_bins.data());
    TH2D* h_t_corr_RP = new TH2D("t_corr_RP", "RP |t| Correlation;Truth |t| [GeV^{2}];Reco |t| [GeV^{2}]",
                                  t_bins.size()-1, t_bins.data(), t_bins.size()-1, t_bins.data());
    
    // Correlation histograms for x_L
    TH2D* h_xL_corr_B0 = new TH2D("xL_corr_B0", "B0 x_L Correlation;Truth x_L;Reco x_L", 100, 0.0, 2.0, 100, 0.0, 2.0);
    TH2D* h_xL_corr_RP = new TH2D("xL_corr_RP", "RP x_L Correlation;Truth x_L;Reco x_L", 100, 0.0, 0.5, 100, 0.0, 0.5);
    
    // Correlation histograms for x_pom
    TH2D* h_xpom_corr_B0 = new TH2D("xpom_corr_B0", "B0 x_pom Correlation (from x_L);Truth x_{pom};Reco x_{pom}",
                                     n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);
    TH2D* h_xpom_corr_RP = new TH2D("xpom_corr_RP", "RP x_pom Correlation (from x_L);Truth x_{pom};Reco x_{pom}",
                                     n_xpom_bins, xpom_bins, n_xpom_bins, xpom_bins);

    // Beta histograms (beta = x_Bj / x_pom, using x_pom from definition)
    TH1D* h_beta_MC = new TH1D("beta_MC", "Truth #beta (from x_{pom} def);#beta;Counts", 10, 0.0, 1.0);
    TH1D* h_beta_B0 = new TH1D("beta_B0", "B0 Reco #beta (from x_{pom} def);#beta;Counts", 10, 0.0, 1.0);
    TH1D* h_beta_RP = new TH1D("beta_RP", "RP Reco #beta (from x_{pom} def);#beta;Counts", 10, 0.0, 1.0);

    // Beta resolution histograms
    TH1D* h_beta_res_B0 = new TH1D("beta_res_B0",
                                    "B0 #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts",
                                    100, -1.0, 1.0);
    TH1D* h_beta_res_RP = new TH1D("beta_res_RP",
                                    "RP #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts",
                                    100, -1.0, 1.0);

    // Beta 2D correlation histograms
    TH2D* h_beta_corr_B0 = new TH2D("beta_corr_B0",
                                     "B0 #beta Correlation;Truth #beta;Reco #beta",
                                     10, 0.0, 1.0, 10, 0.0, 1.0);
    TH2D* h_beta_corr_RP = new TH2D("beta_corr_RP",
                                     "RP #beta Correlation;Truth #beta;Reco #beta",
                                     10, 0.0, 1.0, 10, 0.0, 1.0);

    // Triple differential cross section d³σ/(dQ² dβ dx_pom)
    // Binning: Q² (log), β (linear 0-1), x_pom (log)
    const int n_beta_bins = 5;  // 5 bins from 0 to 1
    double beta_bins[n_beta_bins + 1] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    // Use fewer bins for Q² for the 3D histogram (use coarser binning)
    const int n_Q2_3d_bins = 10;
    std::vector<Double_t> Q2_3d_bins = GetLogBins(1.0, 100.0, n_Q2_3d_bins);

    // Use fewer bins for x_pom for the 3D histogram
    const int n_xpom_3d_bins = 8;
    double xpom_3d_bins[n_xpom_3d_bins + 1];
    double xpom_3d_min = 1e-3;
    double xpom_3d_max = 0.1;
    double log_3d_min = TMath::Log10(xpom_3d_min);
    double log_3d_max = TMath::Log10(xpom_3d_max);
    for(int i = 0; i <= n_xpom_3d_bins; i++){
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
    TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};
    TTreeReaderArray<float>  tsre_px_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float>  tsre_py_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float>  tsre_pz_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float> rp_px_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    TTreeReaderArray<float> rp_py_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    TTreeReaderArray<float> rp_pz_array           = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
    TTreeReaderArray<float> rp_mass_array         = {tree_reader, "ForwardRomanPotRecParticles.mass"};
    TTreeReaderArray<int>   rp_pdg_array          = {tree_reader, "ForwardRomanPotRecParticles.PDG"};

    // Set up SetBranchAddress for Q2, x, y, W from InclusiveKinematics
    Float_t electron_Q2_EM, electron_Q2_DA, electron_Q2_Sigma, electron_Q2_truth;
    Float_t electron_x_EM, electron_x_DA, electron_x_Sigma, electron_x_truth;
    Float_t electron_y_EM, electron_y_DA, electron_y_Sigma, electron_y_truth;
    Float_t electron_W_EM, electron_W_DA, electron_W_Sigma, electron_W_truth;

    events->SetBranchAddress("InclusiveKinematicsElectron.Q2", &electron_Q2_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.Q2"      , &electron_Q2_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.Q2"  , &electron_Q2_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.Q2"   , &electron_Q2_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.x", &electron_x_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.x"      , &electron_x_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.x"  , &electron_x_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.x"   , &electron_x_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.y", &electron_y_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.y"      , &electron_y_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.y"  , &electron_y_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.y"   , &electron_y_truth);

    events->SetBranchAddress("InclusiveKinematicsElectron.W", &electron_W_EM);
    events->SetBranchAddress("InclusiveKinematicsDA.W"      , &electron_W_DA);
    events->SetBranchAddress("InclusiveKinematicsSigma.W"  , &electron_W_Sigma);
    events->SetBranchAddress("InclusiveKinematicsTruth.W"   , &electron_W_truth);

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
    for (Long64_t i = 0; i < nentries; i++) {
        // Update the counter every 1000 events
        if (i % 1000 == 0) {
        printf("\rProcessing event %lld of %lld; %.2f percent done.", i, nentries, 100.0*i/nentries);
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

        // Relative resolutions vs truth (guard against zero)
        if (electron_x_truth != 0.0f) {
            h_RelRes_x_EM->Fill(    (electron_x_EM     - electron_x_truth)    / electron_x_truth);
            h_RelRes_x_DA->Fill(    (electron_x_DA     - electron_x_truth)    / electron_x_truth);
            h_RelRes_x_Sigma->Fill((electron_x_Sigma - electron_x_truth)    / electron_x_truth);
            h_RelRes_x_binned_EM->Fill(    electron_x_truth, (electron_x_EM     - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_DA->Fill(    electron_x_truth, (electron_x_DA     - electron_x_truth)/electron_x_truth);
            h_RelRes_x_binned_Sigma->Fill(electron_x_truth, (electron_x_Sigma - electron_x_truth)/electron_x_truth);
        }
        if (electron_y_truth != 0.0f) {
            h_RelRes_y_EM->Fill(    (electron_y_EM     - electron_y_truth)    / electron_y_truth);
            h_RelRes_y_DA->Fill(    (electron_y_DA     - electron_y_truth)    / electron_y_truth);
            h_RelRes_y_Sigma->Fill((electron_y_Sigma - electron_y_truth)    / electron_y_truth);
            h_RelRes_y_binned_EM->Fill(    electron_y_truth, (electron_y_EM     - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_DA->Fill(    electron_y_truth, (electron_y_DA     - electron_y_truth)/electron_y_truth);
            h_RelRes_y_binned_Sigma->Fill(electron_y_truth, (electron_y_Sigma - electron_y_truth)/electron_y_truth);
        }

        // Correlations
        h_Corr_x_EM->Fill(electron_x_truth, electron_x_EM);
        h_Corr_x_DA->Fill(electron_x_truth, electron_x_DA);
        h_Corr_x_Sigma->Fill(electron_x_truth, electron_x_Sigma);

        h_Corr_y_EM->Fill(electron_y_truth, electron_y_EM);
        h_Corr_y_DA->Fill(electron_y_truth, electron_y_DA);
        h_Corr_y_Sigma->Fill(electron_y_truth, electron_y_Sigma);

        h_Q2_EM->Fill(electron_Q2_EM);
        h_Q2_DA->Fill(electron_Q2_DA);
        h_Q2_Sigma->Fill(electron_Q2_Sigma);

        // Also fill the output tree for x and y
        out_x_truth  = electron_x_truth;
        out_x_EM     = electron_x_EM;
        out_x_DA     = electron_x_DA;
        out_x_Sigma = electron_x_Sigma;
        out_y_truth  = electron_y_truth;
        out_y_EM     = electron_y_EM;
        out_y_DA     = electron_y_DA;
        out_y_Sigma = electron_y_Sigma;
        tree->Fill();
        
        // Fill histograms for truth
        h_Q2_truth->Fill(electron_Q2_truth);

        h_RelRes_Q2_EM->Fill((electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_DA->Fill((electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_Sigma->Fill((electron_Q2_Sigma - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_EM->Fill(electron_Q2_truth, (electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_DA->Fill(electron_Q2_truth, (electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
        h_RelRes_Q2_binned_Sigma->Fill(electron_Q2_truth, (electron_Q2_Sigma - electron_Q2_truth)/electron_Q2_truth);

        // Fill Q2 resolution profile histograms - binned in (x_Bj, y)
        if(electron_Q2_truth > 0.0f && electron_x_truth > 0.0f && electron_y_truth > 0.0f) {
            h_Q2_RelRes_vs_xy_EM->Fill(electron_x_truth, electron_y_truth,
                                       (electron_Q2_EM - electron_Q2_truth)/electron_Q2_truth);
            h_Q2_RelRes_vs_xy_DA->Fill(electron_x_truth, electron_y_truth,
                                       (electron_Q2_DA - electron_Q2_truth)/electron_Q2_truth);
            h_Q2_RelRes_vs_xy_Sigma->Fill(electron_x_truth, electron_y_truth,
                                           (electron_Q2_Sigma - electron_Q2_truth)/electron_Q2_truth);
        }

        // Fill correlation histograms
        h_Corr_Q2_EM->Fill(electron_Q2_truth, electron_Q2_EM);
        h_Corr_Q2_DA->Fill(electron_Q2_truth, electron_Q2_DA);
        h_Corr_Q2_Sigma->Fill(electron_Q2_truth, electron_Q2_Sigma);

        // Calculate and fill E-pz histograms using MATCHED particles only
        double sumEPz_truth_matched, sumEPz_reco_matched;
        CalculateSumEPz_Matched(
            re_px_array, re_py_array, re_pz_array, re_energy_array,
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            assoc_rec_id, assoc_sim_id,
            sumEPz_truth_matched, sumEPz_reco_matched
        );
        h_EPz_truth->Fill(sumEPz_truth_matched);
        h_EPz->Fill(sumEPz_reco_matched);

        // Calculate eta_max for reco and truth (using matched particles via associations)
        double eta_max_reco = -999.0;
        double eta_max_truth = -999.0;
        
        // Loop over reconstructed particles
        for(unsigned int j = 0; j < re_energy_array.GetSize(); j++){
            // Calculate reco eta
            P3MVector particle_reco(re_px_array[j], re_py_array[j], re_pz_array[j], 0.0);
            double eta_reco = particle_reco.Eta();
            if(eta_reco > eta_max_reco){
                eta_max_reco = eta_reco;
            }
            
            // Find corresponding MC particle through associations
            int mc_idx = -1;
            for(unsigned int k = 0; k < assoc_rec_id.GetSize(); k++){
                if(assoc_rec_id[k] == j) {
                    mc_idx = assoc_sim_id[k];
                    break;
                }
            }
            
            // If association found, calculate MC eta
            if(mc_idx >= 0 && mc_idx < mc_px_array.GetSize()){
                P3MVector particle_mc(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
                double eta_mc = particle_mc.Eta();
                if(eta_mc > eta_max_truth){
                    eta_max_truth = eta_mc;
                }
            }
        }
        
        // Fill histograms if valid
        if(eta_max_reco > -999.0){
            h_eta_max->Fill(eta_max_reco);
        }
        if(eta_max_truth > -999.0){
            h_eta_max_truth->Fill(eta_max_truth);
        }

        // Calculate M_X^2 (hadronic invariant mass squared, excluding scattered electron)
        // Get scattered electron index
        int scat_e_idx = -1;
        if(electron_scat_index.GetSize() > 0) {
            scat_e_idx = electron_scat_index[0];
        }

        // Sum all reconstructed particles (excluding scattered electron)
        P3EVector total_hadrons_reco(0.0, 0.0, 0.0, 0.0);
        for(unsigned int j = 0; j < re_energy_array.GetSize(); j++){
            // Skip the scattered electron
            if((int)j == scat_e_idx) continue;
            
            P3EVector particle(re_px_array[j], re_py_array[j], re_pz_array[j], re_energy_array[j]);
            total_hadrons_reco += particle;
        }
        double MX2_reco = total_hadrons_reco.M2();
        
        // Calculate M_X^2 for truth (sum matched MC particles, excluding scattered electron)
        P3EVector total_hadrons_truth(0.0, 0.0, 0.0, 0.0);
        for(unsigned int j = 0; j < re_energy_array.GetSize(); j++){
            // Skip the scattered electron
            if((int)j == scat_e_idx) continue;
            
            // Find corresponding MC particle through associations
            int mc_idx = -1;
            for(unsigned int k = 0; k < assoc_rec_id.GetSize(); k++){
                if(assoc_rec_id[k] == j) {
                    mc_idx = assoc_sim_id[k];
                    break;
                }
            }
            
            // If association found, add to truth sum
            if(mc_idx >= 0 && mc_idx < mc_px_array.GetSize()){
                double px_mc = mc_px_array[mc_idx];
                double py_mc = mc_py_array[mc_idx];
                double pz_mc = mc_pz_array[mc_idx];
                double m_mc = mc_mass_array[mc_idx];
                double E_mc = TMath::Sqrt(px_mc*px_mc + py_mc*py_mc + pz_mc*pz_mc + m_mc*m_mc);
                P3EVector particle_mc(px_mc, py_mc, pz_mc, E_mc);
                total_hadrons_truth += particle_mc;
            }
        }
        double MX2_truth = total_hadrons_truth.M2();
        
        // Fill M_X^2 histograms
        h_MX2->Fill(MX2_reco);
        h_MX2_truth->Fill(MX2_truth);

        //=================================================================
        // PROCESS PROTONS (t analysis part)
        //=================================================================
        
        // Process truth protons
        std::vector<P3MVector> truth_protons;
        for(int j = 0; j < mc_px_array.GetSize(); j++){
            if(mc_genStatus_array[j] == 1 && mc_pdg_array[j] == 2212){
                P3MVector p(mc_px_array[j], mc_py_array[j], mc_pz_array[j], mc_mass_array[j]);
                undoAfterburn(p);

                // Calculate x_L using p'z/pz method
                double xL_pz = p.Pz() / beams.p_beam.Pz();
                double xpom_from_xL = 1.0 - xL_pz;

                // Calculate Mandelstam t
                double t_val = TMath::Abs(CalcT(beams.p_beam, p));
                
                // Calculate x_pom from definition: x_pom = (M_X^2 + Q^2 - t)/(W^2 + Q^2 - m_p^2)
                double m_p_sq = fMass_proton * fMass_proton;
                double W2_truth = electron_W_truth * electron_W_truth;
                double xpom_from_def = -999.0;
                double denominator = W2_truth + electron_Q2_truth - m_p_sq;
                if(denominator > 0) {
                    xpom_from_def = (MX2_truth + electron_Q2_truth + t_val) / denominator;
                }

                truth_protons.push_back(p);
                h_theta_MC->Fill(p.Theta() * 1000.0);
                h_t_MC->Fill(t_val);
                h_dsigma_dt_MC->Fill(t_val);
                h_xL_MC->Fill(xL_pz);
                h_xpom_MC->Fill(xpom_from_xL);
                
                if(xpom_from_def > 0 && xpom_from_def < 1.0) {
                    h_xpom_def_MC->Fill(xpom_from_def);
                    h_xpom_comp_MC->Fill(xpom_from_xL, xpom_from_def);

                    // Calculate beta = x_Bj / x_pom (using x_pom from definition)
                    if(electron_x_truth > 0) {
                        double beta_truth = electron_x_truth / xpom_from_def;
                        if(beta_truth > 0 && beta_truth <= 1.0) {
                            h_beta_MC->Fill(beta_truth);

                            // Fill 3D histogram for triple differential cross section
                            if(electron_Q2_truth > 0) {
                                h_d3sigma_MC->Fill(electron_Q2_truth, beta_truth, xpom_from_def);
                            }
                        }
                    }
                }
            }
        }
        
        // Process B0 protons (truth-seeded)
        for(unsigned int j = 0; j < tsassoc_rec_id.GetSize(); j++){
            auto mc_idx = tsassoc_sim_id[j];
            // std::cout << "TS assoc sim ID: " << mc_idx << std::endl;
            
            if(mc_idx >= (unsigned)mc_pdg_array.GetSize() || mc_genStatus_array[mc_idx] != 1 || mc_pdg_array[mc_idx] != 2212)
                continue;
            
            // std::cout << "Processing B0 proton with MC index: " << mc_idx << std::endl;
            P3MVector p_reco(tsre_px_array[j], tsre_py_array[j], tsre_pz_array[j], mc_mass_array[mc_idx]);
            undoAfterburn(p_reco);

            // Fill histogram for ALL truth-seeded protons (before acceptance cut)
            h_theta_all_TS->Fill(p_reco.Theta() * 1000.0);

            // B0 angular acceptance
            // std::cout << "B0 proton Theta: " << p_reco.Theta() << std::endl;
            if(p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;
            // std::cout << "B0 proton accepted with Theta: " << p_reco.Theta() << std::endl;
            
            h_theta_B0->Fill(p_reco.Theta() * 1000.0);
            
            // Get truth for correlation
            P3MVector p_truth(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
            undoAfterburn(p_truth);
            
            double t_reco = CalcT(beams.p_beam, p_reco);
            double t_truth = CalcT(beams.p_beam, p_truth);
            double t_reco_abs = TMath::Abs(t_reco);
            double t_truth_abs = TMath::Abs(t_truth);

            // Calculate x_L using p'z/pz method
            double xL_reco_pz = p_reco.Pz() / beams.p_beam.Pz();
            double xL_truth_pz = p_truth.Pz() / beams.p_beam.Pz();
            double xpom_reco_from_xL = 1.0 - xL_reco_pz;
            double xpom_truth_from_xL = 1.0 - xL_truth_pz;

            // Calculate x_pom from definition for B0
            double m_p_sq = fMass_proton * fMass_proton;
            double W2_EM = electron_W_EM * electron_W_EM;
            double W2_truth = electron_W_truth * electron_W_truth;
            
            double xpom_reco_from_def = -999.0;
            double denominator_reco = W2_EM + electron_Q2_EM - m_p_sq;
            if(denominator_reco > 0) {
                xpom_reco_from_def = (MX2_reco + electron_Q2_EM + t_reco_abs) / denominator_reco;
            }
            
            double xpom_truth_from_def = -999.0;
            double denominator_truth = W2_truth + electron_Q2_truth - m_p_sq;
            if(denominator_truth > 0) {
                xpom_truth_from_def = (MX2_truth + electron_Q2_truth + t_truth_abs) / denominator_truth;
            }

            // Fill histograms
            h_t_B0->Fill(t_reco_abs);
            h_dsigma_dt_B0->Fill(t_reco_abs);
            h_t_corr_B0->Fill(t_truth_abs, t_reco_abs);
            if(t_truth_abs > 1e-6) {
                h_t_res_B0->Fill((t_reco_abs - t_truth_abs) / t_truth_abs);
            }
            
            h_xL_B0->Fill(xL_reco_pz);
            h_xL_corr_B0->Fill(xL_truth_pz, xL_reco_pz);
            if(xL_truth_pz > 1e-6) {
                h_xL_res_B0->Fill((xL_reco_pz - xL_truth_pz) / xL_truth_pz);
            }
            
            h_xpom_B0->Fill(xpom_reco_from_xL);
            h_xpom_corr_B0->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
            if(xpom_truth_from_xL > 1e-6) {
                h_xpom_res_B0->Fill((xpom_reco_from_xL - xpom_truth_from_xL) / xpom_truth_from_xL);
            }
            
            // Fill x_pom from definition histograms
            if(xpom_reco_from_def > 0 && xpom_reco_from_def < 1.0) {
                h_xpom_def_B0->Fill(xpom_reco_from_def);
                h_xpom_comp_B0->Fill(xpom_reco_from_xL, xpom_reco_from_def);

                // Calculate beta = x_Bj / x_pom (using x_pom from definition) for B0
                if(electron_x_EM > 0) {
                    double beta_reco = electron_x_EM / xpom_reco_from_def;
                    if(beta_reco > 0 && beta_reco <= 1.0) {
                        h_beta_B0->Fill(beta_reco);

                        // Fill 3D histogram for triple differential cross section
                        if(electron_Q2_EM > 0) {
                            h_d3sigma_B0->Fill(electron_Q2_EM, beta_reco, xpom_reco_from_def);
                        }

                        // Calculate truth beta if xpom_truth_from_def is valid
                        if(xpom_truth_from_def > 0 && xpom_truth_from_def < 1.0 && electron_x_truth > 0) {
                            double beta_truth = electron_x_truth / xpom_truth_from_def;
                            if(beta_truth > 0 && beta_truth <= 1.0) {
                                h_beta_corr_B0->Fill(beta_truth, beta_reco);
                                h_beta_res_B0->Fill((beta_reco - beta_truth) / beta_truth);
                            }
                        }
                    }
                }
            }
        }
        
        // Process RP protons
        for(int j = 0; j < rp_px_array.GetSize(); j++){
            if(rp_pdg_array[j] != 2212) continue;
            
            P3MVector p_rp(rp_px_array[j], rp_py_array[j], rp_pz_array[j], rp_mass_array[j]);
            
            h_theta_RP->Fill(p_rp.Theta() * 1000.0);
            
            // Find matching truth proton
            P3MVector* best_match = nullptr;
            double best_dr = 0.1;
            for(auto& p_truth : truth_protons){
                double dr = TMath::Sqrt(TMath::Power(p_rp.Theta() - p_truth.Theta(), 2) + 
                                        TMath::Power(p_rp.Phi() - p_truth.Phi(), 2));
                if(dr < best_dr){
                    best_dr = dr;
                    best_match = &p_truth;
                }
            }
            
            if(best_match){
                double t_reco = CalcT(beams.p_beam, p_rp);
                double t_truth = CalcT(beams.p_beam, *best_match);
                double t_reco_abs = TMath::Abs(t_reco);
                double t_truth_abs = TMath::Abs(t_truth);

                // Calculate x_L using p'z/pz method
                double xL_reco_pz = p_rp.Pz() / beams.p_beam.Pz();
                double xL_truth_pz = best_match->Pz() / beams.p_beam.Pz();
                double xpom_reco_from_xL = 1.0 - xL_reco_pz;
                double xpom_truth_from_xL = 1.0 - xL_truth_pz;

                // Calculate x_pom from definition for RP
                double m_p_sq = fMass_proton * fMass_proton;
                double W2_EM = electron_W_EM * electron_W_EM;
                double W2_truth = electron_W_truth * electron_W_truth;
                
                double xpom_reco_from_def = -999.0;
                double denominator_reco = W2_EM + electron_Q2_EM - m_p_sq;
                if(denominator_reco > 0) {
                    xpom_reco_from_def = (MX2_reco + electron_Q2_EM + t_reco_abs) / denominator_reco;
                }
                
                double xpom_truth_from_def = -999.0;
                double denominator_truth = W2_truth + electron_Q2_truth - m_p_sq;
                if(denominator_truth > 0) {
                    xpom_truth_from_def = (MX2_truth + electron_Q2_truth +  t_truth_abs) / denominator_truth;
                }

                // Fill histograms
                h_t_RP_histo->Fill(t_reco_abs);
                h_dsigma_dt_RP->Fill(t_reco_abs);
                h_t_corr_RP->Fill(t_truth_abs, t_reco_abs);
                if(t_truth_abs > 1e-6) {
                    h_t_res_RP->Fill((t_reco_abs - t_truth_abs) / t_truth_abs);
                }
                
                h_xL_RP->Fill(xL_reco_pz);
                h_xL_corr_RP->Fill(xL_truth_pz, xL_reco_pz);
                if(xL_truth_pz > 1e-6) {
                    h_xL_res_RP->Fill((xL_reco_pz - xL_truth_pz) / xL_truth_pz);
                }
                
                h_xpom_RP->Fill(xpom_reco_from_xL);
                h_xpom_corr_RP->Fill(xpom_truth_from_xL, xpom_reco_from_xL);
                if(xpom_truth_from_xL > 1e-6) {
                    h_xpom_res_RP->Fill((xpom_reco_from_xL - xpom_truth_from_xL) / xpom_truth_from_xL);
                }
                
                // Fill x_pom from definition histograms
                if(xpom_reco_from_def > 0 && xpom_reco_from_def < 1.0) {
                    h_xpom_def_RP->Fill(xpom_reco_from_def);
                    h_xpom_comp_RP->Fill(xpom_reco_from_xL, xpom_reco_from_def);

                    // Calculate beta = x_Bj / x_pom (using x_pom from definition) for RP
                    if(electron_x_EM > 0) {
                        double beta_reco = electron_x_EM / xpom_reco_from_def;
                        if(beta_reco > 0 && beta_reco <= 1.0) {
                            h_beta_RP->Fill(beta_reco);

                            // Fill 3D histogram for triple differential cross section
                            if(electron_Q2_EM > 0) {
                                h_d3sigma_RP->Fill(electron_Q2_EM, beta_reco, xpom_reco_from_def);
                            }

                            // Calculate truth beta if xpom_truth_from_def is valid
                            if(xpom_truth_from_def > 0 && xpom_truth_from_def < 1.0 && electron_x_truth > 0) {
                                double beta_truth = electron_x_truth / xpom_truth_from_def;
                                if(beta_truth > 0 && beta_truth <= 1.0) {
                                    h_beta_corr_RP->Fill(beta_truth, beta_reco);
                                    h_beta_res_RP->Fill((beta_reco - beta_truth) / beta_truth);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout<<"\nDone looping over events.\n"<<std::endl;

    // Calculate and scale differential cross sections d(sigma)/dt
    // dsigma/dt = (sigma_total / N_gen) * (N_events_in_bin / delta_t)
    const double sigma_total = 4.22; // nb for 10x100 GeV configuration
    const double N_gen = 100000.0;   // Number of generated events
    const double scale_factor = sigma_total / N_gen;

    std::cout << "Calculating differential cross sections..." << std::endl;
    std::cout << "Total cross section: " << sigma_total << " nb" << std::endl;
    std::cout << "Generated events: " << N_gen << std::endl;

    // Scale MC dsigma/dt histogram
    for(int i = 1; i <= h_dsigma_dt_MC->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_MC->GetBinContent(i);
        double bin_error = h_dsigma_dt_MC->GetBinError(i);
        double bin_width = h_dsigma_dt_MC->GetBinWidth(i);

        double dsigma_dt = (bin_content * scale_factor) / bin_width;
        double dsigma_dt_error = (bin_error * scale_factor) / bin_width;

        h_dsigma_dt_MC->SetBinContent(i, dsigma_dt);
        h_dsigma_dt_MC->SetBinError(i, dsigma_dt_error);
    }

    // Scale B0 dsigma/dt histogram
    for(int i = 1; i <= h_dsigma_dt_B0->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_B0->GetBinContent(i);
        double bin_error = h_dsigma_dt_B0->GetBinError(i);
        double bin_width = h_dsigma_dt_B0->GetBinWidth(i);

        double dsigma_dt = (bin_content * scale_factor) / bin_width;
        double dsigma_dt_error = (bin_error * scale_factor) / bin_width;

        h_dsigma_dt_B0->SetBinContent(i, dsigma_dt);
        h_dsigma_dt_B0->SetBinError(i, dsigma_dt_error);
    }

    // Scale RP dsigma/dt histogram
    for(int i = 1; i <= h_dsigma_dt_RP->GetNbinsX(); i++) {
        double bin_content = h_dsigma_dt_RP->GetBinContent(i);
        double bin_error = h_dsigma_dt_RP->GetBinError(i);
        double bin_width = h_dsigma_dt_RP->GetBinWidth(i);

        double dsigma_dt = (bin_content * scale_factor) / bin_width;
        double dsigma_dt_error = (bin_error * scale_factor) / bin_width;

        h_dsigma_dt_RP->SetBinContent(i, dsigma_dt);
        h_dsigma_dt_RP->SetBinError(i, dsigma_dt_error);
    }

    std::cout << "Differential cross sections calculated!" << std::endl;

    // Scale triple differential cross section d³σ/(dQ² dβ dx_pom) histograms
    std::cout << "Calculating triple differential cross sections..." << std::endl;

    // Scale MC d3sigma histogram
    for(int i = 1; i <= h_d3sigma_MC->GetNbinsX(); i++) {          // Q² bins
        for(int j = 1; j <= h_d3sigma_MC->GetNbinsY(); j++) {      // β bins
            for(int k = 1; k <= h_d3sigma_MC->GetNbinsZ(); k++) {  // x_pom bins
                double bin_content = h_d3sigma_MC->GetBinContent(i, j, k);
                double bin_error = h_d3sigma_MC->GetBinError(i, j, k);

                // Get bin widths for all three dimensions
                double width_Q2 = h_d3sigma_MC->GetXaxis()->GetBinWidth(i);
                double width_beta = h_d3sigma_MC->GetYaxis()->GetBinWidth(j);
                double width_xpom = h_d3sigma_MC->GetZaxis()->GetBinWidth(k);
                double bin_volume = width_Q2 * width_beta * width_xpom;

                // Scale: d³σ/(dQ² dβ dx_pom) = (σ_total / N_gen) * (N_events / ΔQ² Δβ Δx_pom)
                double d3sigma = (bin_content * scale_factor) / bin_volume;
                double d3sigma_error = (bin_error * scale_factor) / bin_volume;

                h_d3sigma_MC->SetBinContent(i, j, k, d3sigma);
                h_d3sigma_MC->SetBinError(i, j, k, d3sigma_error);
            }
        }
    }

    // Scale B0 d3sigma histogram
    for(int i = 1; i <= h_d3sigma_B0->GetNbinsX(); i++) {
        for(int j = 1; j <= h_d3sigma_B0->GetNbinsY(); j++) {
            for(int k = 1; k <= h_d3sigma_B0->GetNbinsZ(); k++) {
                double bin_content = h_d3sigma_B0->GetBinContent(i, j, k);
                double bin_error = h_d3sigma_B0->GetBinError(i, j, k);

                double width_Q2 = h_d3sigma_B0->GetXaxis()->GetBinWidth(i);
                double width_beta = h_d3sigma_B0->GetYaxis()->GetBinWidth(j);
                double width_xpom = h_d3sigma_B0->GetZaxis()->GetBinWidth(k);
                double bin_volume = width_Q2 * width_beta * width_xpom;

                double d3sigma = (bin_content * scale_factor) / bin_volume;
                double d3sigma_error = (bin_error * scale_factor) / bin_volume;

                h_d3sigma_B0->SetBinContent(i, j, k, d3sigma);
                h_d3sigma_B0->SetBinError(i, j, k, d3sigma_error);
            }
        }
    }

    // Scale RP d3sigma histogram
    for(int i = 1; i <= h_d3sigma_RP->GetNbinsX(); i++) {
        for(int j = 1; j <= h_d3sigma_RP->GetNbinsY(); j++) {
            for(int k = 1; k <= h_d3sigma_RP->GetNbinsZ(); k++) {
                double bin_content = h_d3sigma_RP->GetBinContent(i, j, k);
                double bin_error = h_d3sigma_RP->GetBinError(i, j, k);

                double width_Q2 = h_d3sigma_RP->GetXaxis()->GetBinWidth(i);
                double width_beta = h_d3sigma_RP->GetYaxis()->GetBinWidth(j);
                double width_xpom = h_d3sigma_RP->GetZaxis()->GetBinWidth(k);
                double bin_volume = width_Q2 * width_beta * width_xpom;

                double d3sigma = (bin_content * scale_factor) / bin_volume;
                double d3sigma_error = (bin_error * scale_factor) / bin_volume;

                h_d3sigma_RP->SetBinContent(i, j, k, d3sigma);
                h_d3sigma_RP->SetBinError(i, j, k, d3sigma_error);
            }
        }
    }

    std::cout << "Triple differential cross sections calculated!" << std::endl;

    // Write all histograms and TTree to the output file
    outputFile->cd();
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

    // Write Q2 resolution profile histograms (binned in x_Bj, y)
    h_Q2_RelRes_vs_xy_EM->Write();
    h_Q2_RelRes_vs_xy_DA->Write();
    h_Q2_RelRes_vs_xy_Sigma->Write();

    h_EPz_truth->Write();
    h_EPz->Write();
    h_eta_max->Write();
    h_eta_max_truth->Write();
    h_MX2->Write();
    h_MX2_truth->Write();

    // Write t analysis histograms
    h_t_MC->Write();
    h_t_B0->Write();
    h_t_RP_histo->Write();
    h_dsigma_dt_MC->Write();
    h_dsigma_dt_B0->Write();
    h_dsigma_dt_RP->Write();
    h_theta_MC->Write();
    h_theta_all_TS->Write();
    h_theta_B0->Write();
    h_theta_RP->Write();
    h_xL_MC->Write();
    h_xL_B0->Write();
    h_xL_RP->Write();
    h_xpom_MC->Write();
    h_xpom_B0->Write();
    h_xpom_RP->Write();
    h_xpom_def_MC->Write();
    h_xpom_def_B0->Write();
    h_xpom_def_RP->Write();
    h_xpom_comp_MC->Write();
    h_xpom_comp_B0->Write();
    h_xpom_comp_RP->Write();
    h_t_res_B0->Write();
    h_t_res_RP->Write();
    h_xL_res_B0->Write();
    h_xL_res_RP->Write();
    h_xpom_res_B0->Write();
    h_xpom_res_RP->Write();
    h_t_corr_B0->Write();
    h_t_corr_RP->Write();
    h_xL_corr_B0->Write();
    h_xL_corr_RP->Write();
    h_xpom_corr_B0->Write();
    h_xpom_corr_RP->Write();

    // Write beta histograms
    h_beta_MC->Write();
    h_beta_B0->Write();
    h_beta_RP->Write();
    h_beta_res_B0->Write();
    h_beta_res_RP->Write();
    h_beta_corr_B0->Write();
    h_beta_corr_RP->Write();

    // Write triple differential cross section histograms
    h_d3sigma_MC->Write();
    h_d3sigma_B0->Write();
    h_d3sigma_RP->Write();

    outputFile->Close();
    delete events;
    delete outputFile;
    
    std::cout << "\nAnalysis complete! Output saved to DDIS_Combined_output.root" << std::endl;
    return 0;
}
