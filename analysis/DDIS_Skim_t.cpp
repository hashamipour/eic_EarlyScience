#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "Math/Vector4D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"

#include "Utility.hpp"
#include "RecoMethods.hpp"
#include "Math/VectorUtil.h"



using std::cout; using std::endl;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>, ROOT::Math::DefaultCoordinateSystemTag>;
using ROOT::Math::VectorUtil::boost;

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

// Process truth protons from MC
std::vector<P3MVector> ProcessTruthProtons(
    TTreeReaderArray<double>& mc_px, TTreeReaderArray<double>& mc_py, TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass, TTreeReaderArray<int>& mc_gen_status, TTreeReaderArray<int>& mc_pdg,
    const BeamInfo& beams, TH1D* h_theta, TH1D* h_t, TH1D* h_xL, TH1D* h_xpom, TH1D* h_beta, float xBj
) {
    std::vector<P3MVector> truth_protons;
    
    for(int i = 0; i < mc_px.GetSize(); i++){
        if(mc_gen_status[i] == 1 && mc_pdg[i] == 2212){
            P3MVector p(mc_px[i], mc_py[i], mc_pz[i], mc_mass[i]);
            undoAfterburn(p);

            // Calculate x_L - Method 1: |P'|/|P| (momentum magnitude ratio)
            double xL_P = p.P() / beams.p_beam.P();
            // Calculate x_L - Method 2: Pz/E (light-cone momentum fraction)
            double xL_PzE = CalcXL(p);
            // Calculate x_L - Method 3: p'z/pz (longitudinal momentum ratio)
            double xL_pz = p.Pz() / beams.p_beam.Pz();

            // Calculate x_pom = 1 - x_L (using p'z/pz method)
            double xpom = 1.0 - xL_pz;

            // Calculate beta = x_Bj / x_pom (if both are valid)
            if (xpom > 0 && xBj > 0 && xBj <= 1.0) {
                double beta = xBj / xpom;
                if (beta > 0 && beta <= 1.0) {
                    h_beta->Fill(beta);
                }
            }

            // if (xL_pz < 0.99) continue; // x_L cut to avoid the self imposed effects;

            truth_protons.push_back(p);
            h_theta->Fill(p.Theta() * 1000.0);
            h_t->Fill(TMath::Abs(CalcT(beams.p_beam, p)));
            h_xL->Fill(xL_pz);
            h_xpom->Fill(xpom);
        }
    }
    
    return truth_protons;
}

// Process B0 (truth-seeded) protons
void ProcessB0Protons(
    TTreeReaderArray<unsigned int>& tsassoc_rec_id, TTreeReaderArray<unsigned int>& tsassoc_sim_id,
    TTreeReaderArray<float>& tsre_px, TTreeReaderArray<float>& tsre_py, TTreeReaderArray<float>& tsre_pz,
    TTreeReaderArray<double>& mc_px, TTreeReaderArray<double>& mc_py, TTreeReaderArray<double>& mc_pz,
    TTreeReaderArray<double>& mc_mass, TTreeReaderArray<int>& mc_gen_status, TTreeReaderArray<int>& mc_pdg,
    const BeamInfo& beams, MethodHistograms& hist, MethodHistograms& hist_cutFirstBin,
    int& n_matches, double t_first_bin_edge, float xBj, float Q2,
    TH1D* h_beta, TH1D* h_beta_res, TH2D* h_beta_corr, 
    TH2D* h_beta_vs_Q2, TH2D* h_beta_vs_xpom, TH2D* h_beta_vs_t
) {
    for(unsigned int i = 0; i < tsassoc_rec_id.GetSize(); i++){
        auto mc_idx = tsassoc_sim_id[i];
        
        if(mc_idx >= (unsigned)mc_pdg.GetSize() || mc_gen_status[mc_idx] != 1 || mc_pdg[mc_idx] != 2212)
            continue;
            
        P3MVector p_reco(tsre_px[i], tsre_py[i], tsre_pz[i], mc_mass[mc_idx]);
        undoAfterburn(p_reco);
        
        // B0 angular acceptance
        if(p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;
        
        hist.h_theta->Fill(p_reco.Theta() * 1000.0);
        
        // Get truth for correlation
        P3MVector p_truth(mc_px[mc_idx], mc_py[mc_idx], mc_pz[mc_idx], mc_mass[mc_idx]);
        undoAfterburn(p_truth);
        
        double t_reco = CalcT(beams.p_beam, p_reco);
        double t_truth = CalcT(beams.p_beam, p_truth);

        // Calculate x_L - Method 1: |P'|/|P| (momentum magnitude ratio)
        double xL_reco_P = p_reco.P() / beams.p_beam.P();
        double xL_truth_P = p_truth.P() / beams.p_beam.P();
        // Calculate x_L - Method 2: Pz/E (light-cone momentum fraction)
        double xL_reco_PzE = CalcXL(p_reco);
        double xL_truth_PzE = CalcXL(p_truth);
        // Calculate x_L - Method 3: p'z/pz (longitudinal momentum ratio)
        double xL_reco_pz = p_reco.Pz() / beams.p_beam.Pz();
        double xL_truth_pz = p_truth.Pz() / beams.p_beam.Pz();

        // Calculate x_pom = 1 - x_L (using p'z/pz method)
        double xpom_reco = 1.0 - xL_reco_pz;
        double xpom_truth = 1.0 - xL_truth_pz;

        // Calculate beta = x_Bj / x_pom
        double beta_reco = -999.0;
        double beta_truth = -999.0;
        if (xpom_reco > 0 && xBj > 0 && xBj <= 1.0) {
            beta_reco = xBj / xpom_reco;
        }
        if (xpom_truth > 0 && xBj > 0 && xBj <= 1.0) {
            beta_truth = xBj / xpom_truth;
        }

        // if (xL_reco_pz < 0.99) continue; // x_L cut to avoid the self imposed effects;

        // Fill original histograms (all events)
        hist.FillCorrelation(t_truth, t_reco);
        hist.h_xL->Fill(xL_reco_pz);
        hist.h_xL_corr->Fill(xL_truth_pz, xL_reco_pz);
        if(xL_truth_pz > 1e-6) {
            hist.h_xL_res->Fill((xL_reco_pz - xL_truth_pz) / xL_truth_pz);
        }
        hist.h_xpom->Fill(xpom_reco);
        hist.h_xpom_corr->Fill(xpom_truth, xpom_reco);
        if(xpom_truth > 1e-6) {
            hist.h_xpom_res->Fill((xpom_reco - xpom_truth) / xpom_truth);
        }

        // Fill beta histograms
        if (beta_reco > 0 && beta_reco <= 1.0) {
            h_beta->Fill(beta_reco);
            
            // Fill physics correlations
            if (Q2 > 0) {
                h_beta_vs_Q2->Fill(Q2, beta_reco);
            }
            if (xpom_reco > 0) {
                h_beta_vs_xpom->Fill(xpom_reco, beta_reco);
            }
            if (t_reco > 0) {
                h_beta_vs_t->Fill(t_reco, beta_reco);
            }
            
            // Fill resolution and correlation if truth is available
            if (beta_truth > 0 && beta_truth <= 1.0) {
                double beta_res = (beta_reco - beta_truth) / beta_truth;
                h_beta_res->Fill(beta_res);
                h_beta_corr->Fill(beta_truth, beta_reco);
            }
        }

        // Fill cut histograms (skip first |t| bin)
        if(t_reco > t_first_bin_edge) {
            hist_cutFirstBin.FillCorrelation(t_truth, t_reco);
            hist_cutFirstBin.h_xL->Fill(xL_reco_pz);
            hist_cutFirstBin.h_xL_corr->Fill(xL_truth_pz, xL_reco_pz);
            if(xL_truth_pz > 1e-6) {
                hist_cutFirstBin.h_xL_res->Fill((xL_reco_pz - xL_truth_pz) / xL_truth_pz);
            }
            hist_cutFirstBin.h_xpom->Fill(xpom_reco);
            hist_cutFirstBin.h_xpom_corr->Fill(xpom_truth, xpom_reco);
            if(xpom_truth > 1e-6) {
                hist_cutFirstBin.h_xpom_res->Fill((xpom_reco - xpom_truth) / xpom_truth);
            }
        }

        n_matches++;
    }
}

// Process Roman Pot protons
void ProcessRPProtons(
    TTreeReaderArray<float>& rp_px, TTreeReaderArray<float>& rp_py, TTreeReaderArray<float>& rp_pz,
    TTreeReaderArray<float>& rp_mass, TTreeReaderArray<int>& rp_pdg,
    const std::vector<P3MVector>& truth_protons, const BeamInfo& beams,
    MethodHistograms& hist, int& n_matches, float xBj, 
    TH1D* h_beta, TH1D* h_beta_res, TH2D* h_beta_corr
) {
    for(int i = 0; i < rp_px.GetSize(); i++){
        if(rp_pdg[i] != 2212) continue;
        
        P3MVector p_rp(rp_px[i], rp_py[i], rp_pz[i], rp_mass[i]);
        
        // RP angular acceptance
        if(p_rp.Theta() >= 0.005) continue;
        
        hist.h_theta->Fill(p_rp.Theta() * 1000.0);
        
        // Match to truth by angular proximity
        double min_dtheta = 1e9;
        int best_match = -1;
        
        for(size_t j = 0; j < truth_protons.size(); j++){
            if(truth_protons[j].Theta() >= 0.006) continue;
            
            double dtheta = TMath::Abs(p_rp.Theta() - truth_protons[j].Theta());
            if(dtheta < min_dtheta){
                min_dtheta = dtheta;
                best_match = j;
            }
        }
        
        if(best_match >= 0 && min_dtheta < 0.001){
            double t_reco = CalcT(beams.p_beam, p_rp);
            double t_truth = CalcT(beams.p_beam, truth_protons[best_match]);

            // Calculate x_L - Method 1: |P'|/|P| (momentum magnitude ratio)
            double xL_reco_P = p_rp.P() / beams.p_beam.P();
            double xL_truth_P = truth_protons[best_match].P() / beams.p_beam.P();
            // Calculate x_L - Method 2: Pz/E (light-cone momentum fraction)
            double xL_reco_PzE = CalcXL(p_rp);
            double xL_truth_PzE = CalcXL(truth_protons[best_match]);
            // Calculate x_L - Method 3: p'z/pz (longitudinal momentum ratio)
            double xL_reco_pz = p_rp.Pz() / beams.p_beam.Pz();
            double xL_truth_pz = truth_protons[best_match].Pz() / beams.p_beam.Pz();

            // Calculate x_pom = 1 - x_L (using p'z/pz method)
            double xpom_reco = 1.0 - xL_reco_pz;
            double xpom_truth = 1.0 - xL_truth_pz;

            // Calculate beta = x_Bj / x_pom
            double beta_reco = -999.0;
            double beta_truth = -999.0;
            if (xpom_reco > 0 && xBj > 0 && xBj <= 1.0) {
                beta_reco = xBj / xpom_reco;
            }
            if (xpom_truth > 0 && xBj > 0 && xBj <= 1.0) {
                beta_truth = xBj / xpom_truth;
            }

            // if (xL_reco_pz < 0.99) continue; // x_L cut to avoid the self imposed effects;

            hist.FillCorrelation(t_truth, t_reco);
            hist.h_xL->Fill(xL_reco_pz);
            hist.h_xL_corr->Fill(xL_truth_pz, xL_reco_pz);
            if(xL_truth_pz > 1e-6) {
                hist.h_xL_res->Fill((xL_reco_pz - xL_truth_pz) / xL_truth_pz);
            }
            hist.h_xpom->Fill(xpom_reco);
            hist.h_xpom_corr->Fill(xpom_truth, xpom_reco);
            if(xpom_truth > 1e-6) {
                hist.h_xpom_res->Fill((xpom_reco - xpom_truth) / xpom_truth);
            }

            // Fill beta histograms
            if (beta_reco > 0 && beta_reco <= 1.0) {
                h_beta->Fill(beta_reco);
                
                // Fill resolution and correlation if truth is available
                if (beta_truth > 0 && beta_truth <= 1.0) {
                    double beta_res = (beta_reco - beta_truth) / beta_truth;
                    h_beta_res->Fill(beta_res);
                    h_beta_corr->Fill(beta_truth, beta_reco);
                }
            }

            n_matches++;
        }
    }
}

void analyzeProtonsMandelstamT(TString fileList){
    cout << "----------------------------" << endl;
    cout << "  Proton Mandelstam t Analysis" << endl;
    cout << "----------------------------" << endl;
    cout << "\nInput filelist: " << fileList << endl;

    BeamInfo beams;
    
    // Setup input chain
    TChain* events = new TChain("events");
    Int_t nFiles = 0;
    std::vector<std::string> addedFiles;
    
    std::ifstream fileListStream(fileList.Data());
    std::string fileName;
    while(getline(fileListStream, fileName)){
        events->Add((TString)fileName);
        nFiles++;
        addedFiles.push_back(fileName);
    }
    
    cout << "\nNo. of files: " << nFiles << "; no. of events: " << events->GetEntries() << endl;
    
    // Setup binning
    // std::vector<Double_t> t_bins = GetLinBins(0.0, 1.6, 10);
    std::vector<Double_t> t_bins_low = GetLogBins(1e-3, 0.5, 20);
    std::vector<Double_t> t_bins_high = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.6};

    // Combine bins (skip first element of t_bins_high to avoid duplicate 0.5)
    std::vector<Double_t> t_bins = t_bins_low;
    for(size_t i = 1; i < t_bins_high.size(); i++){
        t_bins.push_back(t_bins_high[i]);
    }

    // cout << "bin edges for t: ";
    // for(const auto& edge : t_bins) cout << edge << " ";
    // cout << endl;
    
    // Create histograms for truth
    TH1D* h_t_MC = new TH1D("t_MC", "Truth Mandelstam t;|t| [GeV^{2}];Counts",
                            t_bins.size()-1, t_bins.data());
    TH1D* h_theta_MC = new TH1D("theta_MC", "MC Proton Scattering Angle;#theta [mrad];Counts",
                                100, 0.0, 25.0);
    TH1D* h_xL_MC = new TH1D("xL_MC", "Truth x_{L};x_{L};Counts", 30, 0.75, 1.05);

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
    TH1D* h_xpom_MC = new TH1D("xpom_MC", "Truth x_{pom};x_{pom};Counts", n_xpom_bins, xpom_bins);

    // Create beta histograms (beta = x_Bj / x_pom)
    TH1D* h_beta_MC = new TH1D("beta_MC", "Truth #beta;#beta;Counts", 100, 0.0, 1.0);
    TH1D* h_beta_B0 = new TH1D("beta_B0", "B0 Reco #beta;#beta;Counts", 100, 0.0, 1.0);
    TH1D* h_beta_RP = new TH1D("beta_RP", "RP Reco #beta;#beta;Counts", 100, 0.0, 1.0);
    
    // Beta resolution histograms
    TH1D* h_beta_res_B0 = new TH1D("beta_res_B0", 
                                    "B0 #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts",
                                    100, -2.0, 2.0);
    TH1D* h_beta_res_RP = new TH1D("beta_res_RP", 
                                    "RP #beta Resolution;(#beta_{reco}-#beta_{truth})/#beta_{truth};Counts",
                                    100, -2.0, 2.0);
    
    // Beta 2D correlation histograms
    TH2D* h_beta_corr_B0 = new TH2D("beta_corr_B0", 
                                     "B0 #beta Correlation;Truth #beta;Reco #beta",
                                     100, 0.0, 1.0, 100, 0.0, 1.0);
    TH2D* h_beta_corr_RP = new TH2D("beta_corr_RP", 
                                     "RP #beta Correlation;Truth #beta;Reco #beta",
                                     100, 0.0, 1.0, 100, 0.0, 1.0);
    
    // Physics correlation: beta vs Q2, xpom, t
    TH2D* h_beta_vs_Q2 = new TH2D("beta_vs_Q2", 
                                   "#beta vs Q^{2};Q^{2} [GeV^{2}];#beta",
                                   50, 0.1, 100.0, 50, 0.0, 1.0);
    TH2D* h_beta_vs_xpom = new TH2D("beta_vs_xpom", 
                                     "#beta vs x_{pom};x_{pom};#beta",
                                     50, 1e-4, 0.4, 50, 0.0, 1.0);
    TH2D* h_beta_vs_t = new TH2D("beta_vs_t", 
                                  "#beta vs |t|;|t| [GeV^{2}];#beta",
                                  50, 0.0, 2.0, 50, 0.0, 1.0);

    // Create histograms for each method
    MethodHistograms hist_B0("B0", t_bins);
    MethodHistograms hist_B0_cutFirstBin("B0_cutFirstBin", t_bins);  // Duplicate with first |t| bin cut
    MethodHistograms hist_RP("RP", t_bins);
    
    // Setup tree reader
    TTreeReader tree_reader(events);
    
    TTreeReaderArray<double> mc_px_array(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py_array(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz_array(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> mc_mass_array(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> mc_genStatus_array(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> mc_pdg_array(tree_reader, "MCParticles.PDG");
    
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
    
    // Read x_Bj from InclusiveKinematicsElectron (for beta calculation)
    TTreeReaderArray<float> xBj_Electron(tree_reader, "InclusiveKinematicsElectron.x");
    TTreeReaderArray<float> Q2_Electron(tree_reader, "InclusiveKinematicsElectron.Q2");
    
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
        cout << "Using beam energies from filename tag (" << eBeamGeV << "x"
             << pBeamGeV << " GeV), first match: " << firstMatchedBeamFile << endl;
    }

    if (!usingFilenameBeams) {
        if (!mismatchBeamFile.empty()) {
            std::cerr << "WARNING: Inconsistent beam-energy tags in file names (first match: "
                      << firstMatchedBeamFile << ", mismatch: " << mismatchBeamFile
                      << "). Falling back to beam-particle scan." << std::endl;
        } else {
            cout << "Beam-energy tag not found in file names; falling back to beam-particle scan." << endl;
        }

        // Find beam particles
        cout << "Finding beam particles..." << endl;
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

        cout << "Found beam energies " << beams.e_beam.E() << "x" << beams.p_beam.E() << " GeV" << endl;
    }
    
    undoAfterburnAndCalc(beams.p_beam, beams.e_beam);
    
    tree_reader.Restart();
    
    // Event loop
    int n_truth_protons = 0, n_b0_matches = 0, n_rp_matches = 0;
    
    while(tree_reader.Next()){
        // Get x_Bj and Q2 from InclusiveKinematicsElectron
        float xBj = -999.0;
        float Q2 = -999.0;
        if (xBj_Electron.GetSize() > 0) {
            xBj = xBj_Electron[0];
        }
        if (Q2_Electron.GetSize() > 0) {
            Q2 = Q2_Electron[0];
        }
        
        // Process truth protons
        auto truth_protons = ProcessTruthProtons(
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array, beams, h_theta_MC, h_t_MC, h_xL_MC, h_xpom_MC, h_beta_MC, xBj
        );
        n_truth_protons += truth_protons.size();
        
        // Process B0 protons
        ProcessB0Protons(
            tsassoc_rec_id, tsassoc_sim_id,
            tsre_px_array, tsre_py_array, tsre_pz_array,
            mc_px_array, mc_py_array, mc_pz_array, mc_mass_array,
            mc_genStatus_array, mc_pdg_array,
            beams, hist_B0, hist_B0_cutFirstBin, n_b0_matches, t_bins[1], xBj, Q2,
            h_beta_B0, h_beta_res_B0, h_beta_corr_B0,
            h_beta_vs_Q2, h_beta_vs_xpom, h_beta_vs_t
        );
        
        // Process RP protons
        ProcessRPProtons(
            rp_px_array, rp_py_array, rp_pz_array, rp_mass_array, rp_pdg_array,
            truth_protons, beams, hist_RP, n_rp_matches, xBj,
            h_beta_RP, h_beta_res_RP, h_beta_corr_RP
        );
    }
    
    // Print results
    cout << "\n=== ANALYSIS RESULTS ===" << endl;
    cout << "Truth protons found: " << n_truth_protons << endl;
    cout << "B0 matched protons: " << n_b0_matches << endl;
    cout << "RP matched protons: " << n_rp_matches << endl;
    if(n_truth_protons > 0){
        cout << "B0 matching efficiency: " << (double)n_b0_matches/n_truth_protons*100 << "%" << endl;
        cout << "RP matching efficiency: " << (double)n_rp_matches/n_truth_protons*100 << "%" << endl;
    }
    
    // Save to file
    gStyle->SetOptStat(0);
    
    TFile* outfile = new TFile("proton_mandelstam_analysis.root", "RECREATE");
    h_t_MC->Write();
    h_theta_MC->Write();
    h_xL_MC->Write();
    h_xpom_MC->Write();
    hist_B0.Write();
    hist_B0_cutFirstBin.Write();  // B0 histograms with first |t| bin cut
    hist_RP.Write();
    
    // Write beta histograms
    h_beta_MC->Write();
    h_beta_B0->Write();
    h_beta_RP->Write();
    h_beta_res_B0->Write();
    h_beta_res_RP->Write();
    h_beta_corr_B0->Write();
    h_beta_corr_RP->Write();
    h_beta_vs_Q2->Write();
    h_beta_vs_xpom->Write();
    h_beta_vs_t->Write();
    
    outfile->Close();
    
    cout << "\nAnalysis complete! Output saved to proton_mandelstam_analysis.root" << endl;
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cerr << "Usage: skim_t <filelist.txt>\n";
        return 1;
    }
    
    analyzeProtonsMandelstamT(argv[1]);
    return 0;
}
