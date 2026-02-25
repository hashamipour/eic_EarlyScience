// DDIS binning input skimmer: computes Q2, x_pom, beta (truth/reco) for binning workflow
// g++ DDIS_Skim_BinningInputs.cpp -o DDIS_Skim_BinningInputs $(root-config --cflags --glibs)
// ./DDIS_Skim_BinningInputs filelist.txt [output.root]

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
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
using P3MVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

const Float_t fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

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

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fileList.txt> [output.root]" << std::endl;
        return 1;
    }

    TString fileList = argv[1];
    TString outputName = (argc > 2) ? argv[2] : "DDIS_BinningInputs.root";

    std::ifstream fileListStream;
    fileListStream.open(fileList);
    std::string fileName;

    TChain* events = new TChain("events");
    Int_t nFiles{0};
    std::vector<std::string> addedFiles;
    while(getline(fileListStream, fileName)){
        TString tmp = fileName;
        const bool is_remote = tmp.BeginsWith("root://") || tmp.BeginsWith("xroot://")
                               || tmp.BeginsWith("http://") || tmp.BeginsWith("https://");
        if (!is_remote && !std::filesystem::exists(tmp.Data())) {
            std::cerr << "Error: File does not exist: " << fileName << std::endl;
            continue;
        }
        TFile* inputRootFile = TFile::Open(tmp);
        if(!inputRootFile || inputRootFile->IsZombie()) {
            std::cerr << "Error: Cannot open file: " << fileName << std::endl;
            if(inputRootFile) inputRootFile->Close();
            continue;
        }
        events->Add(tmp);
        inputRootFile->Close();
        nFiles++;
        addedFiles.push_back(fileName);
    }
    std::cout << "No. of files: " << nFiles << "; no. of events: " << events->GetEntries() << std::endl;

    TFile* outputFile = new TFile(outputName, "RECREATE");
    TTree* outTree = new TTree("binning", "Binning inputs (truth/reco)");

    int method = -1; // 0 = B0, 1 = RP
    float Q2_truth = -999.0f, Q2_reco = -999.0f;
    float x_truth = -999.0f, x_reco = -999.0f;
    float y_truth = -999.0f, y_reco = -999.0f;
    float W_truth = -999.0f, W_reco = -999.0f;
    float MX2_truth = -999.0f, MX2_reco = -999.0f;
    float xpom_truth = -999.0f, xpom_reco = -999.0f;
    float beta_truth = -999.0f, beta_reco = -999.0f;
    float t_truth = -999.0f, t_reco = -999.0f;

    outTree->Branch("method", &method, "method/I");
    outTree->Branch("Q2_truth", &Q2_truth, "Q2_truth/F");
    outTree->Branch("Q2_reco", &Q2_reco, "Q2_reco/F");
    outTree->Branch("x_truth", &x_truth, "x_truth/F");
    outTree->Branch("x_reco", &x_reco, "x_reco/F");
    outTree->Branch("y_truth", &y_truth, "y_truth/F");
    outTree->Branch("y_reco", &y_reco, "y_reco/F");
    outTree->Branch("W_truth", &W_truth, "W_truth/F");
    outTree->Branch("W_reco", &W_reco, "W_reco/F");
    outTree->Branch("MX2_truth", &MX2_truth, "MX2_truth/F");
    outTree->Branch("MX2_reco", &MX2_reco, "MX2_reco/F");
    outTree->Branch("xpom_truth", &xpom_truth, "xpom_truth/F");
    outTree->Branch("xpom_reco", &xpom_reco, "xpom_reco/F");
    outTree->Branch("beta_truth", &beta_truth, "beta_truth/F");
    outTree->Branch("beta_reco", &beta_reco, "beta_reco/F");
    outTree->Branch("t_truth", &t_truth, "t_truth/F");
    outTree->Branch("t_reco", &t_reco, "t_reco/F");

    TTreeReader tree_reader(events);
    TTreeReaderArray<double>  mc_px_array         = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<double>  mc_py_array         = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<double>  mc_pz_array         = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double>  mc_mass_array       = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int>     mc_genStatus_array  = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int>     mc_pdg_array        = {tree_reader, "MCParticles.PDG"};
    TTreeReaderArray<unsigned int> assoc_rec_id   = {tree_reader, "ReconstructedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> assoc_sim_id   = {tree_reader, "ReconstructedParticleAssociations.simID"};
    TTreeReaderArray<float>  re_px_array          = {tree_reader, "ReconstructedParticles.momentum.x"};
    TTreeReaderArray<float>  re_py_array          = {tree_reader, "ReconstructedParticles.momentum.y"};
    TTreeReaderArray<float>  re_pz_array          = {tree_reader, "ReconstructedParticles.momentum.z"};
    TTreeReaderArray<float>  re_energy_array      = {tree_reader, "ReconstructedParticles.energy"};
    TTreeReaderArray<int>    electron_scat_index  = {tree_reader, "ScatteredElectronsTruth_objIdx.index"};
    TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};
    TTreeReaderArray<float>  tsre_px_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float>  tsre_py_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float>  tsre_pz_array        = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float>  rp_px_array          = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    TTreeReaderArray<float>  rp_py_array          = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    TTreeReaderArray<float>  rp_pz_array          = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
    TTreeReaderArray<float>  rp_mass_array        = {tree_reader, "ForwardRomanPotRecParticles.mass"};
    TTreeReaderArray<int>    rp_pdg_array         = {tree_reader, "ForwardRomanPotRecParticles.PDG"};

    TTreeReaderArray<float> kin_Q2_EM(tree_reader, "InclusiveKinematicsElectron.Q2");
    TTreeReaderArray<float> kin_Q2_truth(tree_reader, "InclusiveKinematicsTruth.Q2");
    TTreeReaderArray<float> kin_x_EM(tree_reader, "InclusiveKinematicsElectron.x");
    TTreeReaderArray<float> kin_x_truth(tree_reader, "InclusiveKinematicsTruth.x");
    TTreeReaderArray<float> kin_y_EM(tree_reader, "InclusiveKinematicsElectron.y");
    TTreeReaderArray<float> kin_y_truth(tree_reader, "InclusiveKinematicsTruth.y");
    TTreeReaderArray<float> kin_W_EM(tree_reader, "InclusiveKinematicsElectron.W");
    TTreeReaderArray<float> kin_W_truth(tree_reader, "InclusiveKinematicsTruth.W");

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

    Long64_t nentries = events->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        if (i % 1000 == 0) {
            printf("\rProcessing event %lld of %lld; %.2f percent done.", i, nentries, 100.0*i/nentries);
            fflush(stdout);
        }
        tree_reader.Next();
        events->GetEntry(i);

        const float electron_Q2_EM = (kin_Q2_EM.GetSize() > 0) ? kin_Q2_EM[0] : -999.0f;
        const float electron_Q2_truth = (kin_Q2_truth.GetSize() > 0) ? kin_Q2_truth[0] : -999.0f;
        const float electron_x_EM = (kin_x_EM.GetSize() > 0) ? kin_x_EM[0] : -999.0f;
        const float electron_x_truth = (kin_x_truth.GetSize() > 0) ? kin_x_truth[0] : -999.0f;
        const float electron_y_EM = (kin_y_EM.GetSize() > 0) ? kin_y_EM[0] : -999.0f;
        const float electron_y_truth = (kin_y_truth.GetSize() > 0) ? kin_y_truth[0] : -999.0f;
        const float electron_W_EM = (kin_W_EM.GetSize() > 0) ? kin_W_EM[0] : -999.0f;
        const float electron_W_truth = (kin_W_truth.GetSize() > 0) ? kin_W_truth[0] : -999.0f;

        // Scattered electron index
        int scat_e_idx = -1;
        if(electron_scat_index.GetSize() > 0) {
            scat_e_idx = electron_scat_index[0];
        }

        // Reco hadronic system (exclude scattered electron)
        P3EVector total_hadrons_reco(0.0, 0.0, 0.0, 0.0);
        for(unsigned int j = 0; j < re_energy_array.GetSize(); j++){
            if((int)j == scat_e_idx) continue;
            P3EVector particle(re_px_array[j], re_py_array[j], re_pz_array[j], re_energy_array[j]);
            total_hadrons_reco += particle;
        }
        MX2_reco = total_hadrons_reco.M2();

        // Truth hadronic system (sum matched MC particles, excluding scattered electron)
        P3EVector total_hadrons_truth(0.0, 0.0, 0.0, 0.0);
        for(unsigned int j = 0; j < re_energy_array.GetSize(); j++){
            if((int)j == scat_e_idx) continue;

            int mc_idx = -1;
            for(unsigned int k = 0; k < assoc_rec_id.GetSize(); k++){
                if(assoc_rec_id[k] == j) {
                    mc_idx = assoc_sim_id[k];
                    break;
                }
            }
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
        MX2_truth = total_hadrons_truth.M2();

        // Build truth proton list for RP matching
        std::vector<P3MVector> truth_protons;
        for(int j = 0; j < mc_px_array.GetSize(); j++){
            if(mc_genStatus_array[j] == 1 && mc_pdg_array[j] == 2212){
                P3MVector p(mc_px_array[j], mc_py_array[j], mc_pz_array[j], mc_mass_array[j]);
                undoAfterburn(p);
                truth_protons.push_back(p);
            }
        }

        // Cache event-level kinematics
        Q2_truth = electron_Q2_truth;
        Q2_reco = electron_Q2_EM;
        x_truth = electron_x_truth;
        x_reco = electron_x_EM;
        y_truth = electron_y_truth;
        y_reco = electron_y_EM;
        W_truth = electron_W_truth;
        W_reco = electron_W_EM;

        // Process B0 protons (truth-seeded)
        for(unsigned int j = 0; j < tsassoc_rec_id.GetSize(); j++){
            auto mc_idx = tsassoc_sim_id[j];
            if(mc_idx >= (unsigned)mc_pdg_array.GetSize() || mc_genStatus_array[mc_idx] != 1 || mc_pdg_array[mc_idx] != 2212)
                continue;

            P3MVector p_reco(tsre_px_array[j], tsre_py_array[j], tsre_pz_array[j], mc_mass_array[mc_idx]);
            undoAfterburn(p_reco);

            // B0 angular acceptance
            if(p_reco.Theta() <= 0.0055 || p_reco.Theta() >= 0.02) continue;

            P3MVector p_truth(mc_px_array[mc_idx], mc_py_array[mc_idx], mc_pz_array[mc_idx], mc_mass_array[mc_idx]);
            undoAfterburn(p_truth);

            double t_reco_abs = TMath::Abs(CalcT(beams.p_beam, p_reco));
            double t_truth_abs = TMath::Abs(CalcT(beams.p_beam, p_truth));

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

            if(xpom_reco_from_def > 0 && xpom_reco_from_def < 1.0 && xpom_truth_from_def > 0 && xpom_truth_from_def < 1.0) {
                double beta_reco_val = electron_x_EM / xpom_reco_from_def;
                double beta_truth_val = electron_x_truth / xpom_truth_from_def;
                if(beta_reco_val > 0 && beta_reco_val <= 1.0 && beta_truth_val > 0 && beta_truth_val <= 1.0) {
                    method = 0;
                    xpom_reco = xpom_reco_from_def;
                    xpom_truth = xpom_truth_from_def;
                    beta_reco = beta_reco_val;
                    beta_truth = beta_truth_val;
                    t_reco = t_reco_abs;
                    t_truth = t_truth_abs;
                    outTree->Fill();
                }
            }
        }

        // Process RP protons
        for(int j = 0; j < rp_px_array.GetSize(); j++){
            if(rp_pdg_array[j] != 2212) continue;

            P3MVector p_rp(rp_px_array[j], rp_py_array[j], rp_pz_array[j], rp_mass_array[j]);

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
                double t_reco_abs = TMath::Abs(CalcT(beams.p_beam, p_rp));
                double t_truth_abs = TMath::Abs(CalcT(beams.p_beam, *best_match));

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

                if(xpom_reco_from_def > 0 && xpom_reco_from_def < 1.0 && xpom_truth_from_def > 0 && xpom_truth_from_def < 1.0) {
                    double beta_reco_val = electron_x_EM / xpom_reco_from_def;
                    double beta_truth_val = electron_x_truth / xpom_truth_from_def;
                    if(beta_reco_val > 0 && beta_reco_val <= 1.0 && beta_truth_val > 0 && beta_truth_val <= 1.0) {
                        method = 1;
                        xpom_reco = xpom_reco_from_def;
                        xpom_truth = xpom_truth_from_def;
                        beta_reco = beta_reco_val;
                        beta_truth = beta_truth_val;
                        t_reco = t_reco_abs;
                        t_truth = t_truth_abs;
                        outTree->Fill();
                    }
                }
            }
        }
    }

    std::cout << "\nDone looping over events." << std::endl;
    outputFile->Write();
    outputFile->Close();

    return 0;
}
