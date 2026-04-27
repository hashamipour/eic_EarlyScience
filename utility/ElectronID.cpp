#include "ElectronID.hh"

#include "edm4hep/utils/vector_utils.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"

#include <cmath>
#include <iostream>

ElectronID::ElectronID() {

    mEe = 10.;
    mEh = 100.;
    std::cout << "!!! ElectronID: beam energies not specified; defaulting to 10x100 GeV !!!" << std::endl;

    mEoP_min = 0.9;
    mEoP_max = 1.2;

    mDeltaH_min = 0.0;
    mDeltaH_max = 2.0 * mEe;

    mIsoR = 0.4;
    mIsoE = 0.9;

    rcpart_sum_cluster_E  = 0.0;
    rcpart_lead_cluster_E = 0.0;
    rcpart_others_cone_E  = 0.0;
    rcpart_deltaH         = 0.0;
}

ElectronID::ElectronID(double Ee, double Eh) {

    mEe = Ee;
    mEh = Eh;

    mEoP_min = 0.9;
    mEoP_max = 1.2;

    mDeltaH_min = 0.0;
    mDeltaH_max = 2.0 * mEe;

    mIsoR = 0.4;
    mIsoE = 0.9;

    rcpart_sum_cluster_E  = 0.0;
    rcpart_lead_cluster_E = 0.0;
    rcpart_others_cone_E  = 0.0;
    rcpart_deltaH         = 0.0;
}

ElectronID::~ElectronID() {
}


void ElectronID::SetEvent(const podio::Frame* event) {
    mEvent = event;
}


edm4eic::ReconstructedParticleCollection ElectronID::FindScatteredElectron() {

    auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

    edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
    scatteredElectronCandidates.setSubsetCollection();

    for (const auto& reconPart : rcparts) {

        // Negative curvature (track charge <0)
        if (reconPart.getCharge() >= 0) continue;

        // Must have at least one track and one cluster
        if (reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

        CalculateParticleValues(reconPart, rcparts);

        const double p = edm4hep::utils::magnitude(reconPart.getMomentum());
        if (!(p > 0.0)) continue;
        const double recon_EoP = rcpart_sum_cluster_E / p;

        // Paper A.1: E_clus / Sigma(E within dR<0.4) > 0.9. Sigma includes the
        // matched cluster itself, so the denominator is (own cluster energy) +
        // (other particles' cluster energy in the cone).
        const double cone_total_E = rcpart_sum_cluster_E + rcpart_others_cone_E;
        const double recon_isoE   = (cone_total_E > 0.0) ? rcpart_sum_cluster_E / cone_total_E : 0.0;

        if (recon_EoP < mEoP_min || recon_EoP > mEoP_max) continue;
        if (recon_isoE < mIsoE) continue;

        scatteredElectronCandidates.push_back(reconPart);
    }

    return scatteredElectronCandidates;
}

edm4hep::MCParticleCollection ElectronID::GetMCElectron() {

    edm4hep::MCParticleCollection meMC;
    meMC.setSubsetCollection();

    auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

    for (const auto& mcp : mcparts) {
        if (mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 1) {
            meMC.push_back(mcp);
            break; // first status==1 electron
        }
    }

    return meMC;
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

    edm4hep::MCParticleCollection meMC = GetMCElectron();
    edm4eic::ReconstructedParticleCollection meRecon;
    meRecon.setSubsetCollection();

    if (meMC.size() == 0) return meRecon;

    auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

    for (const auto& assoc : RecoMC) {
        if (assoc.getSim() == meMC[0]) {
            meRecon.push_back(assoc.getRec());
            break;
        }
    }

    return meRecon;
}


void ElectronID::CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
                                         const edm4eic::ReconstructedParticleCollection& rcparts) {

    rcpart_sum_cluster_E  = 0.0;
    rcpart_lead_cluster_E = 0.0;
    rcpart_others_cone_E  = 0.0;

    const edm4eic::Cluster* lead_cluster = nullptr;

    for (const auto& cluster : rcp.getClusters()) {
        rcpart_sum_cluster_E += cluster.getEnergy();
        if (cluster.getEnergy() > rcpart_lead_cluster_E) {
            lead_cluster = &cluster;
            rcpart_lead_cluster_E = cluster.getEnergy();
        }
    }

    if (!lead_cluster) return;

    const auto& lead_pos = lead_cluster->getPosition();
    const double lead_eta = edm4hep::utils::eta(lead_pos);
    const double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);

    for (const auto& other_rcp : rcparts) {
        if (other_rcp == rcp) continue;

        for (const auto& other_cluster : other_rcp.getClusters()) {

            const auto& other_pos = other_cluster.getPosition();
            const double other_eta = edm4hep::utils::eta(other_pos);
            const double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

            double d_eta = other_eta - lead_eta;
            double d_phi = other_phi - lead_phi;
            if (d_phi > M_PI)  d_phi -= 2 * M_PI;
            if (d_phi < -M_PI) d_phi += 2 * M_PI;

            const double dR = std::sqrt(d_eta * d_eta + d_phi * d_phi);
            if (dR < mIsoR) {
                rcpart_others_cone_E += other_cluster.getEnergy();
            }
        }
    }
}


edm4eic::ReconstructedParticle ElectronID::SelectHighestPT(const edm4eic::ReconstructedParticleCollection& ecandidates) {

    edm4eic::ReconstructedParticle erec;
    double max_pT = -1.0;

    for (const auto& ecand : ecandidates) {
        const double e_pT = edm4hep::utils::magnitudeTransverse(ecand.getMomentum());
        if (e_pT > max_pT) {
            erec = ecand;
            max_pT = e_pT;
        }
    }

    return erec;
}


double ElectronID::GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp) {

    double sum_cluster_E = 0.0;
    for (const auto& cluster : rcp.getClusters()) {
        sum_cluster_E += cluster.getEnergy();
    }
    return sum_cluster_E;
}


void ElectronID::ComputeEventDeltaH() {

    rcpart_deltaH = 0.0;
    if (!mEvent) return;

    const auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
    for (const auto& p : rcparts) {
        rcpart_deltaH += (p.getEnergy() - p.getMomentum().z);
    }
}
