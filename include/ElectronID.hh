#ifndef ELECTRONID_HH
#define ELECTRONID_HH

#include "podio/Frame.h"

#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"

class ElectronID {

public:

    ElectronID();
    ElectronID(double Ee, double Eh);
    ~ElectronID();

    inline void SetBeamEnergy(double Ee, double Eh) { mEe = Ee; mEh = Eh; }
    inline void SetEoPRange(double lo, double hi)   { mEoP_min = lo; mEoP_max = hi; }
    inline void SetIsolation(double isor, double isoe) { mIsoR = isor; mIsoE = isoe; }
    inline void SetDeltaHRange(double lo, double hi) { mDeltaH_min = lo; mDeltaH_max = hi; }

    void SetEvent(const podio::Frame* event);

    edm4eic::ReconstructedParticleCollection FindScatteredElectron();
    edm4eic::ReconstructedParticleCollection GetTruthReconElectron();
    edm4hep::MCParticleCollection            GetMCElectron();
    edm4eic::ReconstructedParticle           SelectHighestPT(const edm4eic::ReconstructedParticleCollection& rcparts);
    double GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp);

    void   ComputeEventDeltaH();
    double GetEventDeltaH() const { return rcpart_deltaH; }
    bool   PassEventDeltaH() const { return rcpart_deltaH > mDeltaH_min && rcpart_deltaH < mDeltaH_max; }

private:

    const podio::Frame* mEvent = nullptr;

    double mEe;
    double mEh;

    double mEoP_min;
    double mEoP_max;
    double mDeltaH_min;
    double mDeltaH_max;
    double mIsoR;
    double mIsoE;

    void CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
                                 const edm4eic::ReconstructedParticleCollection& rcparts);

    double rcpart_sum_cluster_E;
    double rcpart_lead_cluster_E;
    double rcpart_others_cone_E;   // energy of OTHER particles' clusters inside the iso cone
    double rcpart_deltaH;          // event-level Sigma(E - p_z) over ReconstructedParticles

};

#endif
