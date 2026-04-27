#pragma once

#include "RecoMethods.hpp"  // BeamInfo

#include <TTreeReaderArray.h>

// Truth-side DIS kinematics from MCParticle arrays plus beam 4-vectors.
// `scat_idx` is the MCParticle index of the scattered electron (genStatus==1).
// Returns true when x and y are finite; Q2/x/y/W are written in-place.
bool ComputeTruthKinematics(const BeamInfo& beams,
                            const TTreeReaderArray<double>& mc_px,
                            const TTreeReaderArray<double>& mc_py,
                            const TTreeReaderArray<double>& mc_pz,
                            const TTreeReaderArray<double>& mc_mass,
                            int   scat_idx,
                            float& Q2,
                            float& x,
                            float& y,
                            float& W);

// E - pz sum over reco<->truth matched particles (scattered electron and
// leading proton are NOT removed; callers that want the DIS E-pz should
// subtract them upstream).
void CalculateSumEPz_Matched(TTreeReaderArray<float>& re_px,
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
                             double& sumEPz_reco);

// Safe accessors for length-0 InclusiveKinematics* branches.
float GetArrayValue(const TTreeReaderArray<float>& arr, float fallback = -999.0f);
float GetArrayValue(const TTreeReaderArray<float>* arr, float fallback = -999.0f);
