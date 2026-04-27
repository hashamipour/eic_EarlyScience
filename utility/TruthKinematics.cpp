#include "TruthKinematics.hpp"

#include <TMath.h>

#include <cmath>

bool ComputeTruthKinematics(const BeamInfo& beams,
                            const TTreeReaderArray<double>& mc_px,
                            const TTreeReaderArray<double>& mc_py,
                            const TTreeReaderArray<double>& mc_pz,
                            const TTreeReaderArray<double>& mc_mass,
                            int   scat_idx,
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
                             double& sumEPz_reco) {
    (void)re_px;
    (void)re_py;
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
        const double E_mc = TMath::Sqrt(px_mc*px_mc + py_mc*py_mc + pz_mc*pz_mc + m_mc*m_mc);
        sumEPz_truth += (E_mc - pz_mc);
    }
}

float GetArrayValue(const TTreeReaderArray<float>& arr, float fallback) {
    return (arr.GetSize() > 0) ? arr[0] : fallback;
}

float GetArrayValue(const TTreeReaderArray<float>* arr, float fallback) {
    if (!arr || arr->GetSize() == 0) return fallback;
    return (*arr)[0];
}
