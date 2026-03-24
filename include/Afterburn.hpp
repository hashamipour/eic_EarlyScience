#pragma once

#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "TMath.h"

using ROOT::Math::VectorUtil::boost;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,
                                                    ROOT::Math::DefaultCoordinateSystemTag>;

// Global afterburner correction state
inline Float_t fXAngle{-0.025};
inline RotationX rotAboutX;
inline RotationY rotAboutY;
inline MomVector vBoostToCoM;
inline MomVector vBoostToHoF;

inline const Float_t fMass_proton{0.938272};
inline const Float_t fMass_electron{0.000511};
inline double MASS_PROTON   = fMass_proton;
inline double MASS_ELECTRON = fMass_electron;

// Compute afterburner correction parameters from beam 4-vectors and apply
inline void undoAfterburnAndCalc(P3MVector& p, P3MVector& k) {
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

// Apply cached afterburner correction to a single 4-vector
inline void undoAfterburn(P3MVector& a) {
    a = boost(a, vBoostToCoM);
    a = rotAboutY(a);
    a = rotAboutX(a);
    a = boost(a, vBoostToHoF);
}
