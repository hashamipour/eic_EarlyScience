# MX² inputs — branches read by the skim

This is a focused inventory of the branches the skim reads for the diffractive
M_X² calculation only — both the hadronic-sum and the kinematic versions, on
truth and reco. Everything else (Q²/W/x/y kinematics, t reconstruction,
electron ID, etc.) is intentionally left out to keep the discussion narrow.

We compute **two independent reconstructions** of MX² and write both:

1. **Hadronic-sum MX²** — the textbook definition
   > M_X² = ( Σ p_i )² , where the sum runs over the hadronic system X.
   Acceptance-limited on the reco side: the central detector misses most of
   the forward-going X.

2. **Kinematic MX²** — the diffractive 4-momentum-conservation identity
   > M_X² = (q + p − p′)²,    q = k − k′
   Resolution comes from the scattered electron and the recoil proton; it does
   not require summing X particles, so it bypasses central acceptance.

---

## Truth side (MC ground truth — no detector)

### Hadronic-sum MX² truth

Source collection: **`MCParticles`**.

Branches read:

| Branch | Field used | Purpose in the MX² sum |
|---|---|---|
| `MCParticles.PDG` | int | Identify the scattered electron (PDG = 11) and the diffractive proton (PDG = 2212); skip neutrinos (|PDG| ∈ {12, 14, 16}). |
| `MCParticles.generatorStatus` | int | Keep only stable final-state particles (`generatorStatus == 1`). |
| `MCParticles.momentum.{x,y,z}` | float, float, float | Four-momentum components of each MC particle. |
| `MCParticles.mass` | float | Particle rest mass; combined with momentum to build E = √(p² + m²) for each MC particle. |

Selection in the sum (per event):
- Iterate all MC particles.
- Skip if `generatorStatus != 1`.
- Skip the scattered MC electron (matched by index `scat_mc_idx`).
- Skip the leading-pz MC proton (the diffractive proton).
- Skip MC neutrinos.
- Sum the remaining four-momenta into `total_hadrons_truth`; `MX²_truth = total_hadrons_truth.M2()`.

### Kinematic MX² truth

Source collections: **`MCParticles`** + **beam 4-vectors** (parsed at startup
from the file name / first event and undoAfterburned).

Inputs read:

| Field | Used as |
|---|---|
| Beam electron 4-vec `k` (from BeamInfo, undoAfterburned) | initial-state electron in `q = k − k′`. |
| Beam proton 4-vec `p` (from BeamInfo, undoAfterburned) | initial-state proton in the formula. |
| `MCParticles.{momentum.{x,y,z}}` at index `scat_mc_idx` (electron mass attached, undoAfterburned) | scattered MC electron `k′`. |
| `MCParticles.{momentum.{x,y,z}}` at index `lead_mc_proton_idx` (proton mass attached, undoAfterburned) | recoil MC proton `p′` — the leading-pz status==1 PDG==2212 particle. |

Computation (per event):
- `q = k − k′`
- `MX²_truth_kin = (q + p − p′)².M2()`

Histogram: `MX2_truth_kin`.

---

## Reco side (detector-level)

### Hadronic-sum reco MX²

Currently a **central-only** sum: just `ReconstructedParticles` minus the
scattered electron. The diffractive recoil proton lives in
`ForwardRomanPotRecParticles` (RP) and so is automatically absent from this
collection — no explicit proton exclusion is needed. This sum is
acceptance-limited and peaks ~10× below truth on diffractive samples; that is
the reason the kinematic MX² below was added.

#### Central detector (PFA-merged tracker + central calorimetry)

Source collection: **`ReconstructedParticles`**.

| Branch | Field used | Purpose in the MX² sum |
|---|---|---|
| `ReconstructedParticles.momentum.{x,y,z}` | float, float, float | 4-momentum components. |
| `ReconstructedParticles.energy` | float | Reconstructed energy (used directly; no recomputation from m, p). |
| `ReconstructedParticles.PDG` | int | Used only for the fallback electron-id (PDG == 11) when the truth-association lookup fails. Not used to filter hadrons in the sum. |
| `ReconstructedParticleAssociations.recID` | unsigned int | Reco ↔ MC association (rec side); used to identify the reco scattered electron index `scat_reco_idx_mx2`. |
| `ReconstructedParticleAssociations.simID` | unsigned int | Reco ↔ MC association (sim side); paired with `recID` above. |

Selection in the sum (per event):
- Iterate all entries in `ReconstructedParticles`.
- Skip the scattered reco electron (index `scat_reco_idx_mx2`).
- Add the rest into `total_hadrons_reco`.

Histogram: `MX2_reco`.

### Kinematic reco MX² (RP-tagged and B0-tagged variants)

Source collections: **`ReconstructedParticles`** (for the scattered electron)
+ recoil-proton collection (RP or B0) + beam 4-vectors (from BeamInfo,
undoAfterburned).

Inputs read:

| Field | Used as |
|---|---|
| Beam electron 4-vec `k` (from BeamInfo, undoAfterburned) | initial-state electron in `q = k − k′`. |
| Beam proton 4-vec `p` (from BeamInfo, undoAfterburned) | initial-state proton. |
| `ReconstructedParticles.{momentum.{x,y,z}}` at index `scat_reco_idx` (electron mass attached, undoAfterburned) | reco scattered electron `k′`. |
| `ForwardRomanPotRecParticles.{momentum.{x,y,z}, mass}` (per RP entry passing PDG==2212 + RP angle acceptance) | recoil proton `p′_RP`. |
| `ReconstructedTruthSeededChargedParticles.{momentum.{x,y,z}}` matched to the truth diffractive proton (proton mass attached, undoAfterburned) for entries passing the B0 angle window 5.5 < θ < 20 mrad | recoil proton `p′_B0`. |

Computation (per event, separately for RP and B0 tags when each is present):
- `q = k − k′`
- `MX²_reco_kin_RP = (q + p − p′_RP)².M2()`
- `MX²_reco_kin_B0 = (q + p − p′_B0)².M2()`

Histograms: `MX2_reco_kin_RP`, `MX2_reco_kin_B0`, and the combined `MX2_reco_kin` (one fill per event — RP value preferred, B0 value as fallback; the two acceptance regions are disjoint so there is no double-counting in practice).

---

## Read but intentionally excluded from the hadronic-sum MX²

| Collection | Why excluded |
|---|---|
| `ForwardRomanPotRecParticles` | Acceptance only catches the diffractive recoil proton (x_L ≈ 1, very close to the beam) — hadronic-system protons cannot reach RP. Used only to provide `p′_RP` for the kinematic MX² and for `t_RP` / `x_pom_RP`. |
| `ForwardOffMRecParticles`, `ReconstructedFarForwardZDCNeutrals`, `ReconstructedFarForwardZDCLambdas` | Available in the input file and physically part of X, but currently not summed into the hadronic-sum reco MX². The kinematic MX² makes the acceptance correction unnecessary, so the central-only hadronic sum is kept as a diagnostic baseline rather than expanded. Open option to revisit. |
| `ReconstructedTruthSeededChargedParticles` | Used elsewhere in the skim for the "B0 proton" identification and as the source of `p′_B0` for the kinematic MX². Not part of the hadronic-sum X sum. |
| All `InclusiveKinematics*` precomputed branches (`Truth/Electron/DA/Sigma/ESigma`) | Event-level Q²/W/x/y inputs — used for the x_pom denominator and for kinematic comparisons, not for building MX². |

---

## Summary

| Histogram | Method | Definition |
|---|---|---|
| `MX2_truth` | Hadronic sum (truth) | (Σ p_i)² over `MCParticles` minus scattered electron, leading-pz proton, neutrinos. |
| `MX2_reco` | Hadronic sum (reco) | (Σ p_i)² over central `ReconstructedParticles` minus the scattered reco electron. Acceptance-limited. |
| `MX2_truth_kin` | Kinematic (truth) | (q + p − p′)² with truth k′ and lead-pz MC proton. Should agree with `MX2_truth` (consistency check). |
| `MX2_reco_kin_RP` | Kinematic (reco, RP) | (q + p − p′)² with reco k′ and RP-reconstructed recoil proton. Acceptance-independent. |
| `MX2_reco_kin_B0` | Kinematic (reco, B0) | (q + p − p′)² with reco k′ and B0-acceptance truth-seeded recoil proton. Acceptance-independent. |
| `MX2_reco_kin` | Kinematic (reco, RP ∪ B0) | One fill per event — RP value if RP-tagged, B0 value otherwise. RP and B0 cover disjoint angular regions so there is no double-counting. This is the curve drawn on the comparison plot. |
