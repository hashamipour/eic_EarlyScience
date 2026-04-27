# ElectronID Bugs and Today's Changes

This document summarizes (a) the bugs found in the attached `ElectronID.cc` / `ElectronID.hh` relative to the scattered-electron selection in paper §A.1, and (b) all other concrete changes made to the `eic_EarlyScience` repository today. It is written as source material for slide generation — each section is self-contained.

---

## Part 1 — ElectronID Bugs (vs. paper §A.1)

Context: the attached `ElectronID` class was provided to implement the paper's §A.1 cuts for scattered-electron identification in ePIC. We reviewed it before integrating it into the DDIS skim. The class implements most of §A.1 correctly, but carries three real bugs and one robustness issue.

### §A.1 cuts and their implementation status

| Paper cut | Implementation | Status |
|---|---|---|
| Track has negative curvature | `if (reconPart.getCharge() >= 0) continue;` | Correct (charge is the right proxy). |
| 0.9 < E_clus / p < 1.2 | `sum_cluster_E / magnitude(momentum)` | Correct. |
| E_clus / Σ(E, ΔR<0.4) > 0.9 | `sum_cluster_E / isolation_E` where `isolation_E` excludes the electron's own cluster | **Bug 1 — isolation denominator.** |
| Highest-p_T tie-break | `SelectHighestPT()` | Correct (kept as a separate method by design). |
| Optional event-level Σ(E − p_z) | `mDeltaH_min`, `mDeltaH_max`, `rcpart_deltaH` declared | **Bug 2 — DeltaH machinery broken.** |

Additional issues:
- **Bug 3 — `GetTruthReconElectron()`** dereferences `meMC[0]` without checking whether any MC electron was found → segfault on events with no status-1 electron.
- Cosmetic: the paper's ΔR formula is written with a minus sign under the square root; this is a typo in the paper. The code correctly computes `sqrt(Δη² + Δφ²)`.

### Bug 1 — Isolation denominator excludes the electron's own cluster

**Paper requirement.** E_clus / Σ(E within ΔR < 0.4) > 0.9, where the Σ in the denominator **includes** the matched electron's cluster itself. This is a near-unity requirement: the electron must dominate the cone.

**What the code did.** `isolation_E` was accumulated only over *other* particles' clusters inside the cone, so the cut took the form

    E_clus / E_others_in_cone > 0.9

which is roughly ten times looser than the paper: almost any cluster passes when there is little neighboring energy, because the denominator is small or zero.

**Fix.** The denominator must be the full cone total, i.e. own cluster + others:

```cpp
const double cone_total_E = rcpart_sum_cluster_E + rcpart_others_cone_E;
const double recon_isoE   = (cone_total_E > 0.0)
                          ? rcpart_sum_cluster_E / cone_total_E
                          : 0.0;
```

To prevent the bug from regressing, the private member was renamed `rcpart_isolation_E → rcpart_others_cone_E` so the variable's name matches what it actually holds.

### Bug 2 — Event-level Σ(E − p_z) machinery was non-functional

The class declared `mDeltaH_min`, `mDeltaH_max`, and `rcpart_deltaH`, and advertised a `SetDeltaHMin()` setter, but the machinery never ran:

1. Both constructors contained a copy-paste typo of the form

   ```cpp
   mDeltaH_min = ...;
   mDeltaH_min = ...;   // intended to be mDeltaH_max
   ```

   so `mDeltaH_max` was never initialized — it contained garbage.
2. `rcpart_deltaH` was declared but never computed — no method iterated over `ReconstructedParticles` to sum `(E − p_z)`.
3. There was no getter or `Pass…` method exposing the value to callers.

**Fix.** Constructors now set sensible defaults (`mDeltaH_min = 0.0`, `mDeltaH_max = 2.0 * mEe`) and initialize `rcpart_deltaH = 0.0`. A new `ComputeEventDeltaH()` method populates the sum, `GetEventDeltaH()` returns it, and `PassEventDeltaH()` applies the range:

```cpp
void ElectronID::ComputeEventDeltaH() {
    rcpart_deltaH = 0.0;
    if (!mEvent) return;
    const auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>(
        "ReconstructedParticles");
    for (const auto& p : rcparts) {
        rcpart_deltaH += (p.getEnergy() - p.getMomentum().z);
    }
}
```

The ill-formed `SetDeltaHMin` setter was replaced by a proper `SetDeltaHRange(lo, hi)`.

### Bug 3 — `GetTruthReconElectron()` dereferences on empty collection

On events with no status-1 MC electron, `meMC[0]` is a null dereference. Fix is a one-line guard:

```cpp
if (meMC.size() == 0) return meRecon;
```

### Integration into the skim (summary)

To let ElectronID live alongside the existing skim, the build picked up podio / EDM4HEP / EDM4EIC dependencies; `utility/ElectronID.cpp` and `include/ElectronID.hh` were added to the `eicplot` static library; the skim opens each input file with both `TTreeReader` and `podio::ROOTReader` in lockstep so ElectronID consumes the `podio::Frame` while the rest of the skim keeps reading flat branches. A new parallel block of `_eid` histograms and a per-event two-method comparison (E, pT, η, φ, ΔE, Δp_T, Δη, Δφ, and a categorical `(old_yes/no, new_yes/no)` counter) were added.

**Important framing.** The comparison is **not** "two cuts-based selectors" — the old method (`ScatteredElectronsTruth_objIdx`) is truth-seeded / MC-matched, so it is effectively an oracle. The comparison is therefore: "MC-matched oracle vs. reco-only §A.1 cuts". `"neither"` in the categorical counter is gated on `passRecoDIS`, so it means "DIS-like event with no electron from either method".

---

## Part 2 — Other changes today

### 2.1 Build: resolved linker errors on this ROOT install

The build was failing with undefined RooFit symbols for `ddis_plot_final`, then with undefined `TTreeReader` symbols for `ddis_skim_final`. Both share one root cause: on this ROOT install, the modern imported targets (`ROOT::RooFit`, `ROOT::Tree`, etc.) exist but do not carry `INTERFACE_LINK_LIBRARIES` — linking against them does not pull in the underlying `.so`.

**Fix.** Go back to `${ROOT_LIBRARIES}` for the bulk, and add a path-based fallback for RooFit via `find_library`:

```cmake
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Graf Gpad RooFit RooFitCore)

find_library(ROOFIT_LIB     NAMES RooFit
             HINTS ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARY_DIRS} REQUIRED)
find_library(ROOFITCORE_LIB NAMES RooFitCore
             HINTS ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARY_DIRS} REQUIRED)
```

Both the `eicplot` library and the executables link against `${ROOT_LIBRARIES} ${ROOFIT_LIB} ${ROOFITCORE_LIB}`. Confirmed to build and run.

### 2.2 ElectronID ↔ old method comparison plots

Added `PlotElectronIDComparison(TFile*, outDir)` to `DDIS_Plot_Final.cpp`. It emits seven plots into `figs/electron_id/`:
- 1D overlays of E, p_T, η, φ, and Σ(E − p_z) for the old vs. ElectronID electron.
- 2D correlation of p_T (old) vs. p_T (ElectronID).
- Difference histograms Δφ and Δp_T between the two electrons on events where both methods fire.

### 2.3 Diff histogram ranges and binning

Initial diff histograms `dphi_e_old_eid` and `dpT_e_old_eid` used wide ranges ([−π, π] and [−5, 5] GeV) with ~50–60 bins, giving one-bin spikes because the two methods agree at the per-mil level. Rebinned to 200 bins over tight, physical ranges:

```cpp
TH1D* h_dphi_e_old_eid = new TH1D("dphi_e_old_eid",
    "#phi_{eid} - #phi_{old};#Delta#phi [rad]; events", 200, -0.02, 0.02);
TH1D* h_dpT_e_old_eid  = new TH1D("dpT_e_old_eid",
    "p_{T}^{eid} - p_{T}^{old};#Delta p_{T} [GeV]; events", 200, -0.03, 0.03);
```

### 2.4 MX² bug (diffractive hadronic-system invariant mass)

The PI flagged two low-MX² spikes in the diffractive MX² distribution. There were two independent problems plus a real detector-acceptance asymmetry that masquerades as a bug:

1. **Truth was not computed from truth.** The "truth" MX² iterated the reco-matched MC particles (MCRecoParticleAssociation) and applied a reco-level lead-proton exclusion — it was biased toward reconstructed particles and inherited reco-level decisions.
2. **The reco-side proton exclusion was wrong for this campaign.** On this dataset the diffractive proton is reconstructed in the separate `ForwardRomanPotRecParticles` / B0 collections, **not** in central `ReconstructedParticles`. Excluding the leading-p_z proton from `ReconstructedParticles` therefore strips a hadronic-system proton (an X constituent) and produces low-MX² spikes at single-hadron mass² (π² ≈ 0.02, K² ≈ 0.25, …) — exactly the shape the PI flagged.
3. **Acceptance shift (not a bug).** Even with the bug fixed, the central-only reco sum under-represents the forward-going X system. The reco MX² peak sits ~10× below the truth peak. A fully apples-to-apples reco MX² would merge `ReconstructedParticles + ForwardRomanPotRecParticles + B0` and exclude the leading-p_z proton from the merged set. Left as a follow-up.

**Fix — split into hadronic-sum and kinematic MX², kept side by side.**

- **Hadronic-sum truth MX²**: iterate `MCParticles` with `genStatus == 1`, skip the scattered MC electron, skip neutrinos, and skip the leading-p_z MC proton (the diffractive proton).
- **Hadronic-sum reco MX²**: sum central `ReconstructedParticles` minus the scattered reco electron. Acceptance-limited — the central detector misses most of the forward-going X, so this peaks ~10× below truth on diffractive samples.
- **Kinematic MX² (new — both truth and reco)**: use the diffractive identity

  > MX² = (q + p − p′)²,    q = k − k′

  where k is the beam electron, k′ the scattered electron, p the beam proton, p′ the recoil proton (RP-tagged for `MX2_reco_kin_RP`, B0-tagged for `MX2_reco_kin_B0`). Resolution is set entirely by the well-measured electron and recoil-proton 4-vectors — no dependence on the X system, so no central-acceptance penalty.

  This is **not** circular when fed into `x_pom = (Q² + MX² − t)/(Q² + W² − m_p²)`: substituting MX² = (q + p − p′)² and t = (p − p′)² collapses the numerator to 2 q·(p − p′), reducing x_pom to the textbook diffractive identity q·(p − p′)/q·p — kinematics, not a tautology.

Histograms produced: `MX2_truth`, `MX2_reco` (hadronic sum); `MX2_truth_kin`, `MX2_reco_kin_RP`, `MX2_reco_kin_B0`, and the combined `MX2_reco_kin` (kinematic). The combined `MX2_reco_kin` is filled once per event — RP value if the event is RP-tagged, B0 value otherwise (the two acceptances are disjoint, so no double-counting).

The overlay plot `figs/diffractive/distributions/MX2_comparison.png` (and `_logy.png`) is now drawn by a dedicated `PlotMX2Comparison()` function (the PlotOptions1D auto-styler collapses every "*truth*" curve to the same color, which makes it impossible to visually separate `MX2_truth` from `MX2_truth_kin`). The hand-styled overlay shows four curves with distinct colors:
  - `MX2_truth` — black line (MC truth, hadronic sum)
  - `MX2_truth_kin` — green line (MC truth, kinematic)
  - `MX2_reco` — blue circles (reco, hadronic sum)
  - `MX2_reco_kin` — red squares (reco, kinematic, RP ∪ B0)

Truth-kin should sit on top of the truth hadronic sum (consistency check); reco-kin should sit on top of the truth peak (it is not acceptance-limited).

See the per-event MX² block in `analysis/DDIS_Skim_Final.cpp` around lines 2474–2640 and the B0 / RP per-tag fills around 2890 and 3245.

### 2.5 Harmonized y-ranges for binned relative-resolution plots

For side-by-side comparison in the analysis note, all binned relative-resolution plots of the same kinematic variable now share a common y-range:

- **Q²** (EM / DA / Sigma): `SetRangeY(-0.15, 0.10)`.
- **x** (EM / DA / Sigma): `SetRangeY(-0.05, 0.05)` (previously unset).
- **y** (EM / DA / Sigma): `SetRangeY(-0.6, 0.15)` (previously unset).
- **W²** (EM / DA / Best / Sigma): `SetRangeY(-0.4, 0.4)` (previously unset).

Edits are in `analysis/DDIS_Plot_Final.cpp` around lines 4228–4465.

### 2.6 |t| binning: explicit "nice" edges

Same motivation as x_pom. The YAML `t_bins` entry was updated from `[0.01, 0.05, 0.1, 0.3, 0.5, 1.0]` to `[0.01, 0.04, 0.1, 0.2, 0.4, 0.7, 1.0]` (6 bins, centers 0.025, 0.07, 0.15, 0.3, 0.55, 0.85 GeV²).

The non-YAML fallback in `DDIS_Skim_Final.cpp` previously stitched together `GetLogBins(1e-3, 0.5, 20)` with a high-|t| tail — producing edges like `0.001, 0.001379, 0.001903, …`. It was replaced with an explicit 26-edge / 25-bin list covering 0.001 → 2.0 GeV² that is denser at low |t| for the d σ/dt fit:

```cpp
t_bins = {
    0.001, 0.002, 0.003, 0.005, 0.007,
    0.01,  0.015, 0.02,  0.03,  0.04,
    0.05,  0.07,  0.1,   0.15,  0.2,
    0.25,  0.3,   0.4,   0.5,   0.6,
    0.7,   0.8,   0.9,   1.0,   1.25,
    1.5,   2.0
};
```

### 2.7 x_pom binning: explicit "nice" edges

The 3D cross-section `d^3σ / (dQ² dβ dx_pom)` previously used

```cpp
std::vector<Double_t> xpom_3d_bins = GetLogBins(1.0e-3, 0.1, 8);
```

which produced edges `{0.001, 0.001778, 0.003162, 0.005623, 0.01, 0.017783, 0.031623, 0.056234, 0.1}` and arithmetic bin centers `0.00139, 0.00247, 0.00439, 0.00781, 0.01389, 0.02470, 0.04393, 0.07812` — ugly in plot titles and filenames.

Replaced with explicit nice edges:

```cpp
std::vector<Double_t> xpom_3d_bins =
    {0.001, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1};
const int n_xpom_3d_bins = static_cast<int>(xpom_3d_bins.size()) - 1;
```

Arithmetic bin centers are now `0.002, 0.004, 0.0075, 0.015, 0.025, 0.04, 0.0625, 0.0875`. Cross-section filenames (built by `makeXpomTag`) become `xpom_0p002`, `xpom_0p004`, `xpom_0p0075`, `xpom_0p015`, `xpom_0p025`, `xpom_0p04`, `xpom_0p0625`, `xpom_0p0875` — replacing the previous `xpom_0p0025`, `xpom_0p0044`, `xpom_0p0078`, ….

The diagnostic startup print was also switched from `std::scientific` to `std::fixed << std::setprecision(5)` so bin edges and centers read as `0.00100`, `0.00300`, … instead of `1.000000e-03`, `3.000000e-03`, ….

---

## One-slide summary (for the deck)

- **ElectronID**: §A.1 is mostly right, but three bugs would silently invalidate it — a loose isolation cut (denominator missed the electron's own cluster), a dead Σ(E − p_z) machinery (uninitialized max, never computed), and a null-deref on events with no MC electron. All three are fixed, integrated into the skim via podio, and cross-checked against the MC-matched `ScatteredElectronsTruth_objIdx` with seven comparison plots.
- **Build**: fixed linker errors on this ROOT install by falling back to `${ROOT_LIBRARIES}` + `find_library` for RooFit.
- **MX²**: low-MX² spike came from excluding a "lead proton" that was never the diffractive proton (which lives in the RP / B0 collections). Fixed.
- **Analysis-note aesthetics**: binned-resolution y-ranges harmonized per variable; x_pom and |t| binning switched to explicit nice edges with clean bin centers and clean filenames.
