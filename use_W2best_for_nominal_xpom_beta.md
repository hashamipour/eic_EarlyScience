# Prompt: Use W2_best as the nominal W2 for xPom and beta reconstruction

## Goal

Modify `DDIS_Skim_Final.cpp` so that the **nominal** diffractive variables
`xPom` and `beta` are computed using `W2_best` (the blended DA/EM estimator)
instead of `W2_EM` alone. The current `W2_EM`-based variant should become a
cross-check, and the current `W2_best`-based variant should become the nominal.

## Current code behaviour

In the reconstruction block (look for the xPom denominator computation):

```cpp
// --- current NOMINAL (uses W2_EM) ---
double xpom_denominator_EM = electron_Q2_EM + W2_EM - m_p_sq;
double xpom_EM = (electron_Q2_EM + MX2_reco - t_reco) / xpom_denominator_EM;
double beta_EM = electron_Q2_EM / (electron_Q2_EM + MX2_reco - t_reco);

// --- current CROSS-CHECK (uses W2_best) ---
double xpom_denominator_W2Best = electron_Q2_EM + W2_best - m_p_sq;
double xpom_W2Best = (electron_Q2_EM + MX2_reco - t_reco) / xpom_denominator_W2Best;
double beta_W2Best = electron_Q2_EM / (electron_Q2_EM + MX2_reco - t_reco);
```

## Required changes

1. **Swap the roles**: the W2_best-based xPom should be stored as the nominal
   diffractive variables used in the 3D binning, correction factors, and all
   downstream histograms. The W2_EM-based xPom should be kept as a cross-check
   variant.

2. **Propagate through the binning/filling logic**: wherever the code assigns
   events to the 3D `(Q2, beta, xPom)` bins, make sure it uses the W2_best
   variant. This includes:
   - The `findBin3D(...)` call (or equivalent) that maps events to YAML bins
   - The migration-matrix / response-matrix filling
   - The Set-A / Set-B correction-factor computation
   - Any 1D projection histograms of xPom and beta

3. **Keep the W2_EM variant** as a named cross-check (e.g. `xpom_EM`,
   `beta_EM`) so it can still be plotted for comparison, but it should NOT be
   used for the nominal binning or correction.

4. **Beta formula**: note that `beta = Q2 / (Q2 + MX2 - t)` does NOT depend
   on W2, so both variants already give the same beta. No change is needed for
   beta itself, but verify that the code is consistent.

5. **Truth-level xPom**: check whether the truth-level xPom also uses W2. If
   the truth xPom is computed from MC truth kinematics (generator-level Q2 and
   W2), leave it unchanged. Only the reco-level xPom denominator needs to
   switch to W2_best.

## Context

The analysis note (LaTeX) already states that the nominal xPom uses W2_best:

> for W2, a smooth blend of DA and EM (W2_best) is used, which propagates
> into xPom = (Q2 + MX2 - t) / (Q2 + W2_best - mp2)

The code should match this description.

## Verification

After the change:
- Re-run the skimming on a test sample and confirm the nominal xPom
  distribution shifts slightly relative to the old W2_EM-based version
- The closure test (Set-A correction applied to Set-B) should still close
- The 2D resolution maps for xPom should use the W2_best variant
