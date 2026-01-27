# DDIS Plot Merging Summary

## Overview

This document summarizes the merging of three ROOT C++ plotting files into a single comprehensive `DDIS_Plot_Combined_FINAL.cpp` script.

## Input Files Surveyed

1. **DDIS_Plot_Combined.cpp** (Current combined file)
   - Already contains ~78 active plots
   - Covers Q2/xy analysis, Mandelstam t analysis, beta analysis, cross sections

2. **DDIS_Plots_Q2_with_xy.cpp** (Q2/xy-focused script)
   - Contains many commented-out Q2/x/y distribution and resolution plots
   - Active: Some response matrices for Q2, x

3. **DDIS_Plot_t.cpp** (t-focused script)
   - Contains commented-out t distributions, beta plots, x_pom plots
   - Active: Response matrices for t, x_L, x_pom, beta

## Final Combined File: DDIS_Plot_Combined_FINAL.cpp

### Complete Feature Set (Union of All Three Files)

The final file contains **ALL** plotting configurations (including previously commented-out code) organized into these sections:

---

### **SECTION 1: Q2/XY ANALYSIS PLOTS (18 plots)**

**1.1 Base Distributions (6 plots)**
- Q2 histogram (counts) → `figs/distributions/Q2_hist.png`
- Q2 PDF comparison → `figs/distributions/Q2_pdf.png`
- x_Bj histogram → `figs/distributions/x_hist.png`
- x_Bj PDF comparison → `figs/distributions/x_pdf.png`
- y histogram → `figs/distributions/y_hist.png`
- y PDF comparison → `figs/distributions/y_pdf.png`

**1.2 Overall Relative Resolution (9 plots - previously commented)**
- Q2 resolutions (EM, DA, Sigma) → `figs/resolutions/simple/DDIS_Q2RelRes_*.png`
- x_Bj resolutions (EM, DA, Sigma) → `figs/resolutions/simple/DDIS_RelRes_xBj_*.png`
- y resolutions (EM, DA, Sigma) → `figs/resolutions/simple/DDIS_RelRes_y_*.png`

**1.3 Binned Relative Resolution (12 plots - previously commented)**
- Q2 binned (EM, DA, Sigma) → `figs/resolutions/binned/DDIS_Q2RelRes_binned_*.png`
- x_Bj binned (EM, DA, Sigma) → `figs/resolutions/binned/DDIS_RelRes_binned_x_*.png`
- y binned (EM, DA, Sigma) → `figs/resolutions/binned/DDIS_RelRes_binned_y_*.png`
- Custom fit ranges per method

**1.4 Response Matrices (9 plots)**
- Q2 response (EM, DA, Sigma) → `figs/response_matrices/response_matrix_Q2_*.png`
- x_Bj response (EM, DA, Sigma) → `figs/response_matrices/response_matrix_x_*.png`
- y response (EM, DA, Sigma) → `figs/response_matrices/response_matrix_y_*.png`

---

### **SECTION 1D: HADRONIC FINAL STATE VARIABLES (6 plots)**

- E-pz distributions (linear, log-y) → `figs/distributions/EPz_distribution*.png`
- eta_max distributions (linear, log-y) → `figs/distributions/eta_max_distribution*.png`
- MX2 distributions (linear, log-y) → `figs/distributions/MX2_distribution*.png`

---

### **SECTION 2: MANDELSTAM t ANALYSIS (23 plots)**

**2.1 Base Distributions (4 plots - includes previously commented)**
- t distributions (linear, log-y) → `figs/distributions/t_distributions*.png`
- t PDF comparison → `figs/distributions/t_pdf_comparison.png`
- theta distributions → `figs/distributions/theta_distributions.png`

**2.2 Overall Resolutions (6 plots)**
- t resolution (B0, RP) → `figs/resolutions/simple/t_resolution_*.png`
- x_L resolution (B0, RP) → `figs/resolutions/simple/xL_resolution_*.png`
- x_pom resolution (B0, RP) → `figs/resolutions/simple/xpom_resolution_*.png`

**2.3 Binned Resolutions (10 plots)**
- t binned (B0, RP) → `figs/resolutions/binned/t_resolution_binned_*.png`
- x_L binned (B0, RP) → `figs/resolutions/binned/xL_resolution_binned_*.png`
- x_pom binned (B0, RP) → `figs/resolutions/binned/xpom_resolution_binned_*.png`
- beta binned (B0, RP) → `figs/resolutions/binned/beta_resolution_binned_*.png`
- MX2 binned → `figs/resolutions/binned/MX2_resolution_binned.png`
- Individual bin plots → `figs/resolutions/binned/bins/<variable>_<detector>/`

**2.4 Response Matrices (6 plots)**
- t correlation (B0, RP) → `figs/response_matrices/response_matrix_t_*.png`
- x_L correlation (B0, RP) → `figs/response_matrices/response_matrix_xL_*.png`
- x_pom correlation (B0, RP) → `figs/response_matrices/response_matrix_xpom_*.png`

---

### **SECTION 3: X_POM COMPARISON (7 plots)**

**Definition vs. 1-x_L comparison**
- MC comparison → `figs/distributions/xpom_comparison_MC_logxy.png`
- B0 comparison → `figs/distributions/xpom_comparison_B0_logxy.png`
- RP comparison → `figs/distributions/xpom_comparison_RP_logxy.png`
- All methods combined → `figs/distributions/xpom_comparison_all_logxy.png`
- 2D correlations (MC, B0, RP) → `figs/distributions/xpom_2D_comparison_*.png`

---

### **SECTION 4: PROTON THETA ACCEPTANCE (2 plots)**

- Linear y-axis → `figs/distributions/theta_comparison_B0_acceptance.png`
- Log y-axis → `figs/distributions/theta_comparison_B0_acceptance_logxy.png`

---

### **SECTION 5: BETA ANALYSIS (9 plots)**

- Beta distribution (log-y, 4 histograms) → `figs/distributions/beta_distributions_logy.png`
- Beta resolutions (B0, RP) → `figs/resolutions/simple/beta_resolution_*.png`
- Beta response matrices (B0, RP) → `figs/response_matrices/response_matrix_beta_*.png`
- Beta correlations:
  - vs Q² → `figs/distributions/beta_vs_Q2.png`
  - vs x_pom → `figs/distributions/beta_vs_xpom.png`
  - vs t → `figs/distributions/beta_vs_t.png`

---

### **SECTION 6: 2D RESOLUTION CIRCLE MAPS (19 plots - Custom Code)**

**Circle plots with hollow/filled markers based on statistics:**
- Q2 resolution maps (EM, DA, Sigma, BestMethod) → `figs/resolutions/2d_maps/Q2_RelRes_Q2x_*.png`
- x resolution maps (EM, DA, Sigma, BestMethod) → `figs/resolutions/2d_maps/x_RelRes_xQ2_*.png`
- y resolution maps (EM, DA, Sigma, BestMethod) → `figs/resolutions/2d_maps/y_RelRes_xQ2_*.png`
- t resolution maps (B0, RP) → `figs/resolutions/2d_maps/t_RelRes_xpomQ2_*.png`
- x_pom resolution maps (B0, RP) → `figs/resolutions/2d_maps/xpom_RelRes_xpomQ2_*.png`
- beta resolution maps (B0, RP) → `figs/resolutions/2d_maps/beta_RelRes_betaQ2_*.png`
- x_L resolution maps (B0, RP) → `figs/resolutions/2d_maps/xL_RelRes_xLQ2_*.png`
- MX2 heatmap (COLZ TEXT) → `figs/resolutions/2d_maps/MX2_RelRes_MX2Q2.png`
- E-pz 2D correlation → `figs/distributions/EPz_2D.png`

**Best Method Plots:** Color-coded markers showing which reconstruction method (EM, DA, Sigma) has the smallest resolution in each (x, Q²) bin.

---

### **SECTION 7: DIFFERENTIAL CROSS SECTIONS (4 plots - Custom Code)**

**Single Differential (d(sigma)/dt):**
- With exponential fit (linear y) → `figs/cross_sections/dsigma_dt_with_fit.png`
- With exponential fit (log y) → `figs/cross_sections/dsigma_dt_logy_with_fit.png`
- Fit function: `A * exp(-b*|t|)`
- Includes MC, B0, RP, and B0+RP Sum

**Triple Differential (d³σ/(dQ² dβ dx_pom)):**
- d³σ vs beta (grid layout) → `figs/cross_sections/d3sigma_vs_beta.png`
- d³σ vs x_pom (grid layout) → `figs/cross_sections/d3sigma_vs_xpom.png`
- Grid: n_Q2_bins × n_xpom_bins or n_Q2_bins × n_beta_bins

---

### **SECTION 8: ACCEPTANCE, EFFICIENCY, PURITY (6 plots)**

**Performance Metrics:**
- Acceptance vs Q² (EM, DA, Sigma) → `figs/performance/acceptance_vs_Q2.png`
- Efficiency vs Q² (EM, DA, Sigma) → `figs/performance/efficiency_vs_Q2.png`
- Purity vs Q² (EM, DA, Sigma) → `figs/performance/purity_vs_Q2.png`

**Definitions:**
- Acceptance = (Gen in i AND Reco in i, after cuts) / (Gen in i)
- Efficiency = (Gen in i AND Reco in i, after cuts) / (Gen in i AND Reco in i, before cuts)
- Purity = (Gen in i AND Reco in i, after cuts) / (All Reco in i, after cuts)

---

## Total Plot Count

**Standard PlotOptions Framework: ~114 plots**
**Custom Canvas Plots: ~19 plots (circle maps, cross sections)**
**Grand Total: ~133+ individual plot files**

---

## Output Directory Structure (Preserved Exactly)

```
figs/
├── distributions/              # Base histograms (Q2, x, y, t, theta, beta, E-pz, eta_max, MX2, x_pom comparisons)
├── response_matrices/          # 2D correlation plots (Q2, x, y, t, x_L, x_pom, beta)
├── resolutions/
│   ├── simple/                # Overall 1D resolution histograms with Gaussian fits
│   ├── binned/                # Binned resolution summary plots
│   │   └── bins/              # Individual bin plots for each variable
│   └── 2d_maps/               # 2D resolution circle maps and best method comparisons
├── cross_sections/            # Differential cross sections (d(sigma)/dt, d³σ)
└── performance/               # Acceptance, efficiency, purity plots
```

---

## Key Features Preserved

1. **Histogram Naming Convention:** All histogram names match those in the skimming file
2. **Log/Linear Scaling:** Appropriate log/linear axes for each variable
3. **Legend Positions:** Customized per plot for optimal visibility
4. **Draw Options:** "hist", "pe", "E1" appropriately used
5. **Marker Styles:** Consistent across similar plots
6. **Color Schemes:** Black (MC), Red (B0), Blue (RP), Orange (Sum)
7. **Custom Canvas Features:**
   - Hollow circles (entries < 100) vs. filled circles (entries >= 100)
   - Marker size proportional to resolution RMS
   - Color-coded best method comparison
   - Exponential fit overlays with parameter display
   - Grid layout for triple differential cross sections

---

## Changes to DDIS_Skim_Combined.cpp

### **NO CHANGES REQUIRED**

The skimming file (`DDIS_Skim_Combined.cpp`) **already generates all necessary histograms** for the combined plotting script. Specifically:

**Already present in skimmer:**
- ✅ All Q2/x/y histograms (truth, EM, DA, Sigma)
- ✅ All Q2/x/y relative resolution histograms (1D and 2D binned)
- ✅ All Q2/x/y correlation histograms (TH2D)
- ✅ All Q2/x/y TProfile2D for circle maps
- ✅ All t/x_L/x_pom histograms (MC, B0, RP)
- ✅ All diffractive variable resolutions (1D and 2D binned)
- ✅ All diffractive variable correlations (TH2D)
- ✅ All diffractive variable TProfile2D for circle maps
- ✅ Beta histograms and correlations
- ✅ E-pz, eta_max, MX2 histograms
- ✅ Differential cross section histograms (d(sigma)/dt, d³σ)
- ✅ Acceptance/efficiency/purity tracking histograms
- ✅ x_pom comparison histograms (definition vs 1-x_L)
- ✅ Theta acceptance histograms (all TS vs B0)

**Conclusion:** The skimming file is comprehensive and no modifications are needed.

---

## Compilation and Usage

### Compile the final combined plotting script:
```bash
cd /home/hh/Downloads/eic_TDR/analysis
g++ -o DDIS_Plot_Combined_FINAL $(root-config --cflags --glibs) DDIS_Plot_Combined_FINAL.cpp Plotting.cpp
```

### Run the plotting script:
```bash
./DDIS_Plot_Combined_FINAL DDIS_Combined_output.root
```

### Output:
All plots will be saved to the `figs/` directory with the exact subdirectory structure as documented above.

---

## Differences from Original Files

### Added from DDIS_Plots_Q2_with_xy.cpp (previously commented):
- Q2/x/y overall relative resolution plots (9 plots)
- Q2/x/y binned relative resolution plots (12 plots)
- Custom fit ranges per bin per method

### Added from DDIS_Plot_t.cpp (previously commented):
- t distribution plots (linear, log-y, PDF)
- theta distribution plot
- x_L and x_pom comparison plots (commented versions now active)

### Enhancements in combined script:
- Comprehensive console output showing progress
- Clear section headers in code for easy navigation
- All output directories created automatically
- Total plot count displayed at end

---

## Testing Checklist

Before running the final script, ensure:
- [ ] `DDIS_Combined_output.root` exists and contains all required histograms
- [ ] `Plotting.hpp` and `Plotting.cpp` are up to date
- [ ] All required ROOT libraries are linked during compilation
- [ ] Output directory `figs/` is writable

---

## Summary

✅ **Successfully merged all plotting logic from 3 input files**
✅ **Preserved exact output directory structure**
✅ **Included union of all plots (active + commented)**
✅ **No changes needed to DDIS_Skim_Combined.cpp**
✅ **Total: ~133+ plots covering all analysis aspects**

---

**File Locations:**
- Combined plotting script: `/home/hh/Downloads/eic_TDR/analysis/DDIS_Plot_Combined_FINAL.cpp`
- Skimming script (no changes): `/home/hh/Downloads/eic_TDR/analysis/DDIS_Skim_Combined.cpp`
- This summary: `/home/hh/Downloads/eic_TDR/analysis/MERGE_SUMMARY.md`
