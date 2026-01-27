# Plot Organization Structure

All plots are now organized into logical subdirectories under `figs/`:

## Directory Structure

```
figs/
├── distributions/              # Base histograms and distributions
│   ├── Q2_hist.png
│   ├── EPz_distribution*.png
│   ├── eta_max_distribution*.png
│   ├── MX2_distribution*.png
│   ├── beta_distributions_logy.png
│   ├── theta_comparison*.png
│   ├── xpom_comparison*.png
│   └── xpom_2D_comparison*.png
│
├── response_matrices/          # Correlation plots (truth vs reco)
│   ├── response_matrix_Q2_*.png
│   ├── response_matrix_x_*.png
│   ├── response_matrix_t_*.png
│   ├── response_matrix_xL_*.png
│   ├── response_matrix_xpom_*.png
│   └── response_matrix_beta_*.png
│
├── resolutions/
│   ├── simple/                 # 1D resolution histograms
│   │   ├── t_resolution*.png
│   │   ├── xL_resolution*.png
│   │   ├── xpom_resolution*.png
│   │   └── beta_resolution*.png
│   │
│   ├── binned/                 # Binned resolution with Gaussian fits
│   │   ├── bins/               # Individual bin projections
│   │   ├── t_resolution_binned*.png
│   │   ├── xL_resolution_binned*.png
│   │   ├── xpom_resolution_binned*.png
│   │   ├── beta_resolution_binned*.png
│   │   └── MX2_resolution_binned.png
│   │
│   └── 2d_maps/                # 2D resolution maps (circle plots)
│       ├── Q2_RelRes_Q2x_*.png
│       ├── x_RelRes_xQ2_*.png
│       ├── y_RelRes_xQ2_*.png
│       ├── t_RelRes_xpomQ2_*.png
│       ├── xpom_RelRes_xpomQ2_*.png
│       ├── beta_RelRes_betaQ2_*.png
│       ├── MX2_RelRes_MX2Q2.png
│       └── xL_RelRes_xLQ2_*.png
│
├── cross_sections/             # Differential cross sections
│   ├── dsigma_dt*.png
│   ├── d3sigma_vs_beta.png
│   └── d3sigma_vs_xpom.png
│
└── performance/                # Acceptance, efficiency, purity
    ├── acceptance_vs_Q2.png
    ├── efficiency_vs_Q2.png
    └── purity_vs_Q2.png
```

## Categories Explained

### Distributions
Basic kinematic distributions showing the data quality and event characteristics.

### Response Matrices
2D correlation plots showing the relationship between truth-level and reconstructed quantities. Used to assess reconstruction quality and migration effects.

### Resolutions
- **Simple**: 1D histograms of relative resolution (reco-truth)/truth
- **Binned**: Resolution in bins of the variable itself, with Gaussian fits to extract mean and sigma
- **2D Maps**: Resolution as a function of two variables (typically the variable itself and Q²), displayed as color-coded bins

### Cross Sections
Differential cross-section measurements:
- dσ/dt: Single differential in Mandelstam t
- d³σ: Triple differential in Q², β, and x_pom

### Performance
Detector/reconstruction performance metrics showing acceptance, efficiency, and purity as a function of Q².

## Usage

After running the plotter, all plots will be automatically organized into these subdirectories. The directory structure is created automatically by the plotting code.
