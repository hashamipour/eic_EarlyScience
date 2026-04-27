# eic_TDR

This repository organizes your EIC DDIS **Q²** analysis into a modular, extendable structure.

## Layout

```
.
├── analysis/                # Entry points / main programs
│   ├── DDIS_Plots_Q2_OOP.cpp
│   └── DDIS_Skim_Q2.cpp
├── include/                 # Public headers exposed to others
│   ├── Plotting.hpp
│   └── Utility.hpp
├── plotting/                # Plotting implementation (links with ROOT)
│   └── Plotting.cpp
├── utility/                 # General helpers (palettes, etc.)
│   └── Utility.cpp
├── data/
│   └── filelist.txt
├── scripts/
│   └── run_ddis.sh
├── docs/
│   └── cmmnds.txt
├── CMakeLists.txt
├── Makefile                 # Wrapper around CMake
└── README.md
```

## Build Requirements

- C++20 compiler
- [ROOT](https://root.cern/) (with `Core`, `Hist`, `Tree`, `RIO`, `Graf`, `Gpad` components)
- CMake ≥ 3.16
- `podio`, `EDM4HEP`, `EDM4EIC` — supplied by the ePIC software stack. The
  simplest way to get them is to build inside the `eic-shell` container:

  ```bash
  eic-shell
  # inside the container:
  make
  ```

Ensure `root-config` is on your `PATH` and `ROOT` environment is set (e.g., `source thisroot.sh`) if you're not using `eic-shell`.

## Configure & Build (cmake + make)

```bash
# From repository root
make            # runs cmake -S . -B build && cmake --build build -j
```

Or manually:

```bash
cmake -S . -B build
cmake --build build -j
```

Executables produced in `build/`:

- `ddis_plots_q2`
- `ddis_skim_q2`

## Run

### Skim executable

```bash
./build/ddis_skim_q2 data/filelist.txt skim_q2.root 100
```
- `data/filelist.txt` — list of input ROOT files (one per line)
- `skim_q2.root` — output skim file name
- `100` — optional max events (example placeholder)

### Plotting executable

```bash
./build/ddis_plots_q2 /path/to/input.root
```

This reads histograms/trees from `input.root` and produces plots using ROOT.


> Adjust arguments to match your program’s `usage` string printed on invalid input.

## Include Hygiene & Modularization

- **Public headers** live in `include/`. Only put declarations needed by other translation units here.  
- **Implementation** stays in `plotting/` and `utility/` as `.cpp` files.
- `analysis/*.cpp` are the **entry points** (contain `main()`), and include only what they need:
  - `#include "Plotting.hpp"` for plotting
  - `#include "Utility.hpp"` only if using shared helpers

### Suggestions applied
- Consolidated headers into `include/` and added a small `utility` module.
- Left `SetCustomPalette()` in `plotting/Plotting.cpp` (no-arg variant) and kept extended overloads in `utility/Utility.cpp` (`SetCustomPalette(std::string)` / `(int)`). This avoids symbol conflicts and lets analysis code call either style.
- Added a static library `eicplot` to link both plotting and utility into analysis executables.
- No wildcard includes; add more headers to `include/` only if you want to expose new APIs.

## Extending

- Add new helpers in `utility/*.cpp` and declare in `include/Utility.hpp`.
- Add new plotting features in `plotting/*.cpp` with declarations in `include/*.hpp`.
- Create a new analysis program as `analysis/MyStudy.cpp` (with `main`) and add an `add_executable` in `CMakeLists.txt` mirroring the others.

## Notes

If you previously used relative includes like `#include "../utilities/Utility.hpp"`, that’s no longer necessary: the compiler now searches `include/` via CMake’s `include_directories`. Use simple includes:

```cpp
#include "Utility.hpp"
#include "Plotting.hpp"
```

TODO : **tighten include hygiene** further (remove unused headers or split large headers).

## Consolidated Final Workflow (Skimmer + Plotter)

This repo now includes consolidated entry points:

- `analysis/DDIS_Skim_Final.cpp` → `ddis_skim_final`
- `analysis/DDIS_Plot_Final.cpp` → `ddis_plot_final`

### Build (CMake)

```bash
cmake -S . -B build
cmake --build build -j
```

### Run

Skim:

```bash
./build/ddis_skim_final data/filelist.txt DDIS_Final_output.root
```

Plot (standard plots only):

```bash
./build/ddis_plot_final DDIS_Final_output.root
```

Plot + binning summary (binning inputs are separate by design):

```bash
./build/ddis_plot_final DDIS_Final_output.root DDIS_BinningInputs.root data/binning_scheme.txt DDIS_BinningPlots.root
```

### Plot Inventory (by category)

**Distributions**
- `figs/distributions/Q2_hist.png`, `figs/distributions/Q2_pdf.png`
- `figs/distributions/x_hist.png`, `figs/distributions/x_pdf.png`
- `figs/distributions/y_hist.png`, `figs/distributions/y_pdf.png`
- `figs/distributions/EPz_distribution.png`, `figs/distributions/EPz_distribution_logY.png`, `figs/distributions/EPz_2D.png`
- `figs/distributions/eta_max_distribution.png`, `figs/distributions/eta_max_distribution_logY.png`
- `figs/distributions/MX2_distribution.png`, `figs/distributions/MX2_distribution_logY.png`, `figs/distributions/MX2_comparison.png`
- `figs/distributions/t_distributions.png`, `figs/distributions/t_distributions_logy.png`, `figs/distributions/t_pdf_comparison.png`
- `figs/distributions/theta_distributions.png`
- `figs/distributions/x_L_comparison.png`, `figs/distributions/x_L_comparison_logy.png`
- `figs/distributions/xpom_comparison_logx.png`, `figs/distributions/xpom_comparison_logxy.png`
- `figs/distributions/xpom_comparison_MC_logxy.png`, `figs/distributions/xpom_comparison_B0_logxy.png`, `figs/distributions/xpom_comparison_RP_logxy.png`, `figs/distributions/xpom_comparison_all_logxy.png`
- `figs/distributions/beta_distributions_logy.png`
- `figs/distributions/beta_vs_Q2.png`, `figs/distributions/beta_vs_xpom.png`, `figs/distributions/beta_vs_t.png`
- `figs/distributions/theta_comparison_B0_acceptance.png`, `figs/distributions/theta_comparison_B0_acceptance_logxy.png`
- `figs/distributions/Q2_comparison.png`
- `figs/distributions/t_B0_comparison_firstBinCut.png`, `figs/distributions/xL_B0_comparison_firstBinCut.png`, `figs/distributions/xpom_B0_comparison_firstBinCut.png`

**Resolutions (1D)**
- `figs/resolutions/simple/DDIS_Q2RelRes_EM.png`, `..._DA.png`, `..._Sigma.png`
- `figs/resolutions/simple/DDIS_RelRes_xBj_EM.png`, `..._DA.png`, `..._Sigma.png`
- `figs/resolutions/simple/DDIS_RelRes_y_EM.png`, `..._DA.png`, `..._Sigma.png`
- `figs/resolutions/simple/t_resolution_B0.png`, `figs/resolutions/simple/t_resolution_RP.png`
- `figs/resolutions/simple/xL_resolution_B0.png`, `figs/resolutions/simple/xL_resolution_RP.png`
- `figs/resolutions/simple/xpom_resolution_B0.png`, `figs/resolutions/simple/xpom_resolution_RP.png`
- `figs/resolutions/simple/beta_resolution_B0.png`, `figs/resolutions/simple/beta_resolution_RP.png`

**Binned Resolutions**
- `figs/resolutions/binned/DDIS_Q2RelRes_binned_EM.png`, `..._DA.png`, `..._Sigma.png`
- `figs/resolutions/binned/DDIS_RelRes_binned_x_EM.png`, `..._x_DA.png`, `..._x_Sigma.png`
- `figs/resolutions/binned/DDIS_RelRes_binned_y_EM.png`, `..._y_DA.png`, `..._y_Sigma.png`
- `figs/resolutions/binned/t_resolution_binned_B0.png`, `figs/resolutions/binned/t_resolution_binned_RP.png`
- `figs/resolutions/binned/xL_resolution_binned_B0.png`, `figs/resolutions/binned/xL_resolution_binned_RP.png`
- `figs/resolutions/binned/xpom_resolution_binned_B0.png`, `figs/resolutions/binned/xpom_resolution_binned_RP.png`
- `figs/resolutions/binned/beta_resolution_binned_B0.png`, `figs/resolutions/binned/beta_resolution_binned_RP.png`
- `figs/resolutions/binned/MX2_resolution_binned.png`

**Response Matrices**
- `figs/response_matrices/response_matrix_Q2_EM.png`, `..._Q2_DA.png`, `..._Q2_Sigma.png`
- `figs/response_matrices/response_matrix_x_EM.png`, `..._x_DA.png`, `..._x_Sigma.png`
- `figs/response_matrices/response_matrix_y_EM.png`, `..._y_DA.png`, `..._y_Sigma.png`
- `figs/response_matrices/response_matrix_t_B0.png`, `figs/response_matrices/response_matrix_t_RP.png`
- `figs/response_matrices/response_matrix_xL_B0.png`, `figs/response_matrices/response_matrix_xL_RP.png`
- `figs/response_matrices/response_matrix_xpom_B0.png`, `figs/response_matrices/response_matrix_xpom_RP.png`
- `figs/response_matrices/response_matrix_beta_B0.png`, `figs/response_matrices/response_matrix_beta_RP.png`
- `figs/response_matrices/response_matrix_B0_cutFirstBin.png`, `figs/response_matrices/response_matrix_xL_B0_cutFirstBin.png`, `figs/response_matrices/response_matrix_xpom_B0_cutFirstBin.png`
- `figs/response_matrices/t_correlation_combined.png`

**2D Resolution Maps**
- `figs/resolutions/2d_maps/Q2_RelRes_Q2x_EM.png`, `..._Q2x_DA.png`, `..._Q2x_Sigma.png`, `..._Q2x_BestMethod.png`
- `figs/resolutions/2d_maps/x_RelRes_xQ2_EM.png`, `..._xQ2_DA.png`, `..._xQ2_Sigma.png`, `..._xQ2_BestMethod.png`
- `figs/resolutions/2d_maps/y_RelRes_xQ2_EM.png`, `..._xQ2_DA.png`, `..._xQ2_Sigma.png`, `..._xQ2_BestMethod.png`
- `figs/resolutions/2d_maps/t_RelRes_xpomQ2_B0.png`, `..._t_RelRes_xpomQ2_RP.png`
- `figs/resolutions/2d_maps/xpom_RelRes_xpomQ2_B0.png`, `..._xpom_RelRes_xpomQ2_RP.png`
- `figs/resolutions/2d_maps/beta_RelRes_betaQ2_B0.png`, `..._beta_RelRes_betaQ2_RP.png`
- `figs/resolutions/2d_maps/xL_RelRes_xLQ2_B0.png`, `..._xL_RelRes_xLQ2_RP.png`
- `figs/resolutions/2d_maps/MX2_RelRes_MX2Q2.png`

**Cross Sections**
- `figs/cross_sections/dsigma_dt.png`, `figs/cross_sections/dsigma_dt_logy.png`
- `figs/cross_sections/dsigma_dt_with_fit.png`, `figs/cross_sections/dsigma_dt_logy_with_fit.png`
- `figs/cross_sections/d3sigma_vs_beta.png`, `figs/cross_sections/d3sigma_vs_xpom.png`

**Performance (if tracking histograms exist)**
- `figs/performance/acceptance_vs_Q2.png`
- `figs/performance/purity_vs_Q2.png`

**Binning Summary (optional)**
- `figs/binning/res_beta_B0.png`, `figs/binning/res_beta_RP.png`
- `figs/binning/res_xpom_B0.png`, `figs/binning/res_xpom_RP.png`
- `figs/binning/res_beta_vs_xpom_B0.png`, `figs/binning/res_beta_vs_xpom_RP.png`
- `figs/binning/res_Q2.png`, `figs/binning/migration.png`, `figs/binning/response_matrix.png`
- `figs/binning/purity.png`, `figs/binning/stability.png`, `figs/binning/efficiency.png`
- `figs/binning/purity_vs_Q2.png`, `figs/binning/stability_vs_Q2.png`, `figs/binning/efficiency_vs_Q2.png`
- `figs/binning/purity_beta_xpom.png`, `figs/binning/efficiency_beta_xpom.png`
- `figs/binning/binning_Q2_vs_xpom.png`, `figs/binning/binning_Q2_vs_beta.png`, `figs/binning/binning_beta_vs_xpom.png`, `figs/binning/binning_Q2_vs_y.png`
