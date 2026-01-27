# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **Diffractive Deep Inelastic Scattering (DDIS) analysis** for the Electron-Ion Collider (EIC). The codebase implements a two-phase analysis pipeline:

1. **Skimming**: Process raw EICRecon ROOT files → create analysis-ready histograms
2. **Plotting**: Read histograms → produce publication-quality figures

The analysis compares multiple reconstruction methods (EM, DA, Sigma for Q²/x/y; B0, RP for diffractive variables) and validates them against Monte Carlo truth data.

## Build & Run Commands

### Build System

```bash
# Standard build (from repository root)
make

# Or manually with CMake
cmake -S . -B build
cmake --build build -j

# Clean rebuild
rm -rf build/*
make
```

**Requirements:**
- C++20 compiler (GCC ≥ 9 or Clang ≥ 10)
- ROOT 6.x with components: Core, Hist, Tree, RIO, Graf, Gpad
- CMake ≥ 3.16
- Ensure `root-config` is on PATH and ROOT environment is sourced

**Build Output:** All executables are placed in `build/` directory

### Running the Analysis Pipeline

#### Phase 1: Data Skimming (Processing)

```bash
# Main combined skimmer (produces all histograms in one pass)
./build/ddis_skim_combined data/filelist.txt

# Output: DDIS_Combined_output.root (~200+ histograms)
```

Alternative specialized skimmers:
```bash
./build/ddis_skim_q2_xy data/filelist.txt    # Q²/x/y analysis only
./build/skim_t data/filelist.txt             # Mandelstam t analysis only
```

**Input:** `data/filelist.txt` contains paths to EICRecon ROOT files (one per line)

#### Phase 2: Plotting (Visualization)

```bash
# Main comprehensive plotter (~133+ plots)
./build/ddis_plot_combined_final ./DDIS_Combined_output.root

# Alternative specialized plotters
./build/ddis_plot_combined ./DDIS_Combined_output.root      # ~78 plots
./build/ddis_plot_q2_xy ./DDIS_Combined_output.root         # Q²/x/y only
./build/ddis_plot_t ./DDIS_Combined_output.root             # t analysis only
```

**Output:** Plots automatically organized in `figs/` subdirectories (see Plot Organization below)

#### Binning Tools

```bash
# Interactive 3D binning explorer (launches TBrowser-like GUI)
./build/ddis_plot_binning_interactive ./DDIS_Combined_output.root

# Visualize binning scheme
./build/ddis_plot_binning ./DDIS_Combined_output.root

# Generate optimal bin edges from data
./build/adaptive_binning <input_root_file>
```

### Quick Development Script

The `run.sh` script provides a convenient wrapper:

```bash
./run.sh  # Cleans build/, rebuilds, runs main plotter
```

Modify `run.sh` to change which executable runs after build.

## Architecture

### Code Organization

```
eic_EarlyScience/
├── analysis/           # Entry points (main programs with main())
│   ├── DDIS_Skim_Combined.cpp         # Main skimmer
│   ├── DDIS_Plot_Combined_FINAL.cpp   # Main plotter
│   └── [other specialized programs]
│
├── include/            # Public API headers (auto-included by CMake)
│   ├── Plotting.hpp    # Plot class hierarchy
│   ├── BinningUtility.hpp
│   ├── RecoMethods.hpp
│   ├── Utility.hpp
│   └── DDIS_Util.hpp
│
├── plotting/           # Plotting implementation
│   ├── PlotOptions.cpp              # Base class
│   ├── PlotOptions1D.cpp            # 1D histogram overlays
│   ├── PlotOptionsRelRes.cpp        # 1D resolution fits
│   ├── PlotOptionsBinnedRelRes.cpp  # Binned resolution
│   ├── PlotOptionsResponseMatrix.cpp # 2D correlation heatmaps
│   └── PlotOptionsCombinedCorrelation.cpp
│
├── utility/            # Support implementations
│   ├── Utility.cpp     # Color palettes, binning functions, Logger
│   └── RecoMethods.cpp # Physics calculations (CalcT, CalcXL, etc.)
│
├── data/
│   ├── filelist.txt    # Input file list for skimmers
│   └── [method-specific subdirs with fit ranges]
│
└── figs/               # Output plots (auto-created, organized by category)
```

### Data Flow Pipeline

```
EICRecon ROOT files
    ↓
[DDIS_Skim_Combined.cpp]
    • Read TTree branches (scattered electron, particles, etc.)
    • Reconstruct kinematics (Q², x, y, t, β, x_pom, M_X²)
    • Apply DIS cuts (Q² > 1, 0.01 < y < 0.95, etc.)
    • Fill histograms: 1D, 2D correlations, 3D for binned resolutions
    ↓
DDIS_Combined_output.root
    (~200+ histograms: distributions, resolutions, response matrices)
    ↓
[DDIS_Plot_Combined_FINAL.cpp]
    • Create PlotOptions* objects for each desired plot
    • Configure histograms, titles, axis ranges, legends
    • Call plot_obj->Plot(inputFile) for each
    • ROOT saves PNG files to organized figs/ subdirectories
    ↓
figs/ directory tree
    (133+ publication-quality plots)
```

### Plotting Framework Architecture

The plotting system uses **template method pattern** with polymorphic plot classes:

```cpp
PlotOptions (abstract base)
├── PlotOptions1D                    // Overlay multiple 1D histograms
├── PlotOptionsRelRes                // 1D resolution with Gaussian fit
├── PlotOptionsBinnedRelRes          // Resolution in variable bins
├── PlotOptionsResponseMatrix        // 2D truth-vs-reco correlation
└── PlotOptionsCombinedCorrelation   // Multiple 2D plots on one canvas
```

**Usage pattern in analysis programs:**

```cpp
// 1. Create plot configuration objects
std::vector<PlotOptions*> plots;

plots.push_back(new PlotOptions1D(
    {"h_Q2_truth", "h_Q2_EM", "h_Q2_DA"},      // histogram names
    {"MC truth", "Reco EM", "Reco DA"},        // legend labels
    {"hist", "pe", "pe"},                       // draw options
    "Q² Reconstruction",                        // title
    "Q² [GeV²]", "# events",                   // axis labels
    "figs/distributions/Q2_hist.png",          // output path
    true, true                                  // logX, logY
));

// 2. Execute all plots
TFile* f = TFile::Open("input.root");
for (auto* plot : plots) {
    plot->Plot(f);  // Reads histograms, draws, saves PNG
}
```

**Key points:**
- Plot classes handle all ROOT drawing, fitting, styling internally
- Each plot is a configuration object, not a standalone function
- Extensible: add new plot type by subclassing PlotOptions
- All output paths in plot objects should be relative to repository root (e.g., `figs/distributions/...`)

### Static Library Pattern

All plotting and utility code is compiled into **`libeicplot.a` static library**, linked by all executables:

```cmake
add_library(eicplot STATIC
    plotting/*.cpp
    utility/*.cpp
)

add_executable(ddis_plot_combined_final analysis/DDIS_Plot_Combined_FINAL.cpp)
target_link_libraries(ddis_plot_combined_final PRIVATE eicplot ${ROOT_LIBRARIES})
```

**Benefits:**
- Faster recompilation (plotting code built once)
- Consistent API across all programs
- Clear separation: analysis programs only contain main() and plot configurations

## Plot Organization

After running a plotter, plots are automatically organized into subdirectories under `figs/`:

```
figs/
├── distributions/           # Kinematic distributions (Q², x, y, t, etc.)
├── response_matrices/       # 2D truth-vs-reco correlation heatmaps
├── resolutions/
│   ├── simple/              # 1D resolution histograms
│   ├── binned/              # Resolution in variable bins (with Gaussian fits)
│   │   └── bins/            # Individual bin projection plots
│   └── 2d_maps/             # Resolution as function of two variables
├── cross_sections/          # dσ/dt, d³σ measurements
└── performance/             # Acceptance, efficiency, purity plots
```

**Important:** When creating PlotOptionsBinnedRelRes objects, the `binSavePrefix` parameter should be the **full relative path** without the `"figs/"` prefix, as it's added internally. Example:

```cpp
// Correct
binned_plot_ptr = new PlotOptionsBinnedRelRes(
    "t_RelRes_binned_B0",
    "resolutions/binned/bins/t_B0",  // Path for individual bin plots
    // ...
);

// Incorrect (creates figs/figs/... double path)
// "figs/resolutions/binned/bins/t_B0"
```

## Physics Variables

**Q²/x/y Analysis (DIS kinematics):**
- Q²: Photon virtuality (GeV²)
- x: Bjorken x (momentum fraction)
- y: Inelasticity (energy transfer fraction)

**Reconstruction Methods:**
- EM: Electron method
- DA: Double angle method
- Sigma: Jacquet-Blondel (E-pz) method

**Diffractive Variables:**
- t: Mandelstam t (momentum transfer, GeV²)
- β: Diffractive inelasticity
- x_L: Pz/E of forward proton
- x_pom: 1 - x_L (Pomeron momentum fraction)
- M_X²: Hadronic system mass squared (GeV²)

**Reconstruction Methods:**
- B0: Forward proton detectors (B0 tracker)
- RP: Roman Pots (far-forward stations)

## Common Development Tasks

### Adding a New Plot

1. In the appropriate plotter program (e.g., `analysis/DDIS_Plot_Combined_FINAL.cpp`), create a new plot object:

```cpp
plots.push_back(new PlotOptions1D(
    {"histogram_name"},
    {"Legend label"},
    {"pe"},
    "Plot Title",
    "X axis", "Y axis",
    "figs/distributions/my_new_plot.png"
));
```

2. Rebuild and run:
```bash
make
./build/ddis_plot_combined_final ./DDIS_Combined_output.root
```

### Adding a New Executable

1. Create source file in `analysis/` directory
2. Add to `CMakeLists.txt`:
```cmake
add_executable(my_program analysis/MyProgram.cpp)
target_link_libraries(my_program PRIVATE eicplot ${ROOT_LIBRARIES})
```
3. Rebuild: `make`

### Modifying Plotting Behavior

- For all plots of a type: Edit the corresponding class in `plotting/PlotOptions*.cpp`
- For a specific plot: Modify its configuration in the analysis program's `main()`
- After changes to `plotting/` or `utility/`: Run `make` to rebuild `libeicplot.a`

### Directory Creation

Plot programs should create output directories automatically using:

```cpp
gSystem->mkdir("figs/distributions", kTRUE);  // kTRUE = recursive
gSystem->mkdir("figs/resolutions/simple", kTRUE);
gSystem->mkdir("figs/resolutions/binned/bins", kTRUE);
```

This is typically done once at the start of the plotter's `main()` function.

## Binning System

The codebase supports multi-dimensional binning via `BinningUtility.hpp`:

```cpp
BinningScheme scheme;
scheme = ReadBinningScheme("data/steering/binning_scheme.txt");

// Get bin indices
int q2_bin = scheme.GetQ2Bin(Q2_value);
int beta_bin = scheme.GetBetaBin(beta_value);
int xpom_bin = scheme.GetXpomBin(xpom_value);

// 3D → linear index mapping
int linear_bin = scheme.GetLinearBin(q2_bin, beta_bin, xpom_bin);
```

Used for multi-differential cross sections (d³σ/dQ²dβdx_pom) and binned resolution studies.

## ROOT Integration

This codebase is **ROOT-centric**:
- All data stored in ROOT TTrees and TFiles
- Histograms: TH1D, TH2D, TH3D
- 4-vector math: ROOT::Math::PxPyPzEVector
- Tree reading: TTreeReader, TTreeReaderArray
- Graphics: TCanvas, TLegend, TLatex, TPaveText

**Important ROOT patterns:**
- Always check `TFile::IsZombie()` after opening
- Use `TFile::Get()` with type cast: `(TH1D*)file->Get("name")`
- Call `gErrorIgnoreLevel = kWarning;` to suppress ROOT info messages
- Use `gSystem->mkdir()` for directory creation (handles path separators cross-platform)

## Debugging

### Build Issues

```bash
# Clean build
rm -rf build/*
make

# Verbose build output
cmake --build build --verbose
```

### Missing Histograms

If a plotter can't find histograms:
1. Open ROOT file interactively: `root DDIS_Combined_output.root`
2. List contents: `.ls` or `file->ls()`
3. Check histogram name matches exactly what's in the code

### Plot Output Issues

- Check `figs/` subdirectories were created: `ls -R figs/`
- Verify write permissions on `figs/` directory
- Look for "cannot open image file" errors in log output
- For binned resolution plots, ensure bin output paths don't start with `"figs/"` (causes double path bug)
