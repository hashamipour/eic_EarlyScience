cmake_minimum_required(VERSION 3.16)
project(eic_TDR_Q2 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type to Debug for development
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif()

# Add debugging flags for Debug and RelWithDebInfo configurations
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2")

find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Graf Gpad RooFit RooFitCore)

# Some ROOT installs ship RooFit headers but don't expose usable
# ROOT::RooFit / ROOT::RooFitCore imported targets (the targets exist but
# don't carry the -lRooFit / -lRooFitCore INTERFACE_LINK_LIBRARIES). Locate
# the libraries explicitly under ROOT_LIBRARY_DIR and fall back to them.
find_library(ROOFIT_LIB     NAMES RooFit     HINTS ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARY_DIRS} REQUIRED)
find_library(ROOFITCORE_LIB NAMES RooFitCore HINTS ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARY_DIRS} REQUIRED)

# EDM dependencies (available inside eic-shell). These are required for
# ElectronID, which consumes podio::Frame and edm4eic collections.
find_package(podio   REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

add_library(eicplot STATIC
    plotting/PlotOptions.cpp
    plotting/PlotOptions1D.cpp
    plotting/PlotOptionsBinnedRelRes.cpp
    plotting/PlotOptionsCombinedCorrelation.cpp
    plotting/PlotOptionsRelRes.cpp
    plotting/PlotOptionsResponseMatrix.cpp
    plotting/PlotOptions2D.cpp

    plotting/PlotDrawing.cpp
    plotting/GridDrawing.cpp

    plotting/plots/PlotXQ2Density.cpp
    plotting/plots/PlotMX2Comparison.cpp
    plotting/plots/PlotEPzWithCuts.cpp
    plotting/plots/PlotDensityFromHist.cpp
    plotting/plots/PlotDensityFromHistWithOverlay.cpp
    plotting/plots/PlotGraphDensity.cpp
    plotting/plots/PlotPhaseSpaceSlices.cpp
    plotting/plots/PlotRecoSetUncorrected.cpp
    plotting/plots/PlotRecoSetComparison.cpp
    plotting/plots/PlotRecoSetComparisonBR.cpp
    plotting/plots/PlotRecoSetComparisonBRWithSetBUncorrected.cpp
    plotting/plots/PlotRecoSetComparisonBRSumOnly.cpp
    plotting/plots/PlotRelResVsK.cpp
    plotting/plots/PlotTBinMetric.cpp
    plotting/plots/PlotResolutionComparison.cpp
    plotting/plots/Plot3DResponseMatrix.cpp
    plotting/plots/PlotXQ2PhaseSpaceWithLines.cpp
    plotting/plots/PlotElectronIDComparison.cpp

    utility/Utility.cpp
    utility/RecoMethods.cpp
    utility/ElectronID.cpp

    utility/Cutflow.cpp
    utility/W2Best.cpp
    utility/ResolutionBinning.cpp
    utility/TruthKinematics.cpp
    utility/SkimFileList.cpp
)
target_link_libraries(eicplot PUBLIC
    ${ROOT_LIBRARIES}
    ${ROOFIT_LIB} ${ROOFITCORE_LIB}
    podio::podio
    podio::podioRootIO
    EDM4HEP::edm4hep
    EDM4EIC::edm4eic
)
target_include_directories(eicplot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)




add_executable(ddis_skim_final analysis/DDIS_Skim_Final.cpp)
target_link_libraries(ddis_skim_final PRIVATE eicplot ${ROOT_LIBRARIES} ${ROOFIT_LIB} ${ROOFITCORE_LIB})

add_executable(ddis_plot_final analysis/DDIS_Plot_Final.cpp)
target_link_libraries(ddis_plot_final PRIVATE eicplot ${ROOT_LIBRARIES} ${ROOFIT_LIB} ${ROOFITCORE_LIB})




