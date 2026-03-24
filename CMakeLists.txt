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

    utility/Utility.cpp
    utility/RecoMethods.cpp
)
target_link_libraries(eicplot PUBLIC ${ROOT_LIBRARIES})
target_include_directories(eicplot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)




add_executable(ddis_skim_final analysis/DDIS_Skim_Final.cpp)
target_link_libraries(ddis_skim_final PRIVATE eicplot ${ROOT_LIBRARIES})

add_executable(ddis_plot_final analysis/DDIS_Plot_Final.cpp)
target_link_libraries(ddis_plot_final PRIVATE eicplot ${ROOT_LIBRARIES})




