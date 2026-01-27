// Adaptive Binning Tool for Q2, beta, and x_pom
// Interactive C++ ROOT program to devise optimal binning scheme
//
// Compile: g++ AdaptiveBinning.cpp -o AdaptiveBinning $(root-config --cflags --glibs)
// Run: ./AdaptiveBinning
//
// Author: Generated for EIC TDR Analysis
// Purpose: Devise optimal binning for triple differential cross section analysis

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPad.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

// Include project headers for physics calculations
#include "Utility.hpp"
#include "RecoMethods.hpp"

// ROOT Math objects
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/GenVector/Boost.h"

using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

// Physical constants
const Float_t fMass_proton = 0.938272;
const Float_t fMass_electron = 0.000511;

// Global afterburner correction parameters
Float_t fXAngle = -0.025;
RotationX rotAboutX;
RotationY rotAboutY;
MomVector vBoostToCoM;
MomVector vBoostToHoF;

// Afterburner correction functions (from DDIS_Skim_Combined.cpp)
void undoAfterburnAndCalc(P3MVector& p, P3MVector& k){
    P3MVector p_beam(fXAngle*p.E(), 0., p.E(), p.M());
    P3MVector e_beam(0., 0., -k.E(), k.M());

    P3MVector CoM_boost = p_beam + e_beam;
    vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());

    p_beam = boost(p_beam, vBoostToCoM);
    e_beam = boost(e_beam, vBoostToCoM);

    Float_t fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
    Float_t fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

    rotAboutY = RotationY(fRotY);
    rotAboutX = RotationX(fRotX);

    p_beam = rotAboutY(p_beam);
    p_beam = rotAboutX(p_beam);
    e_beam = rotAboutY(e_beam);
    e_beam = rotAboutX(e_beam);

    P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
    vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

    p_beam = boost(p_beam, vBoostToHoF);
    e_beam = boost(e_beam, vBoostToHoF);

    p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), p_beam.E());
    k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), e_beam.E());
}

void undoAfterburn(P3MVector& a){
    a = boost(a, vBoostToCoM);
    a = rotAboutY(a);
    a = rotAboutX(a);
    a = boost(a, vBoostToHoF);
}

//==============================================================================
// STRUCTURE TO HOLD BIN EDGE INFORMATION
//==============================================================================
struct BinningScheme {
    std::vector<double> Q2_edges;
    std::vector<double> beta_edges;
    std::vector<double> xpom_edges;

    int nQ2() const { return Q2_edges.size() - 1; }
    int nBeta() const { return beta_edges.size() - 1; }
    int nXpom() const { return xpom_edges.size() - 1; }
    int nTotalBins() const { return nQ2() * nBeta() * nXpom(); }
};

//==============================================================================
// UTILITY FUNCTIONS FOR BINNING
//==============================================================================

// Generate linear bin edges
std::vector<double> GetLinearBins(double xmin, double xmax, int nbins) {
    std::vector<double> edges;
    double step = (xmax - xmin) / nbins;
    for(int i = 0; i <= nbins; i++) {
        edges.push_back(xmin + i * step);
    }
    return edges;
}

// Generate logarithmic bin edges
std::vector<double> GetLogBins(double xmin, double xmax, int nbins) {
    std::vector<double> edges;
    if(xmin <= 0 || xmax <= 0) {
        std::cout << "ERROR: Cannot use log binning with non-positive values!" << std::endl;
        return edges;
    }
    double logmin = TMath::Log10(xmin);
    double logmax = TMath::Log10(xmax);
    double step = (logmax - logmin) / nbins;
    for(int i = 0; i <= nbins; i++) {
        edges.push_back(TMath::Power(10, logmin + i * step));
    }
    return edges;
}

// Get custom bin edges from user input
std::vector<double> GetCustomBins(int nbins, const std::string& varname) {
    std::vector<double> edges;
    std::cout << "\nEnter " << (nbins + 1) << " bin edges for " << varname
              << " (space-separated, in ascending order):" << std::endl;
    std::cout << "Example: 0.0 0.2 0.4 0.6 0.8 1.0" << std::endl;
    std::cout << "> ";

    for(int i = 0; i <= nbins; i++) {
        double edge;
        std::cin >> edge;
        if(i > 0 && edge <= edges.back()) {
            std::cout << "ERROR: Bin edges must be in ascending order!" << std::endl;
            edges.clear();
            return edges;
        }
        edges.push_back(edge);
    }
    return edges;
}

// Round a value to make it "nice" (fewer decimal places)
double RoundToNiceValue(double val, int decimals = 3) {
    double factor = TMath::Power(10, decimals);
    return TMath::Floor(val * factor + 0.5) / factor;
}

// Print bin edges in a formatted way
void PrintBinEdges(const std::vector<double>& edges, const std::string& varname) {
    std::cout << "\n" << varname << " bin edges (" << edges.size() - 1 << " bins):" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for(size_t i = 0; i < edges.size(); i++) {
        std::cout << "  [" << i << "] = " << edges[i];
        if(i < edges.size() - 1) {
            std::cout << "  (width = " << std::setprecision(6) << (edges[i+1] - edges[i]) << ")";
        }
        std::cout << std::endl;
    }
}

//==============================================================================
// FUNCTION TO GET Q2 BINNING FROM USER
//==============================================================================
std::vector<double> GetQ2Binning(int& nQ2bins) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Q² BINNING SETUP" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Enter number of Q² bins: ";
    std::cin >> nQ2bins;

    std::cout << "\nSelect Q² binning type:" << std::endl;
    std::cout << "  1) Logarithmic (recommended for Q²)" << std::endl;
    std::cout << "  2) Linear" << std::endl;
    std::cout << "  3) Custom (enter your own edges)" << std::endl;
    std::cout << "Choice: ";

    int choice;
    std::cin >> choice;

    std::vector<double> edges;

    if(choice == 1) {
        double Q2min, Q2max;
        std::cout << "Enter Q² minimum (e.g., 3.4): ";
        std::cin >> Q2min;
        std::cout << "Enter Q² maximum (e.g., 150.0): ";
        std::cin >> Q2max;
        edges = GetLogBins(Q2min, Q2max, nQ2bins);
        std::cout << "Generated logarithmic Q² bins." << std::endl;
    }
    else if(choice == 2) {
        double Q2min, Q2max;
        std::cout << "Enter Q² minimum: ";
        std::cin >> Q2min;
        std::cout << "Enter Q² maximum: ";
        std::cin >> Q2max;
        edges = GetLinearBins(Q2min, Q2max, nQ2bins);
        std::cout << "Generated linear Q² bins." << std::endl;
    }
    else {
        edges = GetCustomBins(nQ2bins, "Q²");
    }

    PrintBinEdges(edges, "Q²");
    return edges;
}

//==============================================================================
// FUNCTION TO GET BETA BINNING FROM USER
//==============================================================================
std::vector<double> GetBetaBinning(int& nBetaBins) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "BETA BINNING SETUP" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Enter number of beta bins (suggested: 10): ";
    std::cin >> nBetaBins;

    std::cout << "\nSelect beta binning type:" << std::endl;
    std::cout << "  1) Linear from 0.0 to 1.0 (recommended)" << std::endl;
    std::cout << "  2) Custom range (linear)" << std::endl;
    std::cout << "  3) Custom edges" << std::endl;
    std::cout << "Choice: ";

    int choice;
    std::cin >> choice;

    std::vector<double> edges;

    if(choice == 1) {
        edges = GetLinearBins(0.0, 1.0, nBetaBins);
        std::cout << "Generated linear beta bins from 0.0 to 1.0." << std::endl;
    }
    else if(choice == 2) {
        double betaMin, betaMax;
        std::cout << "Enter beta minimum: ";
        std::cin >> betaMin;
        std::cout << "Enter beta maximum: ";
        std::cin >> betaMax;
        edges = GetLinearBins(betaMin, betaMax, nBetaBins);
        std::cout << "Generated linear beta bins." << std::endl;
    }
    else {
        edges = GetCustomBins(nBetaBins, "beta");
    }

    PrintBinEdges(edges, "beta");
    return edges;
}

//==============================================================================
// FUNCTION TO GET INITIAL X_POM BINNING FROM USER
//==============================================================================
std::vector<double> GetXpomBinning(int& nXpomBins) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "X_POM BINNING SETUP (INITIAL)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Enter number of x_pom bins: ";
    std::cin >> nXpomBins;

    std::cout << "\nSelect initial x_pom binning type:" << std::endl;
    std::cout << "  1) Logarithmic (recommended - will be adapted)" << std::endl;
    std::cout << "  2) Linear (will be adapted)" << std::endl;
    std::cout << "  3) Custom starting edges (will be adapted)" << std::endl;
    std::cout << "Choice: ";

    int choice;
    std::cin >> choice;

    std::vector<double> edges;

    if(choice == 1) {
        double xpomMin, xpomMax;
        std::cout << "Enter x_pom minimum (e.g., 0.001): ";
        std::cin >> xpomMin;
        std::cout << "Enter x_pom maximum (e.g., 0.1): ";
        std::cin >> xpomMax;
        edges = GetLogBins(xpomMin, xpomMax, nXpomBins);
        std::cout << "Generated logarithmic x_pom bins (starting point)." << std::endl;
    }
    else if(choice == 2) {
        double xpomMin, xpomMax;
        std::cout << "Enter x_pom minimum: ";
        std::cin >> xpomMin;
        std::cout << "Enter x_pom maximum: ";
        std::cin >> xpomMax;
        edges = GetLinearBins(xpomMin, xpomMax, nXpomBins);
        std::cout << "Generated linear x_pom bins (starting point)." << std::endl;
    }
    else {
        edges = GetCustomBins(nXpomBins, "x_pom");
    }

    PrintBinEdges(edges, "x_pom (initial)");
    return edges;
}

//==============================================================================
// STRUCTURE TO HOLD EVENT DATA
//==============================================================================
struct EventData {
    double Q2;
    double beta;
    double xpom;

    EventData(double q, double b, double x) : Q2(q), beta(b), xpom(x) {}
};

//==============================================================================
// FUNCTION TO READ EVENTS FROM ROOT FILE
//==============================================================================
std::vector<EventData> ReadEventsFromFile(const std::string& filename, Long64_t maxEvents = -1) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "READING EVENTS FROM FILE" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Opening file: " << filename << std::endl;

    TFile* file = TFile::Open(filename.c_str(), "READ");
    if(!file || file->IsZombie()) {
        std::cout << "ERROR: Cannot open file " << filename << std::endl;
        return {};
    }

    TTree* events = (TTree*)file->Get("events");
    if(!events) {
        std::cout << "ERROR: Cannot find 'events' tree in file" << std::endl;
        file->Close();
        return {};
    }

    Long64_t nEntries = events->GetEntries();
    std::cout << "Total entries in tree: " << nEntries << std::endl;

    if(maxEvents > 0 && maxEvents < nEntries) {
        nEntries = maxEvents;
        std::cout << "Reading only first " << nEntries << " events" << std::endl;
    }

    // Set up TTreeReader
    TTreeReader tree_reader(events);
    TTreeReaderArray<double> mc_px_array(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py_array(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz_array(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<double> mc_mass_array(tree_reader, "MCParticles.mass");
    TTreeReaderArray<int> mc_genStatus_array(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<int> mc_pdg_array(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<int> electron_scat_index(tree_reader, "ScatteredElectronsTruth_objIdx.index");

    // Read InclusiveKinematics branches (they are arrays in the file)
    TTreeReaderArray<float> electron_Q2_truth_array(tree_reader, "InclusiveKinematicsTruth.Q2");
    TTreeReaderArray<float> electron_x_truth_array(tree_reader, "InclusiveKinematicsTruth.x");
    TTreeReaderArray<float> electron_y_truth_array(tree_reader, "InclusiveKinematicsTruth.y");
    TTreeReaderArray<float> electron_W_truth_array(tree_reader, "InclusiveKinematicsTruth.W");

    std::vector<EventData> eventList;
    eventList.reserve(nEntries);

    Long64_t validEvents = 0;
    Long64_t processedEvents = 0;

    std::cout << "Reading events..." << std::endl;

    while(tree_reader.Next() && processedEvents < nEntries) {
        processedEvents++;

        if(processedEvents % 10000 == 0) {
            std::cout << "  Processed " << processedEvents << " / " << nEntries
                      << " events (" << validEvents << " valid)" << "\r" << std::flush;
        }

        // Check if InclusiveKinematics arrays have data
        if(electron_Q2_truth_array.GetSize() == 0) continue;

        // Get the first element from each array
        float electron_Q2_truth = electron_Q2_truth_array[0];
        float electron_x_truth = electron_x_truth_array[0];
        float electron_y_truth = electron_y_truth_array[0];
        float electron_W_truth = electron_W_truth_array[0];

        // Basic DIS cuts
        if(electron_Q2_truth < 1.0 || electron_x_truth <= 0.0001 || electron_x_truth >= 1.0) continue;
        if(electron_y_truth < 0.01 || electron_y_truth > 0.95) continue;
        if(electron_W_truth < 2.0) continue;

        // Find scattered electron to set up afterburner corrections
        bool foundElectron = false;
        P3MVector scattered_electron;

        if(electron_scat_index.GetSize() > 0) {
            int scat_idx = electron_scat_index[0];
            if(scat_idx >= 0 && scat_idx < mc_px_array.GetSize()) {
                if(mc_genStatus_array[scat_idx] == 1 && mc_pdg_array[scat_idx] == 11) {
                    scattered_electron.SetPxPyPzE(mc_px_array[scat_idx], mc_py_array[scat_idx],
                                                  mc_pz_array[scat_idx], mc_mass_array[scat_idx]);
                    foundElectron = true;
                }
            }
        }

        if(!foundElectron) continue;

        // Find truth proton and calculate x_pom, beta
        double xpom_from_def = -999.0;
        double beta_truth = -999.0;

        for(unsigned int j = 0; j < mc_genStatus_array.GetSize(); j++) {
            if(mc_genStatus_array[j] == 1 && mc_pdg_array[j] == 2212) {
                P3MVector proton(mc_px_array[j], mc_py_array[j], mc_pz_array[j], mc_mass_array[j]);
                undoAfterburn(proton);

                // Calculate x_pom from definition: x_pom = (M_X^2 + Q^2 - t)/(W^2 + Q^2 - m_p^2)
                BeamInfo beams;
                beams.e_beam.SetPxPyPzE(0, 0, -10.0, 10.0);  // 10 GeV electron beam
                beams.p_beam.SetPxPyPzE(0, 0, 100.0, 100.0); // 100 GeV proton beam

                double t_val = TMath::Abs(CalcT(beams.p_beam, proton));
                double m_p_sq = fMass_proton * fMass_proton;
                double W2_truth = electron_W_truth * electron_W_truth;

                // Calculate MX2 (missing mass squared)
                double MX2_truth = W2_truth + electron_Q2_truth - electron_Q2_truth * (1.0 - electron_y_truth);

                double denominator = W2_truth + electron_Q2_truth - m_p_sq;
                if(denominator > 0) {
                    xpom_from_def = (MX2_truth + electron_Q2_truth - t_val) / denominator;

                    // Calculate beta = x_Bj / x_pom
                    if(xpom_from_def > 0 && xpom_from_def < 1.0 && electron_x_truth > 0) {
                        beta_truth = electron_x_truth / xpom_from_def;

                        if(beta_truth > 0 && beta_truth <= 1.0) {
                            eventList.emplace_back(electron_Q2_truth, beta_truth, xpom_from_def);
                            validEvents++;
                        }
                    }
                }

                break; // Only process first truth proton
            }
        }
    }

    std::cout << std::endl;
    std::cout << "Finished reading events." << std::endl;
    std::cout << "Total events processed: " << processedEvents << std::endl;
    std::cout << "Valid events: " << validEvents << std::endl;

    file->Close();
    return eventList;
}

//==============================================================================
// FUNCTION TO COUNT EVENTS IN BINS
//==============================================================================
int CountEventsInBin(const std::vector<EventData>& events,
                     const BinningScheme& binning,
                     int iQ2, int iBeta, int iXpom) {
    int count = 0;

    double Q2_low = binning.Q2_edges[iQ2];
    double Q2_high = binning.Q2_edges[iQ2 + 1];
    double beta_low = binning.beta_edges[iBeta];
    double beta_high = binning.beta_edges[iBeta + 1];
    double xpom_low = binning.xpom_edges[iXpom];
    double xpom_high = binning.xpom_edges[iXpom + 1];

    for(const auto& evt : events) {
        if(evt.Q2 >= Q2_low && evt.Q2 < Q2_high &&
           evt.beta >= beta_low && evt.beta < beta_high &&
           evt.xpom >= xpom_low && evt.xpom < xpom_high) {
            count++;
        }
    }

    return count;
}

// Count events in a specific x_pom range (for adaptive binning)
int CountEventsInXpomRange(const std::vector<EventData>& events,
                           const BinningScheme& binning,
                           int iQ2, int iBeta,
                           double xpom_low, double xpom_high) {
    int count = 0;

    double Q2_low = binning.Q2_edges[iQ2];
    double Q2_high = binning.Q2_edges[iQ2 + 1];
    double beta_low = binning.beta_edges[iBeta];
    double beta_high = binning.beta_edges[iBeta + 1];

    for(const auto& evt : events) {
        if(evt.Q2 >= Q2_low && evt.Q2 < Q2_high &&
           evt.beta >= beta_low && evt.beta < beta_high &&
           evt.xpom >= xpom_low && evt.xpom < xpom_high) {
            count++;
        }
    }

    return count;
}

//==============================================================================
// ADAPTIVE BINNING ALGORITHM FOR X_POM
//==============================================================================
void AdaptiveXpomBinning(std::vector<EventData>& events, BinningScheme& binning) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "ADAPTIVE X_POM BINNING" << std::endl;
    std::cout << "========================================" << std::endl;

    int nQ2 = binning.nQ2();
    int nBeta = binning.nBeta();
    int nXpom = binning.nXpom();

    // Calculate target average events per bin
    double N_avg = (double)events.size() / (nQ2 * nBeta * nXpom);
    std::cout << "Total events: " << events.size() << std::endl;
    std::cout << "Total bins: " << nQ2 << " x " << nBeta << " x " << nXpom
              << " = " << (nQ2 * nBeta * nXpom) << std::endl;
    std::cout << "Target average events per bin: " << N_avg << std::endl;
    std::cout << "Tolerance: ±5% = [" << (N_avg * 0.95) << ", " << (N_avg * 1.05) << "]" << std::endl;

    std::cout << "\nStarting adaptive binning iteration..." << std::endl;
    std::cout << "This will adjust x_pom bin edges to achieve uniform statistics." << std::endl;

    const int MAX_ITERATIONS = 1000;
    const double TOLERANCE = 0.05; // ±5%
    const double STEP_FRACTION = 0.05; // Adjust bin width by 5% each iteration

    int iteration = 0;
    bool converged = false;

    while(!converged && iteration < MAX_ITERATIONS) {
        converged = true;
        iteration++;

        // For each Q2 and beta bin, adjust x_pom bins
        for(int iQ2 = 0; iQ2 < nQ2; iQ2++) {
            for(int iBeta = 0; iBeta < nBeta; iBeta++) {
                // Process x_pom bins from lowest to highest
                for(int iXpom = 0; iXpom < nXpom; iXpom++) {
                    double xpom_low = binning.xpom_edges[iXpom];
                    double xpom_high = binning.xpom_edges[iXpom + 1];

                    int count = CountEventsInXpomRange(events, binning, iQ2, iBeta,
                                                       xpom_low, xpom_high);

                    // Check if bin is within tolerance
                    double deviation = std::abs(count - N_avg) / N_avg;

                    if(deviation > TOLERANCE) {
                        converged = false;

                        // Calculate adjustment
                        double currentWidth = xpom_high - xpom_low;
                        double adjustment = currentWidth * STEP_FRACTION;

                        // If too many events, shrink the bin
                        if(count > N_avg) {
                            xpom_high -= adjustment;
                        }
                        // If too few events, enlarge the bin
                        else {
                            xpom_high += adjustment;
                        }

                        // Make sure we don't go beyond the overall range
                        if(iXpom < nXpom - 1) {
                            // Not the last bin - ensure we don't exceed next bin's upper edge
                            double maxEdge = binning.xpom_edges[nXpom];
                            if(xpom_high > maxEdge) xpom_high = maxEdge;
                            if(xpom_high < xpom_low + 1e-6) xpom_high = xpom_low + 1e-6;
                        }

                        // Update the edge
                        binning.xpom_edges[iXpom + 1] = xpom_high;

                        // Shift all subsequent bins to maintain continuity
                        if(iXpom < nXpom - 1) {
                            double shift = xpom_high - (binning.xpom_edges[iXpom + 1]);
                            for(int k = iXpom + 2; k <= nXpom; k++) {
                                binning.xpom_edges[k] += shift;
                            }
                        }
                    }
                }
            }
        }

        if(iteration % 100 == 0) {
            std::cout << "  Iteration " << iteration << "..." << std::endl;
        }
    }

    if(converged) {
        std::cout << "Converged after " << iteration << " iterations!" << std::endl;
    } else {
        std::cout << "Reached maximum iterations (" << MAX_ITERATIONS << "). Binning may not be fully optimized." << std::endl;
    }

    PrintBinEdges(binning.xpom_edges, "x_pom (adapted)");
}

//==============================================================================
// FUNCTION TO ALLOW MANUAL ADJUSTMENT OF BIN EDGES
//==============================================================================
void ManuallyAdjustBins(BinningScheme& binning) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "MANUAL BIN EDGE ADJUSTMENT" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Would you like to manually adjust bin edges to 'nice' values? (y/n): ";

    char response;
    std::cin >> response;

    if(response != 'y' && response != 'Y') return;

    // Adjust Q2 bins
    std::cout << "\nAdjusting Q² bin edges..." << std::endl;
    PrintBinEdges(binning.Q2_edges, "Current Q²");
    std::cout << "Round Q² edges? (y/n): ";
    std::cin >> response;
    if(response == 'y' || response == 'Y') {
        std::cout << "Enter number of decimal places (e.g., 1): ";
        int decimals;
        std::cin >> decimals;
        for(auto& edge : binning.Q2_edges) {
            edge = RoundToNiceValue(edge, decimals);
        }
        PrintBinEdges(binning.Q2_edges, "Rounded Q²");
    }

    // Adjust beta bins
    std::cout << "\nAdjusting beta bin edges..." << std::endl;
    PrintBinEdges(binning.beta_edges, "Current beta");
    std::cout << "Round beta edges? (y/n): ";
    std::cin >> response;
    if(response == 'y' || response == 'Y') {
        std::cout << "Enter number of decimal places (e.g., 2): ";
        int decimals;
        std::cin >> decimals;
        for(auto& edge : binning.beta_edges) {
            edge = RoundToNiceValue(edge, decimals);
        }
        PrintBinEdges(binning.beta_edges, "Rounded beta");
    }

    // Adjust x_pom bins
    std::cout << "\nAdjusting x_pom bin edges..." << std::endl;
    PrintBinEdges(binning.xpom_edges, "Current x_pom");
    std::cout << "Round x_pom edges? (y/n): ";
    std::cin >> response;
    if(response == 'y' || response == 'Y') {
        std::cout << "Enter number of decimal places (e.g., 4): ";
        int decimals;
        std::cin >> decimals;
        for(auto& edge : binning.xpom_edges) {
            edge = RoundToNiceValue(edge, decimals);
        }
        PrintBinEdges(binning.xpom_edges, "Rounded x_pom");
    }
}

//==============================================================================
// FUNCTION TO DISPLAY BINNING STATISTICS
//==============================================================================
void DisplayBinningStatistics(const std::vector<EventData>& events, const BinningScheme& binning) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "BINNING STATISTICS" << std::endl;
    std::cout << "========================================" << std::endl;

    int nQ2 = binning.nQ2();
    int nBeta = binning.nBeta();
    int nXpom = binning.nXpom();
    int nTotalBins = binning.nTotalBins();

    std::vector<int> counts(nTotalBins, 0);
    int minCount = INT_MAX;
    int maxCount = 0;
    int totalCount = 0;
    int nonEmptyBins = 0;

    // Count events in each bin
    for(int iQ2 = 0; iQ2 < nQ2; iQ2++) {
        for(int iBeta = 0; iBeta < nBeta; iBeta++) {
            for(int iXpom = 0; iXpom < nXpom; iXpom++) {
                int idx = iQ2 * nBeta * nXpom + iBeta * nXpom + iXpom;
                int count = CountEventsInBin(events, binning, iQ2, iBeta, iXpom);
                counts[idx] = count;
                totalCount += count;

                if(count > 0) {
                    nonEmptyBins++;
                    if(count < minCount) minCount = count;
                }
                if(count > maxCount) maxCount = count;
            }
        }
    }

    double avgCount = (double)totalCount / nTotalBins;

    std::cout << "Total bins: " << nTotalBins << std::endl;
    std::cout << "Non-empty bins: " << nonEmptyBins << std::endl;
    std::cout << "Empty bins: " << (nTotalBins - nonEmptyBins) << std::endl;
    std::cout << "Total events: " << totalCount << std::endl;
    std::cout << "Average events/bin: " << avgCount << std::endl;
    std::cout << "Min events in non-empty bin: " << (minCount == INT_MAX ? 0 : minCount) << std::endl;
    std::cout << "Max events in bin: " << maxCount << std::endl;

    // Show bins with very low statistics
    std::cout << "\nBins with < 10 events:" << std::endl;
    int lowStatBins = 0;
    for(int iQ2 = 0; iQ2 < nQ2; iQ2++) {
        for(int iBeta = 0; iBeta < nBeta; iBeta++) {
            for(int iXpom = 0; iXpom < nXpom; iXpom++) {
                int idx = iQ2 * nBeta * nXpom + iBeta * nXpom + iXpom;
                if(counts[idx] > 0 && counts[idx] < 10) {
                    std::cout << "  Q²[" << iQ2 << "], beta[" << iBeta << "], xpom[" << iXpom
                              << "]: " << counts[idx] << " events" << std::endl;
                    lowStatBins++;
                }
            }
        }
    }
    if(lowStatBins == 0) {
        std::cout << "  None" << std::endl;
    }
}

//==============================================================================
// FUNCTION TO SAVE BINNING TO FILE
//==============================================================================
void SaveBinningToFile(const BinningScheme& binning, const std::string& filename) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "SAVING BINNING TO FILE" << std::endl;
    std::cout << "========================================" << std::endl;

    std::ofstream outfile(filename);
    if(!outfile.is_open()) {
        std::cout << "ERROR: Cannot open output file " << filename << std::endl;
        return;
    }

    // Write header
    outfile << "# Adaptive Binning Scheme for Q², beta, and x_pom" << std::endl;
    outfile << "# Generated by AdaptiveBinning.cpp" << std::endl;
    outfile << "#" << std::endl;
    outfile << "# Format: Each section contains bin edges in ascending order" << std::endl;
    outfile << "#" << std::endl;

    // Write Q2 bins
    outfile << "\n# Q² bins (" << binning.nQ2() << " bins)" << std::endl;
    outfile << "Q2_BINS " << binning.Q2_edges.size() << std::endl;
    outfile << std::fixed << std::setprecision(6);
    for(size_t i = 0; i < binning.Q2_edges.size(); i++) {
        outfile << binning.Q2_edges[i];
        if(i < binning.Q2_edges.size() - 1) outfile << " ";
    }
    outfile << std::endl;

    // Write beta bins
    outfile << "\n# beta bins (" << binning.nBeta() << " bins)" << std::endl;
    outfile << "BETA_BINS " << binning.beta_edges.size() << std::endl;
    for(size_t i = 0; i < binning.beta_edges.size(); i++) {
        outfile << binning.beta_edges[i];
        if(i < binning.beta_edges.size() - 1) outfile << " ";
    }
    outfile << std::endl;

    // Write x_pom bins
    outfile << "\n# x_pom bins (" << binning.nXpom() << " bins)" << std::endl;
    outfile << "XPOM_BINS " << binning.xpom_edges.size() << std::endl;
    for(size_t i = 0; i < binning.xpom_edges.size(); i++) {
        outfile << binning.xpom_edges[i];
        if(i < binning.xpom_edges.size() - 1) outfile << " ";
    }
    outfile << std::endl;

    // Write summary
    outfile << "\n# Summary" << std::endl;
    outfile << "# Total bins: " << binning.nTotalBins()
            << " = " << binning.nQ2() << " x " << binning.nBeta() << " x " << binning.nXpom() << std::endl;

    outfile.close();
    std::cout << "Binning saved to: " << filename << std::endl;
}

//==============================================================================
// FUNCTION TO CREATE BINNING VISUALIZATION PLOTS
//==============================================================================
void CreateBinningVisualization(const std::vector<EventData>& events,
                                 const BinningScheme& binning,
                                 const std::string& outputFile = "binning_visualization.root") {
    std::cout << "\n========================================" << std::endl;
    std::cout << "CREATING BINNING VISUALIZATION" << std::endl;
    std::cout << "========================================" << std::endl;

    // Create output ROOT file
    TFile* outfile = new TFile(outputFile.c_str(), "RECREATE");
    if(!outfile || outfile->IsZombie()) {
        std::cout << "ERROR: Cannot create output file " << outputFile << std::endl;
        return;
    }

    // Set ROOT style for better plots
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // Determine ranges for histograms based on actual data
    double Q2_plot_min = binning.Q2_edges.front() * 0.9;
    double Q2_plot_max = binning.Q2_edges.back() * 1.1;
    double beta_plot_min = binning.beta_edges.front();
    double beta_plot_max = binning.beta_edges.back();
    double xpom_plot_min = binning.xpom_edges.front();
    double xpom_plot_max = binning.xpom_edges.back();

    // Create 2D histograms
    // Q2 vs x_pom (using log10(x_pom) for x-axis)
    TH2D* h_Q2_vs_xpom = new TH2D("Q2_vs_xpom",
                                   "Q^{2} vs x_{pom};log_{10}(x_{pom});Q^{2} [GeV^{2}]",
                                   100, TMath::Log10(xpom_plot_min), TMath::Log10(xpom_plot_max),
                                   100, Q2_plot_min, Q2_plot_max);

    // Q2 vs beta
    TH2D* h_Q2_vs_beta = new TH2D("Q2_vs_beta",
                                   "Q^{2} vs #beta;#beta;Q^{2} [GeV^{2}]",
                                   100, beta_plot_min, beta_plot_max,
                                   100, Q2_plot_min, Q2_plot_max);

    // beta vs x_pom (using log10(x_pom) for x-axis)
    TH2D* h_beta_vs_xpom = new TH2D("beta_vs_xpom",
                                     "#beta vs x_{pom};log_{10}(x_{pom});#beta",
                                     100, TMath::Log10(xpom_plot_min), TMath::Log10(xpom_plot_max),
                                     100, beta_plot_min, beta_plot_max);

    // Fill histograms with event data
    std::cout << "Filling histograms with " << events.size() << " events..." << std::endl;
    for(const auto& evt : events) {
        if(evt.xpom > 0) {
            h_Q2_vs_xpom->Fill(TMath::Log10(evt.xpom), evt.Q2);
            h_beta_vs_xpom->Fill(TMath::Log10(evt.xpom), evt.beta);
        }
        h_Q2_vs_beta->Fill(evt.beta, evt.Q2);
    }

    // Create canvases and draw histograms with bin lines
    std::cout << "Creating visualization plots..." << std::endl;

    // Canvas 1: Q2 vs x_pom
    TCanvas* c1 = new TCanvas("c_Q2_vs_xpom", "Q2 vs x_pom", 800, 600);
    c1->SetRightMargin(0.15);
    h_Q2_vs_xpom->Draw("COLZ");

    // Draw x_pom bin edges (vertical lines)
    for(size_t i = 0; i < binning.xpom_edges.size(); i++) {
        TLine* line = new TLine(TMath::Log10(binning.xpom_edges[i]), Q2_plot_min,
                                TMath::Log10(binning.xpom_edges[i]), Q2_plot_max);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    // Draw Q2 bin edges (horizontal lines)
    for(size_t i = 0; i < binning.Q2_edges.size(); i++) {
        TLine* line = new TLine(TMath::Log10(xpom_plot_min), binning.Q2_edges[i],
                                TMath::Log10(xpom_plot_max), binning.Q2_edges[i]);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    c1->Write();

    // Canvas 2: Q2 vs beta
    TCanvas* c2 = new TCanvas("c_Q2_vs_beta", "Q2 vs beta", 800, 600);
    c2->SetRightMargin(0.15);
    h_Q2_vs_beta->Draw("COLZ");

    // Draw beta bin edges (vertical lines)
    for(size_t i = 0; i < binning.beta_edges.size(); i++) {
        TLine* line = new TLine(binning.beta_edges[i], Q2_plot_min,
                                binning.beta_edges[i], Q2_plot_max);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    // Draw Q2 bin edges (horizontal lines)
    for(size_t i = 0; i < binning.Q2_edges.size(); i++) {
        TLine* line = new TLine(beta_plot_min, binning.Q2_edges[i],
                                beta_plot_max, binning.Q2_edges[i]);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    c2->Write();

    // Canvas 3: beta vs x_pom
    TCanvas* c3 = new TCanvas("c_beta_vs_xpom", "beta vs x_pom", 800, 600);
    c3->SetRightMargin(0.15);
    h_beta_vs_xpom->Draw("COLZ");

    // Draw x_pom bin edges (vertical lines)
    for(size_t i = 0; i < binning.xpom_edges.size(); i++) {
        TLine* line = new TLine(TMath::Log10(binning.xpom_edges[i]), beta_plot_min,
                                TMath::Log10(binning.xpom_edges[i]), beta_plot_max);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    // Draw beta bin edges (horizontal lines)
    for(size_t i = 0; i < binning.beta_edges.size(); i++) {
        TLine* line = new TLine(TMath::Log10(xpom_plot_min), binning.beta_edges[i],
                                TMath::Log10(xpom_plot_max), binning.beta_edges[i]);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    c3->Write();

    // Save histograms
    h_Q2_vs_xpom->Write();
    h_Q2_vs_beta->Write();
    h_beta_vs_xpom->Write();

    // Also save as PNG files for easy viewing
    c1->SaveAs("binning_Q2_vs_xpom.png");
    c2->SaveAs("binning_Q2_vs_beta.png");
    c3->SaveAs("binning_beta_vs_xpom.png");

    outfile->Close();

    std::cout << "Visualization plots created:" << std::endl;
    std::cout << "  ROOT file: " << outputFile << std::endl;
    std::cout << "  PNG files:" << std::endl;
    std::cout << "    - binning_Q2_vs_xpom.png" << std::endl;
    std::cout << "    - binning_Q2_vs_beta.png" << std::endl;
    std::cout << "    - binning_beta_vs_xpom.png" << std::endl;

    delete c1;
    delete c2;
    delete c3;
    delete outfile;
}

//==============================================================================
// MAIN PROGRAM
//==============================================================================
int main(int argc, char** argv) {
    std::cout << "\n================================================" << std::endl;
    std::cout << "  ADAPTIVE BINNING TOOL FOR Q², beta, x_pom" << std::endl;
    std::cout << "  EIC TDR Inclusive Diffractive DIS Analysis" << std::endl;
    std::cout << "================================================\n" << std::endl;

    // Get input filename
    std::string inputFile;
    if(argc > 1) {
        inputFile = argv[1];
    } else {
        std::cout << "Enter ROOT file name (or path): ";
        std::cin >> inputFile;
    }

    // Ask if user wants to limit number of events (for testing)
    std::cout << "Limit number of events? Enter max events or -1 for all: ";
    Long64_t maxEvents;
    std::cin >> maxEvents;

    // Read events from file
    std::vector<EventData> events = ReadEventsFromFile(inputFile, maxEvents);

    if(events.empty()) {
        std::cout << "No valid events found. Exiting." << std::endl;
        return 1;
    }

    // Set up binning scheme
    BinningScheme binning;

    // Get Q2 binning from user
    int nQ2bins;
    binning.Q2_edges = GetQ2Binning(nQ2bins);

    // Get beta binning from user
    int nBetaBins;
    binning.beta_edges = GetBetaBinning(nBetaBins);

    // Get initial x_pom binning from user
    int nXpomBins;
    binning.xpom_edges = GetXpomBinning(nXpomBins);

    // Display initial statistics
    std::cout << "\nInitial binning statistics (before adaptation):" << std::endl;
    DisplayBinningStatistics(events, binning);

    // Ask if user wants to perform adaptive binning
    std::cout << "\nPerform adaptive binning on x_pom? (y/n): ";
    char doAdaptive;
    std::cin >> doAdaptive;

    if(doAdaptive == 'y' || doAdaptive == 'Y') {
        AdaptiveXpomBinning(events, binning);

        // Display statistics after adaptation
        std::cout << "\nBinning statistics after adaptation:" << std::endl;
        DisplayBinningStatistics(events, binning);
    }

    // Allow manual adjustment of bins
    ManuallyAdjustBins(binning);

    // Final statistics
    std::cout << "\nFinal binning statistics:" << std::endl;
    DisplayBinningStatistics(events, binning);

    // Save to file
    std::string outputFile = "binning_scheme.txt";
    std::cout << "\nEnter output filename (default: binning_scheme.txt): ";
    std::string userOutput;
    std::cin >> userOutput;
    if(!userOutput.empty()) {
        outputFile = userOutput;
    }

    SaveBinningToFile(binning, outputFile);

    // Create visualization plots
    std::cout << "\nCreate visualization plots? (y/n): ";
    char createViz;
    std::cin >> createViz;

    if(createViz == 'y' || createViz == 'Y') {
        CreateBinningVisualization(events, binning);
    }

    std::cout << "\n================================================" << std::endl;
    std::cout << "  BINNING COMPLETE!" << std::endl;
    std::cout << "================================================\n" << std::endl;

    return 0;
}
