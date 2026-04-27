#pragma once

#include <limits>

// DIS / forward-proton selection cuts used by the DDIS skim.
struct CutflowConfig {
    bool apply_dis_cuts = true;                 // all cuts
    bool apply_proton_tag_cuts = true;
    bool use_em_for_dis_selection = false;      // y-method switch can override later

    double q2_min = 1.0;
    double q2_max = std::numeric_limits<double>::infinity();
    double y_min = 0.01;
    double y_max = 0.95;

    bool apply_epz_cut = true;
    bool apply_epz_cut_to_truth = false;
    double epz_min = 16.0;                      // 18x275 GeV2; for 10x100 GeV2 try 16-24 GeV
    double epz_max = 24.0;

    double theta_rp_max_mrad = 5.0;
    double theta_b0_min_mrad = 5.5;
    double theta_b0_max_mrad = 20.0;
};

bool PassDISKinematicCuts(const CutflowConfig& cfg,
                          bool validQ2,
                          bool validY,
                          double q2,
                          double y,
                          double epz,
                          bool isTruth);

bool PassRPAcceptance(const CutflowConfig& cfg, double theta_mrad);
bool PassB0Acceptance(const CutflowConfig& cfg, double theta_mrad);
