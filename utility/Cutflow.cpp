#include "Cutflow.hpp"

#include <cmath>

bool PassDISKinematicCuts(const CutflowConfig& cfg,
                          bool validQ2,
                          bool validY,
                          double q2,
                          double y,
                          double epz,
                          bool isTruth) {
    if (!validQ2 || !validY) return false;
    if (!cfg.apply_dis_cuts) return true;

    if (!(q2 > cfg.q2_min && q2 < cfg.q2_max)) return false;
    if (!(y > cfg.y_min && y < cfg.y_max)) return false;

    const bool applyEPz = cfg.apply_epz_cut && (!isTruth || cfg.apply_epz_cut_to_truth);
    if (applyEPz) {
        if (!(std::isfinite(epz) && epz > cfg.epz_min && epz < cfg.epz_max)) return false;
    }

    return true;
}

bool PassRPAcceptance(const CutflowConfig& cfg, double theta_mrad) {
    if (!std::isfinite(theta_mrad)) return false;
    if (!cfg.apply_proton_tag_cuts) return true;
    return theta_mrad <= cfg.theta_rp_max_mrad;
}

bool PassB0Acceptance(const CutflowConfig& cfg, double theta_mrad) {
    if (!std::isfinite(theta_mrad)) return false;
    if (!cfg.apply_proton_tag_cuts) return true;
    return theta_mrad > cfg.theta_b0_min_mrad && theta_mrad < cfg.theta_b0_max_mrad;
}
