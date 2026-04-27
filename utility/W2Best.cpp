#include "W2Best.hpp"

#include <algorithm>
#include <cmath>

double ComputeW2BestAlpha(double w2_proxy) {
    if (!std::isfinite(w2_proxy) || w2_proxy <= 0.0) return 0.5;
    const double alpha = 1.0 / (1.0 + std::exp(-(w2_proxy - kW2BestW0) / kW2BestSlope));
    return std::clamp(alpha, 0.0, 1.0);
}

bool ComputeW2Best(double w2_em,
                   double w2_da,
                   double& w2_best,
                   double& alpha,
                   double& w2_proxy) {
    const bool validEM = std::isfinite(w2_em) && w2_em > 0.0;
    const bool validDA = std::isfinite(w2_da) && w2_da > 0.0;

    if (!validEM && !validDA) return false;

    if (validEM && validDA) {
        w2_proxy = std::sqrt(w2_em * w2_da);
        alpha = ComputeW2BestAlpha(w2_proxy);
        w2_best = (1.0 - alpha) * w2_da + alpha * w2_em;
        return std::isfinite(w2_best) && w2_best > 0.0;
    }

    if (validDA) {
        w2_proxy = w2_da;
        alpha = 0.0;
        w2_best = w2_da;
        return true;
    }

    w2_proxy = w2_em;
    alpha = 1.0;
    w2_best = w2_em;
    return true;
}
