#pragma once

// Sigmoid DA<->EM blend for W^2 reconstruction. The proxy is the geometric
// mean of the two methods; the sigmoid centres the hand-off between them.

constexpr double kW2BestW0    = 400.0;
constexpr double kW2BestSlope = 80.0;

double ComputeW2BestAlpha(double w2_proxy);

// Returns true if the blended W^2 is finite and positive. Outputs are set
// even when only one of the two inputs is valid.
bool ComputeW2Best(double w2_em,
                   double w2_da,
                   double& w2_best,
                   double& alpha,
                   double& w2_proxy);
