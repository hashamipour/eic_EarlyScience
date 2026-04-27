#pragma once

#include <string>
#include <vector>

// Small drawing/formatting primitives shared by the DDIS plotter and any other
// analyses that need log/linear bin edge vectors, range labels, or the
// green-yellow-red colour palette for 2D maps.

std::vector<double> BuildLogEdges(double minVal, double maxVal, int nSlices);
std::vector<double> BuildLinEdges(double minVal, double maxVal, int nSlices);
std::string         FormatRange(double lo, double hi);
void                SetGreenYellowRedPalette();
