#include "Plotting.hpp"
#include <string>

static bool HasSuffix(const std::string& value, const std::string& suffix) {
    if (value.size() < suffix.size()) return false;
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void SaveCanvas(TCanvas* canvas, const char* filename) {
    if (!canvas || !filename) return;
    std::string name(filename);
    if (HasSuffix(name, ".pdf")) {
        name = name.substr(0, name.size() - 4) + ".png";
    } else if (!HasSuffix(name, ".png")) {
        name += ".png";
    }
    canvas->SaveAs(name.c_str());
}

PlotOptions::~PlotOptions() {}
