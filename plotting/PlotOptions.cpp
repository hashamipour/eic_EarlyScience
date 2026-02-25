#include "Plotting.hpp"
#include <string>
#include <sstream>
#include <filesystem>
#include <TParameter.h>
#include <TTree.h>
#include <TLatex.h>

static bool HasSuffix(const std::string& value, const std::string& suffix) {
    if (value.size() < suffix.size()) return false;
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static bool HasPrefix(const std::string& value, const std::string& prefix) {
    return value.size() >= prefix.size() &&
           value.compare(0, prefix.size(), prefix) == 0;
}

void SaveCanvas(TCanvas* canvas, const char* filename) {
    if (!canvas || !filename) return;
    std::string name(filename);
    if (HasSuffix(name, ".pdf")) {
        name = name.substr(0, name.size() - 4) + ".png";
    } else if (!HasSuffix(name, ".png")) {
        name += ".png";
    }
    if (HasPrefix(name, "./figs/") && !HasPrefix(name, "./figs/organized/")) {
        name = "./figs/organized/" + name.substr(7);
    } else if (HasPrefix(name, "figs/") && !HasPrefix(name, "figs/organized/")) {
        name = "figs/organized/" + name.substr(5);
    }

    const std::filesystem::path outPath(name);
    if (outPath.has_parent_path()) {
        std::filesystem::create_directories(outPath.parent_path());
    }
    canvas->SaveAs(name.c_str());
}

PlotOptions::~PlotOptions() {}

static long long ReadProcessedEvents(TFile* inputFile) {
    if (!inputFile) return -1;
    auto* param = dynamic_cast<TParameter<Long64_t>*>(inputFile->Get("nEventsProcessed"));
    if (param) return static_cast<long long>(param->GetVal());
    auto* tree = dynamic_cast<TTree*>(inputFile->Get("Q2_tree"));
    if (tree) return static_cast<long long>(tree->GetEntries());
    return -1;
}

static std::string FormatEventCount(long long nEvents) {
    if (nEvents <= 0) return std::string();
    const long long rounded = ((nEvents + 500) / 1000) * 1000;
    const long long thousands = rounded / 1000;
    std::ostringstream oss;
    oss << thousands << "K";
    return oss.str();
}

std::string BuildSimLabel(TFile* inputFile) {
    std::string label = "#bf{ePIC} Simulation 25.12.0";
    const long long nEvents = ReadProcessedEvents(inputFile);
    const std::string count = FormatEventCount(nEvents);
    if (!count.empty()) {
        label += " (" + count + " events)";
    }
    return label;
}

void DrawSimLabels(TFile* inputFile) {
    TLatex latex;
    latex.SetTextSize(0.033);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.12, 0.92, simLabel.c_str());
    latex.DrawLatex(0.65, 0.92, "#bf{Diff. DIS} 10x100 GeV");
}
