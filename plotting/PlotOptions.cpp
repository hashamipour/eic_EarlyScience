#include "Plotting.hpp"
#include "RecoMethods.hpp"
#include <string>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <algorithm>
#include <TParameter.h>
#include <TTree.h>
#include <TLatex.h>

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
    const std::filesystem::path outPath(name);
    if (outPath.has_parent_path()) {
        std::filesystem::create_directories(outPath.parent_path());
    }
    canvas->SaveAs(name.c_str());
}

PlotOptions::~PlotOptions() {}

static bool InferBeamEnergiesFromInputPath(const std::string& inputPath,
                                           double& eBeamGeV,
                                           double& pBeamGeV) {
    std::vector<std::string> tokens;
    tokens.reserve(16);
    tokens.push_back(inputPath);

    std::string token;
    token.reserve(inputPath.size());
    for (char c : inputPath) {
        if (c == '/' || c == '\\') {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            continue;
        }
        token.push_back(c);
    }
    if (!token.empty()) tokens.push_back(token);

    bool found = false;
    double eRef = 0.0;
    double pRef = 0.0;
    for (const auto& tok : tokens) {
        double eCand = 0.0;
        double pCand = 0.0;
        if (!ParseBeamEnergiesFromFilename(tok, eCand, pCand)) continue;

        if (!found) {
            found = true;
            eRef = eCand;
            pRef = pCand;
            continue;
        }

        const bool sameE = std::fabs(eCand - eRef) < 1e-9;
        const bool sameP = std::fabs(pCand - pRef) < 1e-9;
        if (!sameE || !sameP) {
            return false;
        }
    }

    if (!found) return false;
    eBeamGeV = eRef;
    pBeamGeV = pRef;
    return true;
}

static bool EstimateDISsFromQ2Tree(TFile* inputFile, double& sGeV2) {
    if (!inputFile) return false;
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Q2_tree"));
    if (!tree) return false;
    if (!tree->GetBranch("Q2_truth") || !tree->GetBranch("x_truth") || !tree->GetBranch("y_truth")) {
        return false;
    }

    float q2 = -1.0f;
    float x = -1.0f;
    float y = -1.0f;
    tree->SetBranchAddress("Q2_truth", &q2);
    tree->SetBranchAddress("x_truth", &x);
    tree->SetBranchAddress("y_truth", &y);

    const Long64_t nEntries = tree->GetEntries();
    if (nEntries <= 0) return false;

    const Long64_t maxSample = 200000;
    const Long64_t stride = std::max<Long64_t>(1, nEntries / maxSample);

    std::vector<double> sValues;
    sValues.reserve(static_cast<size_t>(nEntries / stride + 1));
    for (Long64_t i = 0; i < nEntries; i += stride) {
        tree->GetEntry(i);
        if (!std::isfinite(q2) || !std::isfinite(x) || !std::isfinite(y)) continue;
        if (q2 <= 0.0f || x <= 0.0f || y <= 0.0f) continue;

        const double sCand = static_cast<double>(q2) /
                             (static_cast<double>(x) * static_cast<double>(y));
        if (!std::isfinite(sCand) || sCand <= 0.0) continue;
        sValues.push_back(sCand);
    }

    if (sValues.empty()) return false;
    auto mid = sValues.begin() + static_cast<long>(sValues.size() / 2);
    std::nth_element(sValues.begin(), mid, sValues.end());
    sGeV2 = *mid;
    return std::isfinite(sGeV2) && sGeV2 > 0.0;
}

std::string BuildSimLabel(TFile* inputFile) {
    (void)inputFile;
    return "#bf{ePIC} Simulation 25.12.0";
}

static std::string BuildSqrtSLabel(TFile* inputFile) {
    if (!inputFile) return std::string();

    double sGeV2 = -1.0;
    double eBeamGeV = 0.0;
    double pBeamGeV = 0.0;
    const std::string inputPath = inputFile->GetName() ? inputFile->GetName() : "";
    if (InferBeamEnergiesFromInputPath(inputPath, eBeamGeV, pBeamGeV)) {
        sGeV2 = 4.0 * eBeamGeV * pBeamGeV;
    } else if (!EstimateDISsFromQ2Tree(inputFile, sGeV2)) {
        return std::string();
    }

    if (!std::isfinite(sGeV2) || sGeV2 <= 0.0) return std::string();

    const long long sqrtSRounded = std::llround(std::sqrt(sGeV2));
    std::ostringstream oss;
    oss << "#sqrt{s} = " << sqrtSRounded << " GeV";
    return oss.str();
}

void DrawSimLabels(TFile* inputFile) {
    TLatex latex;
    latex.SetTextSize(0.033);
    latex.SetNDC();
    latex.SetTextColor(kBlack);
    const std::string simLabel = BuildSimLabel(inputFile);
    latex.DrawLatex(0.15, 0.92, simLabel.c_str());
    const std::string sqrtSLabel = BuildSqrtSLabel(inputFile);
    if (!sqrtSLabel.empty()) {
        latex.DrawLatex(0.65, 0.92, sqrtSLabel.c_str());
    }
}
