#include "SkimFileList.hpp"
#include "Utility.hpp"

#include <fstream>

std::string TrimLine(const std::string& line) {
    const auto start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    const auto end = line.find_last_not_of(" \t\r\n");
    return line.substr(start, end - start + 1);
}

std::vector<std::string> ReadFileList(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::string> files;
    if (!in.is_open()) {
        Logger::error("Could not open file list " + path);
        return files;
    }
    std::string line;
    while (std::getline(in, line)) {
        line = TrimLine(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        files.push_back(line);
    }
    return files;
}
