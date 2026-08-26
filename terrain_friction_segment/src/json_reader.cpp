#include "json_reader.hpp"
#include <iostream>
#include <regex>
#include <fstream>
#include <sstream>

bool JSONReader::readRBFParams(const std::string& filename,
                               std::vector<double>& centers_x,
                               std::vector<double>& centers_y,
                               std::vector<double>& weights,
                               double& sigma) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string json_content = buffer.str();
    file.close();

    if (!parse2DArray(json_content, "cur_rbf_grid", centers_x, centers_y)) {
        std::cerr << "Error: Failed to parse cur_rbf_grid" << std::endl;
        return false;
    }
    weights = parseDoubleArray(json_content, "rbf_weight");
    if (weights.empty()) {
        std::cerr << "Error: Failed to parse rbf_weight" << std::endl;
        return false;
    }
    if (centers_x.size() != weights.size()) {
        std::cerr << "Error: Mismatch between centers and weights" << std::endl;
        return false;
    }
    sigma = 0.14;
    std::cout << "Loaded RBF: " << centers_x.size() << " centers, sigma=" << sigma << std::endl;
    return true;
}

std::vector<double> JSONReader::parseDoubleArray(const std::string& json_str,
                                                 const std::string& key) {
    std::vector<double> result;
    std::string key_pattern = "\"" + key + "\"\\s*:\\s*\\[";
    std::regex key_regex(key_pattern);
    std::smatch match;
    if (!std::regex_search(json_str, match, key_regex)) return result;
    size_t start_pos = match.position() + match.length();
    int bracket_count = 1;
    size_t pos = start_pos;
    while (pos < json_str.length() && bracket_count > 0) {
        if (json_str[pos] == '[') bracket_count++;
        else if (json_str[pos] == ']') bracket_count--;
        pos++;
    }
    if (bracket_count != 0) return result;
    std::string array_content = json_str.substr(start_pos, pos - start_pos - 1);
    std::regex number_regex(R"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)");
    std::sregex_iterator iter(array_content.begin(), array_content.end(), number_regex);
    std::sregex_iterator end;
    for (; iter != end; ++iter) result.push_back(std::stod(iter->str()));
    return result;
}

bool JSONReader::parse2DArray(const std::string& json_str,
                              const std::string& key,
                              std::vector<double>& x,
                              std::vector<double>& y) {
    x.clear();
    y.clear();
    std::string key_pattern = "\"" + key + "\"\\s*:\\s*\\[";
    std::regex key_regex(key_pattern);
    std::smatch match;
    if (!std::regex_search(json_str, match, key_regex)) return false;
    size_t start_pos = match.position() + match.length();
    int bracket_count = 1;
    size_t pos = start_pos;
    while (pos < json_str.length() && bracket_count > 0) {
        if (json_str[pos] == '[') bracket_count++;
        else if (json_str[pos] == ']') bracket_count--;
        pos++;
    }
    if (bracket_count != 0) return false;
    std::string array_content = json_str.substr(start_pos, pos - start_pos - 1);
    std::regex pair_regex(R"(\[\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*,\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\])");
    std::sregex_iterator iter(array_content.begin(), array_content.end(), pair_regex);
    std::sregex_iterator end;
    for (; iter != end; ++iter) {
        x.push_back(std::stod(iter->str(1)));
        y.push_back(std::stod(iter->str(2)));
    }
    return !x.empty();
}
