#include "json_reader.hpp"
#include <iostream>
#include <regex>
#include <algorithm>

bool JSONReader::readRBFParams(const std::string& filename,
                               std::vector<double>& centers_x,
                               std::vector<double>& centers_y,
                               std::vector<double>& weights,
                               double& sigma) {
    // 读取文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string json_content = buffer.str();
    file.close();
    
    // 解析cur_rbf_grid (二维数组)
    if (!parse2DArray(json_content, "cur_rbf_grid", centers_x, centers_y)) {
        std::cerr << "Error: Failed to parse cur_rbf_grid" << std::endl;
        return false;
    }
    
    // 解析rbf_weight (一维数组)
    weights = parseDoubleArray(json_content, "rbf_weight");
    if (weights.empty()) {
        std::cerr << "Error: Failed to parse rbf_weight" << std::endl;
        return false;
    }
    
    // 检查数组长度是否匹配
    if (centers_x.size() != weights.size()) {
        std::cerr << "Error: Mismatch between centers (" << centers_x.size() 
                  << ") and weights (" << weights.size() << ")" << std::endl;
        return false;
    }
    
    // sigma使用默认值（JSON中没有这个字段）
    sigma = 0.14;
    
    std::cout << "Loaded RBF parameters:" << std::endl;
    std::cout << "  Number of centers: " << centers_x.size() << std::endl;
    std::cout << "  Sigma: " << sigma << std::endl;
    
    return true;
}

std::vector<double> JSONReader::parseDoubleArray(const std::string& json_str, 
                                                  const std::string& key) {
    std::vector<double> result;
    
    // 查找key的位置
    std::string key_pattern = "\"" + key + "\"\\s*:\\s*\\[";
    std::regex key_regex(key_pattern);
    std::smatch match;
    
    if (!std::regex_search(json_str, match, key_regex)) {
        return result;
    }
    
    size_t start_pos = match.position() + match.length();
    
    // 找到对应的结束括号
    int bracket_count = 1;
    size_t pos = start_pos;
    while (pos < json_str.length() && bracket_count > 0) {
        if (json_str[pos] == '[') bracket_count++;
        else if (json_str[pos] == ']') bracket_count--;
        pos++;
    }
    
    if (bracket_count != 0) {
        return result;
    }
    
    // 提取数组内容
    std::string array_content = json_str.substr(start_pos, pos - start_pos - 1);
    
    // 使用正则表达式提取所有数字
    std::regex number_regex(R"(-?\d+\.?\d*(?:[eE][+-]?\d+)?)");
    std::sregex_iterator iter(array_content.begin(), array_content.end(), number_regex);
    std::sregex_iterator end;
    
    for (; iter != end; ++iter) {
        result.push_back(std::stod(iter->str()));
    }
    
    return result;
}

bool JSONReader::parse2DArray(const std::string& json_str,
                              const std::string& key,
                              std::vector<double>& x,
                              std::vector<double>& y) {
    x.clear();
    y.clear();
    
    // 查找key的位置
    std::string key_pattern = "\"" + key + "\"\\s*:\\s*\\[";
    std::regex key_regex(key_pattern);
    std::smatch match;
    
    if (!std::regex_search(json_str, match, key_regex)) {
        return false;
    }
    
    size_t start_pos = match.position() + match.length();
    
    // 找到对应的结束括号
    int bracket_count = 1;
    size_t pos = start_pos;
    while (pos < json_str.length() && bracket_count > 0) {
        if (json_str[pos] == '[') bracket_count++;
        else if (json_str[pos] == ']') bracket_count--;
        pos++;
    }
    
    if (bracket_count != 0) {
        return false;
    }
    
    // 提取数组内容
    std::string array_content = json_str.substr(start_pos, pos - start_pos - 1);
    
    // 使用正则表达式匹配所有 [x, y] 对
    std::regex pair_regex(R"(\[\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*,\s*(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\])");
    std::sregex_iterator iter(array_content.begin(), array_content.end(), pair_regex);
    std::sregex_iterator end;
    
    for (; iter != end; ++iter) {
        x.push_back(std::stod(iter->str(1)));
        y.push_back(std::stod(iter->str(2)));
    }
    
    return !x.empty();
}

