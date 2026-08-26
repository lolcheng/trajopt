#ifndef JSON_READER_HPP
#define JSON_READER_HPP

#include <string>
#include <vector>

/**
 * 简单的 JSON 读取工具，用于读取 RBF 参数（与 trajopt_cpp 兼容）
 */
class JSONReader {
public:
    static bool readRBFParams(const std::string& filename,
                              std::vector<double>& centers_x,
                              std::vector<double>& centers_y,
                              std::vector<double>& weights,
                              double& sigma);

    static std::vector<double> parseDoubleArray(const std::string& json_str,
                                                const std::string& key);

    static bool parse2DArray(const std::string& json_str,
                             const std::string& key,
                             std::vector<double>& x,
                             std::vector<double>& y);
};

#endif
