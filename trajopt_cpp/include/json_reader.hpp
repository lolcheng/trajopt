#ifndef JSON_READER_HPP
#define JSON_READER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

/**
 * 简单的JSON读取工具
 * 用于读取RBF参数
 */
class JSONReader {
public:
    /**
     * 从JSON文件读取RBF参数
     * @param filename JSON文件路径
     * @param centers_x 输出的RBF中心x坐标
     * @param centers_y 输出的RBF中心y坐标
     * @param weights 输出的RBF权重
     * @param sigma 输出的RBF sigma参数（如果JSON中没有，使用默认值）
     * @return 是否成功读取
     */
    static bool readRBFParams(const std::string& filename,
                               std::vector<double>& centers_x,
                               std::vector<double>& centers_y,
                               std::vector<double>& weights,
                               double& sigma);
    
    /**
     * 简单的JSON解析辅助函数
     * 提取数组中的数值
     */
    static std::vector<double> parseDoubleArray(const std::string& json_str, 
                                                const std::string& key);
    
    /**
     * 解析二维数组（用于cur_rbf_grid）
     */
    static bool parse2DArray(const std::string& json_str,
                             const std::string& key,
                             std::vector<double>& x,
                             std::vector<double>& y);
};

#endif // JSON_READER_HPP

