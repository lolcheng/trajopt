#ifndef TRAJECTORY_VISUALIZER_HPP
#define TRAJECTORY_VISUALIZER_HPP

#include "rbf_terrain.hpp"
#include <vector>
#include <string>

/**
 * 轨迹可视化类
 * 用于可视化优化后的轨迹和地形
 */
class TrajectoryVisualizer {
public:
    /**
     * 可视化轨迹和地形
     * @param terrain RBF地形对象
     * @param Xb 基座x轨迹
     * @param Yb 基座y轨迹
     * @param Zb 基座z轨迹
     * @param Psib 基座偏航角轨迹
     * @param Xl 左轮x轨迹
     * @param Yl 左轮y轨迹
     * @param Zl 左轮z轨迹
     * @param Xr 右轮x轨迹
     * @param Yr 右轮y轨迹
     * @param Zr 右轮z轨迹
     * @param terrain_data_file 地形数据文件路径（可选，如果已存在）
     * @param auto_close 是否自动关闭（默认false）
     * @param close_time 自动关闭时间（秒，默认5.0）
     */
    static void visualize(const RBFTerrain& terrain,
                         const std::vector<double>& Xb,
                         const std::vector<double>& Yb,
                         const std::vector<double>& Zb,
                         const std::vector<double>& Psib,
                         const std::vector<double>& Xl,
                         const std::vector<double>& Yl,
                         const std::vector<double>& Zl,
                         const std::vector<double>& Xr,
                         const std::vector<double>& Yr,
                         const std::vector<double>& Zr,
                         const std::string& terrain_data_file = "",
                         bool auto_close = false, double close_time = 5.0);
};

#endif // TRAJECTORY_VISUALIZER_HPP

