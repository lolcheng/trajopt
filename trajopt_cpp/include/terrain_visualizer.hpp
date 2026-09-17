#ifndef TERRAIN_VISUALIZER_HPP
#define TERRAIN_VISUALIZER_HPP

#include "rbf_terrain.hpp"
#include <vector>
#include <memory>

/**
 * RBF地形可视化类
 * 将地形网格写入临时文件，并调用 Python 脚本进行可视化
 */
class TerrainVisualizer {
public:
    /**
     * 可视化RBF地形
     * @param terrain RBF地形对象
     * @param nx 网格x方向分辨率（默认161）
     * @param ny 网格y方向分辨率（默认161）
     * @param margin 边界扩展量（默认0.2）
     * @param auto_close 是否自动关闭（默认false）
     * @param close_time 自动关闭时间（秒，默认5.0）
     */
    static void visualize(const RBFTerrain& terrain, 
                         int nx = 161, int ny = 161, double margin = 0.2,
                         bool auto_close = false, double close_time = 5.0);
    
    /**
     * 可视化RBF地形（指定范围）
     * @param terrain RBF地形对象
     * @param xmin x最小值
     * @param xmax x最大值
     * @param ymin y最小值
     * @param ymax y最大值
     * @param nx 网格x方向分辨率
     * @param ny 网格y方向分辨率
     * @param auto_close 是否自动关闭（默认false）
     * @param close_time 自动关闭时间（秒，默认5.0）
     */
    static void visualize(const RBFTerrain& terrain,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         int nx = 161, int ny = 161,
                         bool auto_close = false, double close_time = 5.0);
};

#endif // TERRAIN_VISUALIZER_HPP
