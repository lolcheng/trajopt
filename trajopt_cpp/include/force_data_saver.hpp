#ifndef FORCE_DATA_SAVER_HPP
#define FORCE_DATA_SAVER_HPP

#include "rbf_terrain.hpp"
#include <vector>
#include <string>

/**
 * 力数据保存类
 * 用于保存接触力、轨迹、法向量等数据供可视化使用
 */
class ForceDataSaver {
public:
    /**
     * 可视化力数据（不保存文件，直接传递给Python脚本）
     * @param terrain RBF地形对象（用于计算法向量）
     * @param fLx 左轮x方向力
     * @param fLy 左轮y方向力
     * @param fLz 左轮z方向力
     * @param fRx 右轮x方向力
     * @param fRy 右轮y方向力
     * @param fRz 右轮z方向力
     * @param Xb 基座x轨迹
     * @param Yb 基座y轨迹
     * @param Zb 基座z轨迹
     * @param Xl 左轮x轨迹
     * @param Yl 左轮y轨迹
     * @param Zl 左轮z轨迹
     * @param Xr 右轮x轨迹
     * @param Yr 右轮y轨迹
     * @param Zr 右轮z轨迹
     * @param MASS 质量
     * @param GRAVITY 重力加速度
     * @param MU_FRICTION 摩擦系数
     * @param auto_close 是否自动关闭
     * @param close_time 自动关闭时间（秒）
     */
    static void visualize(const RBFTerrain& terrain,
                         const std::vector<double>& fLx, const std::vector<double>& fLy,
                         const std::vector<double>& fLz,
                         const std::vector<double>& fRx, const std::vector<double>& fRy,
                         const std::vector<double>& fRz,
                         const std::vector<double>& Xb, const std::vector<double>& Yb,
                         const std::vector<double>& Zb,
                         const std::vector<double>& Xl, const std::vector<double>& Yl,
                         const std::vector<double>& Zl,
                         const std::vector<double>& Xr, const std::vector<double>& Yr,
                         const std::vector<double>& Zr,
                         double MASS, double GRAVITY, double MU_FRICTION,
                         bool auto_close = false, double close_time = 5.0);
};

#endif // FORCE_DATA_SAVER_HPP

