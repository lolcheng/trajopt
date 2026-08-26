#ifndef WHEEL_AIRBORNE_PENALTY_OBJECTIVE_HPP
#define WHEEL_AIRBORNE_PENALTY_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include "polynomial.hpp"
#include <vector>
#include <memory>

/**
 * 轮子离地惩罚项
 * 惩罚轮子离地的情况，让轮子尽量贴地运行
 * 当轮子离地超过阈值时（Zl - Hl > threshold），惩罚离地距离的平方
 * 
 * 目标：minimize sum(max(0, Zl - Hl - threshold)² + max(0, Zr - Hr - threshold)²)
 */
class WheelAirbornePenaltyObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param terrain RBF地形对象
     * @param weight 权重（默认值会根据n_sample自动设置）
     * @param threshold 离地判据阈值（默认0.05m = 5cm）
     */
    WheelAirbornePenaltyObjective(int n_coeff, int n_sample,
                                  std::shared_ptr<RBFTerrain> terrain,
                                  double weight = -1.0,
                                  double threshold = 0.05);
    
    bool isEnabled() const override { return enabled_; }
    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    std::shared_ptr<RBFTerrain> terrain_;
    double weight_;
    double threshold_;  // 离地判据阈值（默认5cm）
    bool enabled_;
    
    // 变量索引
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：平滑正函数（smooth_pos）
    double smoothPos(double x) const;
};

#endif // WHEEL_AIRBORNE_PENALTY_OBJECTIVE_HPP

