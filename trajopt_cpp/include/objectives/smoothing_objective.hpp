#ifndef SMOOTHING_OBJECTIVE_HPP
#define SMOOTHING_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 平滑目标项
 * 对轨迹的二阶差分进行惩罚: (v3 - 2*v2 + v1)²
 * 应用于: Xl, Yl, Zl, Xr, Yr, Zr, Xb, Yb, Zb, Psi
 * 权重: 1.0 / n_sample
 */
class SmoothingObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param weight 权重（默认 1.0 / n_sample）
     */
    SmoothingObjective(int n_coeff, int n_sample, double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double weight_;
    bool enabled_;
    
    // 变量索引
    int idx_cbx_, idx_cby_, idx_cbz_, idx_cpsi_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：计算平滑惩罚
    double smoothPenalty(const double* traj) const;
};

#endif // SMOOTHING_OBJECTIVE_HPP

