#ifndef TORQUE_BALANCE_OBJECTIVE_HPP
#define TORQUE_BALANCE_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 力矩平衡目标项
 * 惩罚合矩不为零: ||rL × fL + rR × fR||²
 * 权重: W_TORQUE_BAL / n_sample = 10.0 / n_sample
 */
class TorqueBalanceObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param weight 权重（默认 10.0 / n_sample）
     */
    TorqueBalanceObjective(int n_coeff, int n_sample, double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double weight_;
    bool enabled_;
    
    // 变量索引
    int idx_cbx_, idx_cby_, idx_cbz_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    int idx_fLx_, idx_fLy_, idx_fLz_;
    int idx_fRx_, idx_fRy_, idx_fRz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：计算叉积
    void crossProduct(double ax, double ay, double az,
                      double bx, double by, double bz,
                      double& cx, double& cy, double& cz) const;
};

#endif // TORQUE_BALANCE_OBJECTIVE_HPP

