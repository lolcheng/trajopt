#ifndef FORCE_BALANCE_OBJECTIVE_HPP
#define FORCE_BALANCE_OBJECTIVE_HPP

#include "objective_base.hpp"
#include <vector>

/**
 * 力平衡目标项
 * 惩罚合力不为零: ||fL + fR + [0, 0, -mg]||²
 * 权重: W_FORCE_BAL / n_sample = 10.0 / n_sample
 */
class ForceBalanceObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param mass 质量（默认 1.0）
     * @param gravity 重力加速度（默认 9.81）
     * @param weight 权重（默认 10.0 / n_sample）
     */
    ForceBalanceObjective(int n_coeff, int n_sample,
                          double mass = 1.0, double gravity = 9.81,
                          double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double mg_;  // MASS * GRAVITY
    double weight_;
    bool enabled_;
    
    // 变量索引
    int idx_fLx_, idx_fLy_, idx_fLz_;
    int idx_fRx_, idx_fRy_, idx_fRz_;
};

#endif // FORCE_BALANCE_OBJECTIVE_HPP

