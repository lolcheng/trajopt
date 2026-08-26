#ifndef REGULARIZATION_OBJECTIVE_HPP
#define REGULARIZATION_OBJECTIVE_HPP

#include "objective_base.hpp"

/**
 * 正则化目标项
 * 对多项式系数的L2范数惩罚: ||coeff||²
 * 权重: 1.0
 */
class RegularizationObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param weight 权重（默认 1.0）
     */
    RegularizationObjective(int n_coeff, double weight = 1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    double weight_;
    bool enabled_;
    
    // 变量索引（所有10组系数）
    int idx_cbx_, idx_cby_, idx_cbz_, idx_cpsi_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
};

#endif // REGULARIZATION_OBJECTIVE_HPP

