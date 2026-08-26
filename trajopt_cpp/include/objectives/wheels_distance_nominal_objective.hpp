#ifndef WHEELS_DISTANCE_NOMINAL_OBJECTIVE_HPP
#define WHEELS_DISTANCE_NOMINAL_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 左右轮距离标称目标项
 * 惩罚偏离标称距离: (dist² - WHEELS_NOMINAL²)²
 * 参数: WHEELS_NOMINAL = 0.5
 * 权重: 1.0 / n_sample
 */
class WheelsDistanceNominalObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param nominal_dist 标称距离（默认 0.5）
     * @param weight 权重（默认 1.0 / n_sample）
     */
    WheelsDistanceNominalObjective(int n_coeff, int n_sample,
                                   double nominal_dist = 0.5,
                                   double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double nominal_dist_sq_;  // WHEELS_NOMINAL²
    double weight_;
    bool enabled_;
    
    // 变量索引
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：计算距离平方
    double computeDistSq(double xl, double yl, double zl,
                         double xr, double yr, double zr) const;
};

#endif // WHEELS_DISTANCE_NOMINAL_OBJECTIVE_HPP

