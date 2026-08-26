#ifndef WHEELS_DISTANCE_CONSTRAINT_HPP
#define WHEELS_DISTANCE_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 左右轮距离约束
 * 左右轮距离: WHEELS_MIN² ≤ dist² ≤ WHEELS_MAX²
 * 其中: dist² = (Xr - Xl)² + (Yr - Yl)² + (Zr - Zl)²
 * 
 * 共 n_sample 个不等式约束
 */
class WheelsDistanceConstraint : public ConstraintBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param min_dist 最小距离（默认 0.1）
     * @param max_dist 最大距离（默认 0.8）
     */
    WheelsDistanceConstraint(int n_coeff, int n_sample,
                             double min_dist = 0.1, double max_dist = 0.8);
    
    virtual ~WheelsDistanceConstraint() override;
    
    int getNumConstraints() const override;
    
    void eval_g(const double* x, double* g) override;
    
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;
    double min_dist_sq_;  // WHEELS_MIN²
    double max_dist_sq_;  // WHEELS_MAX²
    
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

#endif // WHEELS_DISTANCE_CONSTRAINT_HPP

