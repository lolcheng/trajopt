#ifndef WHEEL_BASE_DISTANCE_CONSTRAINT_HPP
#define WHEEL_BASE_DISTANCE_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 轮子-中心距离约束
 * 左轮-中心距离: WHEEL_BASE_MIN² ≤ dist²_L ≤ WHEEL_BASE_MAX²
 * 右轮-中心距离: WHEEL_BASE_MIN² ≤ dist²_R ≤ WHEEL_BASE_MAX²
 * 其中: dist² = (Xb - Xw)² + (Yb - Yw)² + (Zb - Zw)²
 * 
 * 共 2×n_sample 个不等式约束
 */
class WheelBaseDistanceConstraint : public ConstraintBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param min_dist 最小距离（默认 0.4）
     * @param max_dist 最大距离（默认 1.2）
     */
    WheelBaseDistanceConstraint(int n_coeff, int n_sample,
                                double min_dist = 0.4, double max_dist = 1.2);
    
    virtual ~WheelBaseDistanceConstraint() override;
    
    int getNumConstraints() const override;
    
    void eval_g(const double* x, double* g) override;
    
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;
    double min_dist_sq_;  // WHEEL_BASE_MIN²
    double max_dist_sq_;  // WHEEL_BASE_MAX²
    
    // 变量索引
    int idx_cbx_, idx_cby_, idx_cbz_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：计算距离平方
    double computeDistSq(double xb, double yb, double zb,
                         double xw, double yw, double zw) const;
};

#endif // WHEEL_BASE_DISTANCE_CONSTRAINT_HPP

