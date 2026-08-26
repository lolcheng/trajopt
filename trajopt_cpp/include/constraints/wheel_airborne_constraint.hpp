#ifndef WHEEL_AIRBORNE_CONSTRAINT_HPP
#define WHEEL_AIRBORNE_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "rbf_terrain.hpp"
#include "polynomial.hpp"
#include <vector>
#include <memory>

/**
 * 轮子离地距离约束
 * 限制轮子离地距离不超过阈值（默认10cm）
 * 
 * 约束: Zl - Hl <= max_airborne_distance
 *       Zr - Hr <= max_airborne_distance
 * 
 * 共 2×n_sample 个不等式约束
 */
class WheelAirborneConstraint : public ConstraintBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param terrain RBF地形对象
     * @param min_ground_distance 最小离地距离（默认0.0m，防止下陷）
     * @param max_airborne_distance 最大离地距离（默认0.1m = 10cm）
     */
    WheelAirborneConstraint(int n_coeff, int n_sample,
                           std::shared_ptr<RBFTerrain> terrain,
                           double min_ground_distance = 0.0,
                           double max_airborne_distance = 0.1);
    
    virtual ~WheelAirborneConstraint() override;
    
    int getNumConstraints() const override;
    
    void eval_g(const double* x, double* g) override;
    
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;
    std::shared_ptr<RBFTerrain> terrain_;
    double min_ground_distance_;   // 最小离地距离（防止下陷）
    double max_airborne_distance_; // 最大离地距离（防止跳起）
    
    // 变量索引
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和Vandermonde矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    
    void evalTrajectory(const double* coeff, double* traj) const;
};

#endif // WHEEL_AIRBORNE_CONSTRAINT_HPP

