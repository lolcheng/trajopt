#ifndef NORMAL_FORCE_UPPER_BOUND_CONSTRAINT_HPP
#define NORMAL_FORCE_UPPER_BOUND_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "rbf_terrain.hpp"
#include "polynomial.hpp"
#include <vector>
#include <memory>

/**
 * 法向力上限约束
 * fLn ≤ fN_max = MASS × (GRAVITY + 1.0)
 * fRn ≤ fN_max
 * 其中: fLn = fL · nL, fRn = fR · nR (法向分量)
 * 
 * 共 2×n_sample 个不等式约束
 */
class NormalForceUpperBoundConstraint : public ConstraintBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param terrain RBF地形（用于计算法向量）
     * @param mass 质量（默认 1.0）
     * @param gravity 重力加速度（默认 9.81）
     * @param max_accel 最大加速度（默认 1.0）
     */
    NormalForceUpperBoundConstraint(int n_coeff, int n_sample,
                                    std::shared_ptr<RBFTerrain> terrain,
                                    double mass = 1.0, double gravity = 9.81,
                                    double max_accel = 1.0);
    
    virtual ~NormalForceUpperBoundConstraint() override;
    
    int getNumConstraints() const override;
    
    void eval_g(const double* x, double* g) override;
    
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;
    double fN_max_;  // MASS × (GRAVITY + max_accel)
    std::shared_ptr<RBFTerrain> terrain_;
    
    // 变量索引
    int idx_fLx_, idx_fLy_, idx_fLz_;
    int idx_fRx_, idx_fRy_, idx_fRz_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：计算法向力
    double computeNormalForce(double fx, double fy, double fz,
                              double nx, double ny, double nz) const;
};

#endif // NORMAL_FORCE_UPPER_BOUND_CONSTRAINT_HPP

