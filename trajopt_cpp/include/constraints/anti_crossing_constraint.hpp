#ifndef ANTI_CROSSING_CONSTRAINT_HPP
#define ANTI_CROSSING_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 防交叉约束
 * 在base坐标系下: yL_base ≥ yR_base
 * 其中:
 *   u_lat = [-sin(psi), cos(psi)] (横向单位向量)
 *   yL_base = (Xl - Xb) · u_lat_x + (Yl - Yb) · u_lat_y
 *   yR_base = (Xr - Xb) · u_lat_x + (Yr - Yb) · u_lat_y
 * 
 * 共 n_sample 个不等式约束
 */
class AntiCrossingConstraint : public ConstraintBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     */
    AntiCrossingConstraint(int n_coeff, int n_sample);
    
    virtual ~AntiCrossingConstraint() override;
    
    int getNumConstraints() const override;
    
    void eval_g(const double* x, double* g) override;
    
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;
    
    // 变量索引
    int idx_cbx_, idx_cby_;
    int idx_clx_, idx_cly_;
    int idx_crx_, idx_cry_;
    int idx_cpsi_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
};

#endif // ANTI_CROSSING_CONSTRAINT_HPP

