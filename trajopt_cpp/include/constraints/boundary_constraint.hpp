#ifndef BOUNDARY_CONSTRAINT_HPP
#define BOUNDARY_CONSTRAINT_HPP

#include "../constraint_base.hpp"
#include "../rbf_terrain.hpp"
#include "../basis_function.hpp"
#include <memory>

/**
 * 边界约束
 * 约束起点和终点的位置和姿态
 */
class BoundaryConstraint : public ConstraintBase {
public:
    /**
     * 构造函数
     * @param n_coeff 基函数系数数量
     * @param x0, y0, psi0 起点位置和姿态
     * @param x1, y1, psi1 终点位置和姿态
     * @param clear_z 中心相对地面高度
     * @param terrain RBF地形（用于计算起点终点高度）
     * @param basis 基函数指针（如果为nullptr，则默认使用切比雪夫基）
     * @param wheels_nominal 轮子与中心的标称距离（默认0.5）
     */
    BoundaryConstraint(int n_coeff,
                      double x0, double y0, double psi0,
                      double x1, double y1, double psi1,
                      double clear_z,
                      const RBFTerrain& terrain,
                      std::shared_ptr<BasisFunction> basis = nullptr,
                      double wheels_nominal = 0.5);
    
    ~BoundaryConstraint() override;
    
    int getNumConstraints() const override { return 20; }  // 8个base约束 + 12个轮子约束（6个起点+6个终点）
    
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* jac_g, int offset = 0) override;
    void getBounds(double* g_lb, double* g_ub, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    double x0_, y0_, psi0_;
    double x1_, y1_, psi1_;
    double clear_z_;
    double wheels_nominal_;  // 轮子与中心的标称距离
    const RBFTerrain& terrain_;
    
    // 基函数指针
    std::shared_ptr<BasisFunction> basis_;
    
    // 预计算的基矩阵值（s=0和s=1）
    double* Phi0_;  // (1 × n_coeff)
    double* Phi1_;  // (1 × n_coeff)
    
    // 变量索引（假设变量按顺序排列）
    int idx_cbx_, idx_cby_, idx_cbz_, idx_cpsi_;
    int idx_clx_, idx_cly_, idx_clz_;  // 左轮
    int idx_crx_, idx_cry_, idx_crz_;  // 右轮
};

#endif // BOUNDARY_CONSTRAINT_HPP

