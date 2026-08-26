#ifndef TRAJECTORY_BOUNDS_CONSTRAINT_HPP
#define TRAJECTORY_BOUNDS_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 轨迹平面边界约束（硬约束）
 * 将 base / left / right 三条轨迹的 x,y 坐标严格限制在指定范围内：
 *   x_min <= x <= x_max
 *   y_min <= y <= y_max
 */
class TrajectoryBoundsConstraint : public ConstraintBase {
public:
    TrajectoryBoundsConstraint(int n_coeff, int n_sample,
                               double x_min, double x_max,
                               double y_min, double y_max);

    ~TrajectoryBoundsConstraint() override;

    // 6组轨迹量：Xb,Yb,Xl,Yl,Xr,Yr，每组 n_sample 个约束
    int getNumConstraints() const override { return 6 * n_sample_; }

    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;

    double x_min_;
    double x_max_;
    double y_min_;
    double y_max_;

    // 变量索引
    int idx_cbx_, idx_cby_;
    int idx_clx_, idx_cly_;
    int idx_crx_, idx_cry_;

    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample x n_coeff)

    void evalTrajectory(const double* coeff, double* traj) const;
};

#endif // TRAJECTORY_BOUNDS_CONSTRAINT_HPP


