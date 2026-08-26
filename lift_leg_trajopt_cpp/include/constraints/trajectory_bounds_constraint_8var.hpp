#ifndef TRAJECTORY_BOUNDS_CONSTRAINT_8VAR_HPP
#define TRAJECTORY_BOUNDS_CONSTRAINT_8VAR_HPP

#include "constraint_base.hpp"
#include <vector>

class TrajectoryBoundsConstraint8Var : public ConstraintBase {
public:
    TrajectoryBoundsConstraint8Var(int n_coeff, int n_sample,
                                   double x_min, double x_max,
                                   double y_min, double y_max);
    ~TrajectoryBoundsConstraint8Var() override = default;

    int getNumConstraints() const override { return 6 * n_sample_; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_, n_sample_;
    double x_min_, x_max_, y_min_, y_max_;
    int idx_cbx_, idx_cby_, idx_clx_, idx_cly_, idx_crx_, idx_cry_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
};

#endif
