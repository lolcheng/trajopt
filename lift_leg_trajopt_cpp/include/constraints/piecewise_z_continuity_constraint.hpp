#ifndef PIECEWISE_Z_CONTINUITY_CONSTRAINT_HPP
#define PIECEWISE_Z_CONTINUITY_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "basis_function.hpp"
#include <memory>
#include <vector>

class PiecewiseZContinuityConstraint : public ConstraintBase {
public:
    PiecewiseZContinuityConstraint(int n_coeff, int n_seg_z,
                                   std::shared_ptr<BasisFunction> basis,
                                   double tol = 0.0);
    ~PiecewiseZContinuityConstraint() override = default;

    int getNumConstraints() const override { return 2 * std::max(0, n_seg_z_ - 1); }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_seg_z_;
    int n_coeff_z_;
    double tol_;
    int idx_clz_;
    int idx_crz_;
    std::shared_ptr<BasisFunction> basis_;
    std::vector<double> phi0_;
    std::vector<double> phi1_;
};

#endif
