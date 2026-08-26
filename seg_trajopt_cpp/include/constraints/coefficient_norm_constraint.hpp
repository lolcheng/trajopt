#ifndef COEFFICIENT_NORM_CONSTRAINT_HPP
#define COEFFICIENT_NORM_CONSTRAINT_HPP

#include "constraint_base.hpp"

class CoefficientNormConstraint : public ConstraintBase {
public:
    CoefficientNormConstraint(int n_vars, double max_l2_norm);
    ~CoefficientNormConstraint() override = default;

    int getNumConstraints() const override { return 1; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_lb, double* g_ub, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_vars_;
    double max_l2_norm_sq_;
};

#endif
