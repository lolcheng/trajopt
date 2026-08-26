#include "constraints/coefficient_norm_constraint.hpp"
#include <algorithm>

CoefficientNormConstraint::CoefficientNormConstraint(int n_vars, double max_l2_norm)
    : n_vars_(n_vars),
      max_l2_norm_sq_(std::max(1e-12, max_l2_norm * max_l2_norm)) {}

void CoefficientNormConstraint::eval_g(const double* x, double* g) {
    double norm_sq = 0.0;
    for (int i = 0; i < n_vars_; ++i) {
        norm_sq += x[i] * x[i];
    }
    g[0] = norm_sq / static_cast<double>(n_vars_);
}

void CoefficientNormConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)offset;
    const double scale = 2.0 / static_cast<double>(n_vars_);
    for (int i = 0; i < n_vars_; ++i) {
        values[i] = scale * x[i];
    }
}

void CoefficientNormConstraint::getBounds(double* g_lb, double* g_ub, int offset) {
    (void)offset;
    g_lb[0] = -1e20;
    g_ub[0] = max_l2_norm_sq_;
}

int CoefficientNormConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol) return n_vars_;
    for (int i = 0; i < n_vars_; ++i) {
        iRow[i] = offset;
        jCol[i] = i;
    }
    return n_vars_;
}
