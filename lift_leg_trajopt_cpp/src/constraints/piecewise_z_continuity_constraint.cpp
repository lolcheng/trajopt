#include "constraints/piecewise_z_continuity_constraint.hpp"
#include "chebyshev_basis.hpp"
#include <algorithm>
#include <cmath>

PiecewiseZContinuityConstraint::PiecewiseZContinuityConstraint(
    int n_coeff, int n_seg_z, std::shared_ptr<BasisFunction> basis, double tol)
    : n_coeff_(n_coeff),
      n_seg_z_(std::max(1, n_seg_z)),
      n_coeff_z_(n_coeff_ * n_seg_z_),
      tol_(std::abs(tol)),
      idx_clz_(6 * n_coeff_),
      idx_crz_(0) {
    const int idx_crx = idx_clz_ + n_coeff_z_;
    const int idx_cry = idx_crx + n_coeff_;
    idx_crz_ = idx_cry + n_coeff_;
    basis_ = basis ? basis : std::make_shared<ChebyshevBasis>();
    phi0_.resize(n_coeff_);
    phi1_.resize(n_coeff_);
    const double s0 = 0.0;
    const double s1 = 1.0;
    basis_->computeBasisMatrix(&s0, 1, n_coeff_, phi0_.data());
    basis_->computeBasisMatrix(&s1, 1, n_coeff_, phi1_.data());
}

void PiecewiseZContinuityConstraint::eval_g(const double* x, double* g) {
    const int n_link = std::max(0, n_seg_z_ - 1);
    const double* clz = x + idx_clz_;
    const double* crz = x + idx_crz_;
    for (int seg = 0; seg < n_link; ++seg) {
        const double* clz_curr = clz + seg * n_coeff_;
        const double* clz_next = clz + (seg + 1) * n_coeff_;
        const double* crz_curr = crz + seg * n_coeff_;
        const double* crz_next = crz + (seg + 1) * n_coeff_;
        const double zl_end = basis_->evalPoint(1.0, clz_curr, n_coeff_);
        const double zl_next_start = basis_->evalPoint(0.0, clz_next, n_coeff_);
        const double zr_end = basis_->evalPoint(1.0, crz_curr, n_coeff_);
        const double zr_next_start = basis_->evalPoint(0.0, crz_next, n_coeff_);
        g[seg] = zl_end - zl_next_start;
        g[n_link + seg] = zr_end - zr_next_start;
    }
}

void PiecewiseZContinuityConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)x;
    (void)offset;
    const int n_link = std::max(0, n_seg_z_ - 1);
    int idx = 0;
    for (int seg = 0; seg < n_link; ++seg) {
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi1_[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = -phi0_[j];
    }
    for (int seg = 0; seg < n_link; ++seg) {
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi1_[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = -phi0_[j];
    }
}

void PiecewiseZContinuityConstraint::getBounds(double* g_l, double* g_u, int offset) {
    (void)offset;
    const int m = getNumConstraints();
    for (int i = 0; i < m; ++i) {
        g_l[i] = -tol_;
        g_u[i] = tol_;
    }
}

int PiecewiseZContinuityConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    const int n_link = std::max(0, n_seg_z_ - 1);
    if (!iRow || !jCol) return 4 * n_link * n_coeff_;
    int idx = 0;
    for (int seg = 0; seg < n_link; ++seg) {
        const int row = offset + seg;
        const int col_curr = idx_clz_ + seg * n_coeff_;
        const int col_next = idx_clz_ + (seg + 1) * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row;
            jCol[idx++] = col_curr + j;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row;
            jCol[idx++] = col_next + j;
        }
    }
    for (int seg = 0; seg < n_link; ++seg) {
        const int row = offset + n_link + seg;
        const int col_curr = idx_crz_ + seg * n_coeff_;
        const int col_next = idx_crz_ + (seg + 1) * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row;
            jCol[idx++] = col_curr + j;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row;
            jCol[idx++] = col_next + j;
        }
    }
    return idx;
}
