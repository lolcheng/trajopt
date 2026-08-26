#include "constraints/trajectory_bounds_constraint_8var.hpp"
#include "polynomial.hpp"
#include <vector>

TrajectoryBoundsConstraint8Var::TrajectoryBoundsConstraint8Var(int n_coeff, int n_sample,
                                                               double x_min, double x_max,
                                                               double y_min, double y_max)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max) {
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_crx_ = 6 * n_coeff_;
    idx_cry_ = 7 * n_coeff_;
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void TrajectoryBoundsConstraint8Var::eval_g(const double* x, double* g) {
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Xl(n_sample_), Yl(n_sample_), Xr(n_sample_), Yr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cbx_, n_coeff_, Xb.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cby_, n_coeff_, Yb.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    for (int i = 0; i < n_sample_; ++i) {
        g[i] = Xb[i];
        g[n_sample_ + i] = Yb[i];
        g[2 * n_sample_ + i] = Xl[i];
        g[3 * n_sample_ + i] = Yl[i];
        g[4 * n_sample_ + i] = Xr[i];
        g[5 * n_sample_ + i] = Yr[i];
    }
}

void TrajectoryBoundsConstraint8Var::eval_jac_g(const double* x, double* values, int offset) {
    (void)x;
    (void)offset;
    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        const double* phi = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
    }
}

void TrajectoryBoundsConstraint8Var::getBounds(double* g_l, double* g_u, int offset) {
    (void)offset;
    for (int i = 0; i < n_sample_; ++i) {
        g_l[i] = x_min_; g_u[i] = x_max_;
        g_l[n_sample_ + i] = y_min_; g_u[n_sample_ + i] = y_max_;
        g_l[2 * n_sample_ + i] = x_min_; g_u[2 * n_sample_ + i] = x_max_;
        g_l[3 * n_sample_ + i] = y_min_; g_u[3 * n_sample_ + i] = y_max_;
        g_l[4 * n_sample_ + i] = x_min_; g_u[4 * n_sample_ + i] = x_max_;
        g_l[5 * n_sample_ + i] = y_min_; g_u[5 * n_sample_ + i] = y_max_;
    }
}

int TrajectoryBoundsConstraint8Var::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol) return 6 * n_sample_ * n_coeff_;
    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        int rows[] = {offset + i, offset + n_sample_ + i, offset + 2 * n_sample_ + i,
                      offset + 3 * n_sample_ + i, offset + 4 * n_sample_ + i, offset + 5 * n_sample_ + i};
        int cols[] = {idx_cbx_, idx_cby_, idx_clx_, idx_cly_, idx_crx_, idx_cry_};
        for (int k = 0; k < 6; ++k) {
            for (int j = 0; j < n_coeff_; ++j) {
                iRow[idx] = rows[k];
                jCol[idx] = cols[k] + j;
                ++idx;
            }
        }
    }
    return idx;
}
