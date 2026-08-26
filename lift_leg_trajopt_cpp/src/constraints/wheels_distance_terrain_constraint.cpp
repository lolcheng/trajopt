#include "constraints/wheels_distance_terrain_constraint.hpp"
#include "polynomial.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

WheelsDistanceTerrainConstraint::WheelsDistanceTerrainConstraint(
    int n_coeff, int n_sample, int n_seg_z, double min_dist, double max_dist)
    : n_coeff_(n_coeff), n_seg_z_(std::max(1, n_seg_z)), n_coeff_z_(n_coeff_ * n_seg_z_), n_sample_(n_sample),
      min_dist_sq_(min_dist * min_dist),
      max_dist_sq_(max_dist * max_dist) {
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = idx_clz_ + n_coeff_z_;
    idx_cry_ = idx_crx_ + n_coeff_;
    idx_crz_ = idx_cry_ + n_coeff_;
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
    Phi_z_.assign(n_sample_ * n_coeff_z_, 0.0);
    for (int i = 0; i < n_sample_; ++i) {
        const double s = s_samples_[i];
        const int seg = std::min(static_cast<int>(std::floor(s * n_seg_z_)), n_seg_z_ - 1);
        const double local = std::clamp(s * n_seg_z_ - seg, 0.0, 1.0);
        std::vector<double> phi_local(n_coeff_);
        Polynomial::vandermonde(&local, 1, n_coeff_, phi_local.data());
        const int row = i * n_coeff_z_ + seg * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) Phi_z_[row + j] = phi_local[j];
    }
}

void WheelsDistanceTerrainConstraint::eval_g(const double* x, double* g) {
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_), Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_clz_, n_coeff_z_, Zl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_crz_, n_coeff_z_, Zr.data());
    for (int i = 0; i < n_sample_; ++i) {
        double dx = Xr[i] - Xl[i], dy = Yr[i] - Yl[i], dz = Zr[i] - Zl[i];
        g[i] = dx * dx + dy * dy + dz * dz;
    }
}

void WheelsDistanceTerrainConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)offset;
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_), Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_clz_, n_coeff_z_, Zl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_crz_, n_coeff_z_, Zr.data());
    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        const double* phi = Phi_.data() + i * n_coeff_;
        const double* phi_z = Phi_z_.data() + i * n_coeff_z_;
        double dx = Xr[i] - Xl[i], dy = Yr[i] - Yl[i], dz = Zr[i] - Zl[i];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = -2.0 * dx * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = -2.0 * dy * phi[j];
        for (int j = 0; j < n_coeff_z_; ++j) values[idx++] = -2.0 * dz * phi_z[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = 2.0 * dx * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = 2.0 * dy * phi[j];
        for (int j = 0; j < n_coeff_z_; ++j) values[idx++] = 2.0 * dz * phi_z[j];
    }
}

void WheelsDistanceTerrainConstraint::getBounds(double* g_l, double* g_u, int offset) {
    (void)offset;
    for (int i = 0; i < n_sample_; ++i) {
        g_l[i] = min_dist_sq_;
        g_u[i] = max_dist_sq_;
    }
}

int WheelsDistanceTerrainConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol) return n_sample_ * (4 * n_coeff_ + 2 * n_coeff_z_);
    int idx = 0;
    int cols[] = {idx_clx_, idx_cly_, idx_clz_, idx_crx_, idx_cry_, idx_crz_};
    for (int i = 0; i < n_sample_; ++i) {
        for (int k = 0; k < 6; ++k) {
            const int width = (k == 2 || k == 5) ? n_coeff_z_ : n_coeff_;
            for (int j = 0; j < width; ++j) {
                iRow[idx] = offset + i;
                jCol[idx] = cols[k] + j;
                ++idx;
            }
        }
    }
    return idx;
}
