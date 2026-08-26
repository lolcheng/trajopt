#include "constraints/wheels_distance_terrain_constraint.hpp"
#include "polynomial.hpp"
#include <vector>

WheelsDistanceTerrainConstraint::WheelsDistanceTerrainConstraint(
    int n_coeff, int n_sample, std::shared_ptr<RBFTerrain> terrain,
    double min_dist, double max_dist)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      min_dist_sq_(min_dist * min_dist),
      max_dist_sq_(max_dist * max_dist),
      terrain_(std::move(terrain)) {
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_crx_ = 6 * n_coeff_;
    idx_cry_ = 7 * n_coeff_;
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void WheelsDistanceTerrainConstraint::eval_g(const double* x, double* g) {
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Xr(n_sample_), Yr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    for (int i = 0; i < n_sample_; ++i) {
        double zl = terrain_->height(Xl[i], Yl[i]);
        double zr = terrain_->height(Xr[i], Yr[i]);
        double dx = Xr[i] - Xl[i], dy = Yr[i] - Yl[i], dz = zr - zl;
        g[i] = dx * dx + dy * dy + dz * dz;
    }
}

void WheelsDistanceTerrainConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)offset;
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Xr(n_sample_), Yr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        const double* phi = Phi_.data() + i * n_coeff_;
        double ghx_l, ghy_l, ghx_r, ghy_r;
        terrain_->gradient(Xl[i], Yl[i], ghx_l, ghy_l);
        terrain_->gradient(Xr[i], Yr[i], ghx_r, ghy_r);
        double zl = terrain_->height(Xl[i], Yl[i]);
        double zr = terrain_->height(Xr[i], Yr[i]);
        double dx = Xr[i] - Xl[i], dy = Yr[i] - Yl[i], dz = zr - zl;
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = (-2.0 * dx - 2.0 * dz * ghx_l) * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = (-2.0 * dy - 2.0 * dz * ghy_l) * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = (2.0 * dx + 2.0 * dz * ghx_r) * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = (2.0 * dy + 2.0 * dz * ghy_r) * phi[j];
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
    if (!iRow || !jCol) return n_sample_ * 4 * n_coeff_;
    int idx = 0;
    int cols[] = {idx_clx_, idx_cly_, idx_crx_, idx_cry_};
    for (int i = 0; i < n_sample_; ++i) {
        for (int k = 0; k < 4; ++k) {
            for (int j = 0; j < n_coeff_; ++j) {
                iRow[idx] = offset + i;
                jCol[idx] = cols[k] + j;
                ++idx;
            }
        }
    }
    return idx;
}
