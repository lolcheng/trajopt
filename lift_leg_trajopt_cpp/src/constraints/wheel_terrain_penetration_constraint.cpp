#include "constraints/wheel_terrain_penetration_constraint.hpp"
#include "polynomial.hpp"
#include <algorithm>
#include <cmath>

WheelTerrainPenetrationConstraint::WheelTerrainPenetrationConstraint(
    int n_coeff, int n_sample, int n_seg_z, std::shared_ptr<RBFTerrain> terrain, double max_penetration)
    : n_coeff_(n_coeff), n_seg_z_(std::max(1, n_seg_z)), n_coeff_z_(n_coeff_ * n_seg_z_), n_sample_(n_sample),
      max_penetration_(max_penetration), terrain_(std::move(terrain)) {
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

void WheelTerrainPenetrationConstraint::eval_g(const double* x, double* g) {
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_), Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_clz_, n_coeff_z_, Zl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_crz_, n_coeff_z_, Zr.data());
    for (int i = 0; i < n_sample_; ++i) {
        g[i] = terrain_->height(Xl[i], Yl[i]) - Zl[i];
        g[n_sample_ + i] = terrain_->height(Xr[i], Yr[i]) - Zr[i];
    }
}

void WheelTerrainPenetrationConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)offset;
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Xr(n_sample_), Yr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());

    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        double hx, hy;
        terrain_->gradient(Xl[i], Yl[i], hx, hy);
        const double* phi = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = hx * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = hy * phi[j];
        const double* phi_z = Phi_z_.data() + i * n_coeff_z_;
        for (int j = 0; j < n_coeff_z_; ++j) values[idx++] = -phi_z[j];
    }
    for (int i = 0; i < n_sample_; ++i) {
        double hx, hy;
        terrain_->gradient(Xr[i], Yr[i], hx, hy);
        const double* phi = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = hx * phi[j];
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = hy * phi[j];
        const double* phi_z = Phi_z_.data() + i * n_coeff_z_;
        for (int j = 0; j < n_coeff_z_; ++j) values[idx++] = -phi_z[j];
    }
}

void WheelTerrainPenetrationConstraint::getBounds(double* g_l, double* g_u, int offset) {
    (void)offset;
    for (int i = 0; i < 2 * n_sample_; ++i) {
        g_l[i] = -1e20;
        g_u[i] = max_penetration_;
    }
}

int WheelTerrainPenetrationConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol) return 2 * n_sample_ * (2 * n_coeff_ + n_coeff_z_);
    int idx = 0;
    const int left_cols[] = {idx_clx_, idx_cly_, idx_clz_};
    const int right_cols[] = {idx_crx_, idx_cry_, idx_crz_};
    for (int i = 0; i < n_sample_; ++i) {
        for (int k = 0; k < 3; ++k) {
            const int width = (k == 2) ? n_coeff_z_ : n_coeff_;
            for (int j = 0; j < width; ++j) {
                iRow[idx] = offset + i;
                jCol[idx] = left_cols[k] + j;
                ++idx;
            }
        }
    }
    for (int i = 0; i < n_sample_; ++i) {
        for (int k = 0; k < 3; ++k) {
            const int width = (k == 2) ? n_coeff_z_ : n_coeff_;
            for (int j = 0; j < width; ++j) {
                iRow[idx] = offset + n_sample_ + i;
                jCol[idx] = right_cols[k] + j;
                ++idx;
            }
        }
    }
    return idx;
}
