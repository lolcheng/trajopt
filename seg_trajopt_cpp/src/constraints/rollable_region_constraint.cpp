#include "constraints/rollable_region_constraint.hpp"
#include "polynomial.hpp"
#include <cmath>

RollableRegionConstraint::RollableRegionConstraint(
    int n_coeff, int n_sample,
    std::shared_ptr<SegmentationReader> seg,
    double threshold,
    bool include_wheels,
    double buffer)
    : n_coeff_(n_coeff), n_sample_(n_sample), threshold_(threshold),
      buffer_(buffer), include_wheels_(include_wheels), seg_(std::move(seg)) {
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_crx_ = 6 * n_coeff_;
    idx_cry_ = 7 * n_coeff_;
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i)
        s_samples_[i] = (n_sample_ <= 1) ? 0.5 : static_cast<double>(i) / (n_sample_ - 1);
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void RollableRegionConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void RollableRegionConstraint::eval_g(const double* x, double* g) {
    if (!seg_) return;
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    std::vector<double> Xl(n_sample_), Yl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_);
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    for (int i = 0; i < n_sample_; ++i)
        g[i] = threshold_ - seg_->interp(Xb[i], Yb[i]);
    if (include_wheels_) {
        for (int i = 0; i < n_sample_; ++i) {
            g[n_sample_ + i] = threshold_ - seg_->interp(Xl[i], Yl[i]);
            g[2 * n_sample_ + i] = threshold_ - seg_->interp(Xr[i], Yr[i]);
        }
    }
}

void RollableRegionConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)offset;
    if (!seg_) return;
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    std::vector<double> Xl(n_sample_), Yl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_);
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        double ddx_b, ddy_b, ddx_l = 0.0, ddy_l = 0.0, ddx_r = 0.0, ddy_r = 0.0;
        seg_->interpGrad(Xb[i], Yb[i], ddx_b, ddy_b);
        if (include_wheels_) {
            seg_->interpGrad(Xl[i], Yl[i], ddx_l, ddy_l);
            seg_->interpGrad(Xr[i], Yr[i], ddx_r, ddy_r);
        }
        const double* phi = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j)
            values[idx++] = -ddx_b * phi[j];
        for (int j = 0; j < n_coeff_; ++j)
            values[idx++] = -ddy_b * phi[j];
        if (include_wheels_) {
            for (int j = 0; j < n_coeff_; ++j)
                values[idx++] = -ddx_l * phi[j];
            for (int j = 0; j < n_coeff_; ++j)
                values[idx++] = -ddy_l * phi[j];
            for (int j = 0; j < n_coeff_; ++j)
                values[idx++] = -ddx_r * phi[j];
            for (int j = 0; j < n_coeff_; ++j)
                values[idx++] = -ddy_r * phi[j];
        }
    }
}

void RollableRegionConstraint::getBounds(double* g_lb, double* g_ub, int offset) {
    (void)offset;
    const int n = include_wheels_ ? 3 * n_sample_ : n_sample_;
    for (int i = 0; i < n; ++i) {
        g_lb[i] = -1e20;
        g_ub[i] = buffer_;
    }
}

int RollableRegionConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol)
        return include_wheels_ ? (3 * n_sample_ * 2 * n_coeff_) : (n_sample_ * 2 * n_coeff_);
    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        int r0 = offset + i;
        for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = r0; jCol[idx] = idx_cbx_ + j; ++idx; }
        for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = r0; jCol[idx] = idx_cby_ + j; ++idx; }
        if (include_wheels_) {
            int r1 = offset + n_sample_ + i;
            int r2 = offset + 2 * n_sample_ + i;
            for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = r1; jCol[idx] = idx_clx_ + j; ++idx; }
            for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = r1; jCol[idx] = idx_cly_ + j; ++idx; }
            for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = r2; jCol[idx] = idx_crx_ + j; ++idx; }
            for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = r2; jCol[idx] = idx_cry_ + j; ++idx; }
        }
    }
    return idx;
}
