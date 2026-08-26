#include "objectives/swing_wheel_lift_objective.hpp"
#include "polynomial.hpp"
#include <algorithm>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SwingWheelLiftObjective::SwingWheelLiftObjective(
    int n_coeff, int n_sample, int n_seg_z, std::shared_ptr<RBFTerrain> terrain,
    bool swing_left, double lift_peak, double stance_clearance_target, double weight)
    : n_coeff_(n_coeff), n_seg_z_(std::max(1, n_seg_z)), n_coeff_z_(n_coeff_ * n_seg_z_), n_sample_(n_sample), swing_left_(swing_left),
      lift_peak_(lift_peak), stance_clearance_target_(stance_clearance_target),
      weight_(weight), enabled_(true), terrain_(std::move(terrain)) {
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

double SwingWheelLiftObjective::profile(double s) const {
    const double ss = std::sin(M_PI * s);
    return lift_peak_ * ss * ss;
}

double SwingWheelLiftObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_), Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_clz_, n_coeff_z_, Zl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_crz_, n_coeff_z_, Zr.data());

    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        const double hl = terrain_->height(Xl[i], Yl[i]);
        const double hr = terrain_->height(Xr[i], Yr[i]);
        const double cl = Zl[i] - hl;
        const double cr = Zr[i] - hr;
        const double cl_ref = swing_left_ ? profile(s_samples_[i]) : stance_clearance_target_;
        const double cr_ref = swing_left_ ? stance_clearance_target_ : profile(s_samples_[i]);
        const double el = cl - cl_ref;
        const double er = cr - cr_ref;
        obj += el * el + er * er;
    }
    return weight_ * obj;
}

void SwingWheelLiftObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_), Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_clx_, n_coeff_, Xl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cly_, n_coeff_, Yl.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_clz_, n_coeff_z_, Zl.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_crx_, n_coeff_, Xr.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cry_, n_coeff_, Yr.data());
    Polynomial::eval(Phi_z_.data(), n_sample_, x + idx_crz_, n_coeff_z_, Zr.data());

    for (int i = 0; i < n_sample_; ++i) {
        const double* phi = Phi_.data() + i * n_coeff_;
        const double* phi_z = Phi_z_.data() + i * n_coeff_z_;
        double hxl, hyl, hxr, hyr;
        terrain_->gradient(Xl[i], Yl[i], hxl, hyl);
        terrain_->gradient(Xr[i], Yr[i], hxr, hyr);
        const double hl = terrain_->height(Xl[i], Yl[i]);
        const double hr = terrain_->height(Xr[i], Yr[i]);
        const double cl_ref = swing_left_ ? profile(s_samples_[i]) : stance_clearance_target_;
        const double cr_ref = swing_left_ ? stance_clearance_target_ : profile(s_samples_[i]);
        const double el = (Zl[i] - hl) - cl_ref;
        const double er = (Zr[i] - hr) - cr_ref;
        const double fl = 2.0 * weight_ * el;
        const double fr = 2.0 * weight_ * er;
        for (int j = 0; j < n_coeff_; ++j) {
            grad_f[idx_clx_ + j] += fl * (-hxl) * phi[j];
            grad_f[idx_cly_ + j] += fl * (-hyl) * phi[j];
            grad_f[idx_crx_ + j] += fr * (-hxr) * phi[j];
            grad_f[idx_cry_ + j] += fr * (-hyr) * phi[j];
        }
        for (int j = 0; j < n_coeff_z_; ++j) grad_f[idx_clz_ + j] += fl * phi_z[j];
        for (int j = 0; j < n_coeff_z_; ++j) grad_f[idx_crz_ + j] += fr * phi_z[j];
    }
}
