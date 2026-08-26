#include "objectives/base_height_stability_objective.hpp"
#include "polynomial.hpp"

BaseHeightStabilityObjective::BaseHeightStabilityObjective(
    int n_coeff, int n_sample, std::shared_ptr<RBFTerrain> terrain,
    double clear_z_nominal, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      clear_z_nominal_(clear_z_nominal),
      weight_(weight < 0 ? 1.0 / n_sample : weight),
      enabled_(true), terrain_(std::move(terrain)) {
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

double BaseHeightStabilityObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Zb(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cbx_, n_coeff_, Xb.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cby_, n_coeff_, Yb.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cbz_, n_coeff_, Zb.data());

    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        double h = terrain_->height(Xb[i], Yb[i]);
        double e = (Zb[i] - h) - clear_z_nominal_;
        obj += e * e;
    }
    return weight_ * obj;
}

void BaseHeightStabilityObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Zb(n_sample_);
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cbx_, n_coeff_, Xb.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cby_, n_coeff_, Yb.data());
    Polynomial::eval(Phi_.data(), n_sample_, x + idx_cbz_, n_coeff_, Zb.data());

    for (int i = 0; i < n_sample_; ++i) {
        const double* phi = Phi_.data() + i * n_coeff_;
        double hx, hy;
        terrain_->gradient(Xb[i], Yb[i], hx, hy);
        double h = terrain_->height(Xb[i], Yb[i]);
        double e = (Zb[i] - h) - clear_z_nominal_;
        double fac = 2.0 * weight_ * e;
        for (int j = 0; j < n_coeff_; ++j) {
            grad_f[idx_cbx_ + j] += fac * (-hx) * phi[j];
            grad_f[idx_cby_ + j] += fac * (-hy) * phi[j];
            grad_f[idx_cbz_ + j] += fac * phi[j];
        }
    }
}
