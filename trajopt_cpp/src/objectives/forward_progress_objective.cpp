#include "objectives/forward_progress_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

ForwardProgressObjective::ForwardProgressObjective(int n_coeff, int n_sample,
                                                   double xb0, double yb0,
                                                   double xb1, double yb1,
                                                   double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      weight_(weight < 0 ? 50.0 / n_sample : weight),
      enabled_(true),
      dir_x_(0.0), dir_y_(0.0) {
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;

    double dx = xb1 - xb0;
    double dy = yb1 - yb0;
    double norm = std::sqrt(dx * dx + dy * dy);
    if (norm > 1e-12) {
        dir_x_ = dx / norm;
        dir_y_ = dy / norm;
    } else {
        dir_x_ = 0.0;
        dir_y_ = -1.0;
    }

    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }

    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void ForwardProgressObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double ForwardProgressObjective::smoothPos(double z, double eps) {
    return 0.5 * (z + std::sqrt(z * z + eps));
}

double ForwardProgressObjective::smoothPosGrad(double z, double eps) {
    return 0.5 * (1.0 + z / std::sqrt(z * z + eps));
}

double ForwardProgressObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;

    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;

    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());

    double obj = 0.0;
    for (int i = 0; i < n_sample_ - 1; ++i) {
        double step_x = Xb[i + 1] - Xb[i];
        double step_y = Yb[i + 1] - Yb[i];
        double progress = step_x * dir_x_ + step_y * dir_y_;
        double backtrack = smoothPos(-progress);
        obj += backtrack * backtrack;
    }

    return weight_ * obj;
}

void ForwardProgressObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;

    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    double* grad_cbx = grad_f + idx_cbx_;
    double* grad_cby = grad_f + idx_cby_;

    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());

    for (int i = 0; i < n_sample_ - 1; ++i) {
        double step_x = Xb[i + 1] - Xb[i];
        double step_y = Yb[i + 1] - Yb[i];
        double progress = step_x * dir_x_ + step_y * dir_y_;
        double z = -progress;
        double backtrack = smoothPos(z);
        double dbacktrack_dz = smoothPosGrad(z);

        // d(backtrack^2)/dprogress = 2*backtrack*dbacktrack/dz * dz/dprogress
        // z = -progress => dz/dprogress = -1
        double dterm_dprogress = -2.0 * backtrack * dbacktrack_dz;
        double factor = weight_ * dterm_dprogress;

        const double* phi_i = Phi_.data() + i * n_coeff_;
        const double* phi_ip1 = Phi_.data() + (i + 1) * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            double dstepx = (phi_ip1[j] - phi_i[j]);
            double dstepy = (phi_ip1[j] - phi_i[j]);
            double dprogress_dcbx = dir_x_ * dstepx;
            double dprogress_dcby = dir_y_ * dstepy;
            grad_cbx[j] += factor * dprogress_dcbx;
            grad_cby[j] += factor * dprogress_dcby;
        }
    }
}


