#include "constraints/trajectory_bounds_constraint.hpp"

TrajectoryBoundsConstraint::TrajectoryBoundsConstraint(int n_coeff, int n_sample,
                                                       double x_min, double x_max,
                                                       double y_min, double y_max)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      x_min_(x_min), x_max_(x_max),
      y_min_(y_min), y_max_(y_max) {
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_crx_ = 7 * n_coeff_;
    idx_cry_ = 8 * n_coeff_;

    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }

    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

TrajectoryBoundsConstraint::~TrajectoryBoundsConstraint() = default;

void TrajectoryBoundsConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void TrajectoryBoundsConstraint::eval_g(const double* x, double* g) {
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

    // 约束顺序：Xb,Yb,Xl,Yl,Xr,Yr（每组 n_sample）
    for (int i = 0; i < n_sample_; ++i) {
        g[i] = Xb[i];
        g[n_sample_ + i] = Yb[i];
        g[2 * n_sample_ + i] = Xl[i];
        g[3 * n_sample_ + i] = Yl[i];
        g[4 * n_sample_ + i] = Xr[i];
        g[5 * n_sample_ + i] = Yr[i];
    }
}

void TrajectoryBoundsConstraint::eval_jac_g(const double* x, double* values, int offset) {
    (void)x;
    (void)offset;

    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        const double* phi = Phi_.data() + i * n_coeff_;

        // Xb
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        // Yb
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        // Xl
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        // Yl
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        // Xr
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
        // Yr
        for (int j = 0; j < n_coeff_; ++j) values[idx++] = phi[j];
    }
}

void TrajectoryBoundsConstraint::getBounds(double* g_l, double* g_u, int offset) {
    (void)offset;
    for (int i = 0; i < n_sample_; ++i) {
        g_l[i] = x_min_;
        g_u[i] = x_max_;

        g_l[n_sample_ + i] = y_min_;
        g_u[n_sample_ + i] = y_max_;

        g_l[2 * n_sample_ + i] = x_min_;
        g_u[2 * n_sample_ + i] = x_max_;

        g_l[3 * n_sample_ + i] = y_min_;
        g_u[3 * n_sample_ + i] = y_max_;

        g_l[4 * n_sample_ + i] = x_min_;
        g_u[4 * n_sample_ + i] = x_max_;

        g_l[5 * n_sample_ + i] = y_min_;
        g_u[5 * n_sample_ + i] = y_max_;
    }
}

int TrajectoryBoundsConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        return 6 * n_sample_ * n_coeff_;
    }

    int idx = 0;
    for (int i = 0; i < n_sample_; ++i) {
        int row_xb = offset + i;
        int row_yb = offset + n_sample_ + i;
        int row_xl = offset + 2 * n_sample_ + i;
        int row_yl = offset + 3 * n_sample_ + i;
        int row_xr = offset + 4 * n_sample_ + i;
        int row_yr = offset + 5 * n_sample_ + i;

        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row_xb; jCol[idx] = idx_cbx_ + j; ++idx;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row_yb; jCol[idx] = idx_cby_ + j; ++idx;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row_xl; jCol[idx] = idx_clx_ + j; ++idx;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row_yl; jCol[idx] = idx_cly_ + j; ++idx;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row_xr; jCol[idx] = idx_crx_ + j; ++idx;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row_yr; jCol[idx] = idx_cry_ + j; ++idx;
        }
    }
    return idx;
}


