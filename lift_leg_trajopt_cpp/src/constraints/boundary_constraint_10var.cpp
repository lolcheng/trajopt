#include "constraints/boundary_constraint_10var.hpp"
#include "chebyshev_basis.hpp"
#include <algorithm>
#include <cmath>

BoundaryConstraint10Var::BoundaryConstraint10Var(int n_coeff, int n_seg_z,
                                                 double x0, double y0, double psi0,
                                                 double x1, double y1, double psi1,
                                                 double clear_z,
                                                 const RBFTerrain& terrain,
                                                 std::shared_ptr<BasisFunction> basis,
                                                 double wheels_nominal,
                                                 double base_xy_tol,
                                                 double base_z_tol,
                                                 double base_yaw_tol,
                                                 double wheel_xy_tol,
                                                 double wheel_z_tol,
                                                 double base_z_min,
                                                 double base_z_max,
                                                 double wheel_z_min,
                                                 double wheel_z_max)
    : n_coeff_(n_coeff),
      n_seg_z_(std::max(1, n_seg_z)),
      n_coeff_z_(n_coeff_ * n_seg_z_),
      x0_(x0), y0_(y0), psi0_(psi0),
      x1_(x1), y1_(y1), psi1_(psi1),
      clear_z_(clear_z),
      wheels_nominal_(wheels_nominal),
      base_xy_tol_(std::abs(base_xy_tol)),
      base_z_tol_(std::abs(base_z_tol)),
      base_yaw_tol_(std::abs(base_yaw_tol)),
      wheel_xy_tol_(std::abs(wheel_xy_tol)),
      wheel_z_tol_(std::abs(wheel_z_tol)),
      base_z_min_(std::min(base_z_min, base_z_max)),
      base_z_max_(std::max(base_z_min, base_z_max)),
      wheel_z_min_(std::min(wheel_z_min, wheel_z_max)),
      wheel_z_max_(std::max(wheel_z_min, wheel_z_max)),
      terrain_(terrain),
      Phi0_(new double[n_coeff_]),
      Phi1_(new double[n_coeff_]),
      Phi0z_(new double[n_coeff_z_]),
      Phi1z_(new double[n_coeff_z_]) {
    basis_ = basis ? basis : std::make_shared<ChebyshevBasis>();
    double s0 = 0.0, s1 = 1.0;
    basis_->computeBasisMatrix(&s0, 1, n_coeff_, Phi0_);
    basis_->computeBasisMatrix(&s1, 1, n_coeff_, Phi1_);
    std::fill(Phi0z_, Phi0z_ + n_coeff_z_, 0.0);
    std::fill(Phi1z_, Phi1z_ + n_coeff_z_, 0.0);
    basis_->computeBasisMatrix(&s0, 1, n_coeff_, Phi0z_);
    basis_->computeBasisMatrix(&s1, 1, n_coeff_, Phi1z_ + (n_seg_z_ - 1) * n_coeff_);
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = idx_clz_ + n_coeff_z_;
    idx_cry_ = idx_crx_ + n_coeff_;
    idx_crz_ = idx_cry_ + n_coeff_;
}

BoundaryConstraint10Var::~BoundaryConstraint10Var() {
    delete[] Phi0_;
    delete[] Phi1_;
    delete[] Phi0z_;
    delete[] Phi1z_;
}

void BoundaryConstraint10Var::eval_g(const double* x, double* g) {
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    const double* cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;

    double xb0 = basis_->evalPoint(0.0, cbx, n_coeff_);
    double yb0 = basis_->evalPoint(0.0, cby, n_coeff_);
    double zb0 = basis_->evalPoint(0.0, cbz, n_coeff_);
    double yaw0 = basis_->evalPoint(0.0, cpsi, n_coeff_);
    double xb1 = basis_->evalPoint(1.0, cbx, n_coeff_);
    double yb1 = basis_->evalPoint(1.0, cby, n_coeff_);
    double zb1 = basis_->evalPoint(1.0, cbz, n_coeff_);
    double yaw1 = basis_->evalPoint(1.0, cpsi, n_coeff_);
    double xl0 = basis_->evalPoint(0.0, clx, n_coeff_);
    double yl0 = basis_->evalPoint(0.0, cly, n_coeff_);
    double zl0 = basis_->evalPoint(0.0, clz, n_coeff_);
    double xr0 = basis_->evalPoint(0.0, crx, n_coeff_);
    double yr0 = basis_->evalPoint(0.0, cry, n_coeff_);
    double zr0 = basis_->evalPoint(0.0, crz, n_coeff_);
    double xl1 = basis_->evalPoint(1.0, clx, n_coeff_);
    double yl1 = basis_->evalPoint(1.0, cly, n_coeff_);
    const double* clz_last = clz + (n_seg_z_ - 1) * n_coeff_;
    const double* crz_last = crz + (n_seg_z_ - 1) * n_coeff_;
    double zl1 = basis_->evalPoint(1.0, clz_last, n_coeff_);
    double xr1 = basis_->evalPoint(1.0, crx, n_coeff_);
    double yr1 = basis_->evalPoint(1.0, cry, n_coeff_);
    double zr1 = basis_->evalPoint(1.0, crz_last, n_coeff_);

    double z0_target = terrain_.height(x0_, y0_) + clear_z_;
    double z1_target = terrain_.height(x1_, y1_) + clear_z_;
    double ux0 = -std::sin(yaw0), uy0 = std::cos(yaw0);
    double ux1 = -std::sin(yaw1), uy1 = std::cos(yaw1);
    double xl0_target = xb0 + 0.5 * wheels_nominal_ * ux0;
    double yl0_target = yb0 + 0.5 * wheels_nominal_ * uy0;
    double xr0_target = xb0 - 0.5 * wheels_nominal_ * ux0;
    double yr0_target = yb0 - 0.5 * wheels_nominal_ * uy0;
    double xl1_target = xb1 + 0.5 * wheels_nominal_ * ux1;
    double yl1_target = yb1 + 0.5 * wheels_nominal_ * uy1;
    double xr1_target = xb1 - 0.5 * wheels_nominal_ * ux1;
    double yr1_target = yb1 - 0.5 * wheels_nominal_ * uy1;

    g[0] = xb0 - x0_;
    g[1] = yb0 - y0_;
    g[2] = zb0 - z0_target;
    g[3] = yaw0 - psi0_;
    g[4] = xb1 - x1_;
    g[5] = yb1 - y1_;
    g[6] = zb1 - z1_target;
    g[7] = yaw1 - psi1_;

    g[8] = xl0 - xl0_target;
    g[9] = yl0 - yl0_target;
    g[10] = zl0 - terrain_.height(xl0, yl0);
    g[11] = xr0 - xr0_target;
    g[12] = yr0 - yr0_target;
    g[13] = zr0 - terrain_.height(xr0, yr0);

    g[14] = xl1 - xl1_target;
    g[15] = yl1 - yl1_target;
    g[16] = zl1 - terrain_.height(xl1, yl1);
    g[17] = xr1 - xr1_target;
    g[18] = yr1 - yr1_target;
    g[19] = zr1 - terrain_.height(xr1, yr1);

    // Absolute z bounds at endpoints.
    g[20] = zb0;
    g[21] = zb1;
    g[22] = zl0;
    g[23] = zr0;
    g[24] = zl1;
    g[25] = zr1;
}

void BoundaryConstraint10Var::eval_jac_g(const double* x, double* jac_g, int offset) {
    (void)offset;
    const double* cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    double yaw0 = basis_->evalPoint(0.0, cpsi, n_coeff_);
    double yaw1 = basis_->evalPoint(1.0, cpsi, n_coeff_);
    double xl0 = basis_->evalPoint(0.0, clx, n_coeff_);
    double yl0 = basis_->evalPoint(0.0, cly, n_coeff_);
    double xr0 = basis_->evalPoint(0.0, crx, n_coeff_);
    double yr0 = basis_->evalPoint(0.0, cry, n_coeff_);
    double xl1 = basis_->evalPoint(1.0, clx, n_coeff_);
    double yl1 = basis_->evalPoint(1.0, cly, n_coeff_);
    double xr1 = basis_->evalPoint(1.0, crx, n_coeff_);
    double yr1 = basis_->evalPoint(1.0, cry, n_coeff_);
    double hxl0, hyl0, hxr0, hyr0, hxl1, hyl1, hxr1, hyr1;
    terrain_.gradient(xl0, yl0, hxl0, hyl0);
    terrain_.gradient(xr0, yr0, hxr0, hyr0);
    terrain_.gradient(xl1, yl1, hxl1, hyl1);
    terrain_.gradient(xr1, yr1, hxr1, hyr1);

    double dux0 = -std::cos(yaw0), duy0 = -std::sin(yaw0);
    double dux1 = -std::cos(yaw1), duy1 = -std::sin(yaw1);
    int idx = 0;

    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi0z_[j];
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi0z_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi1z_[j];
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi1z_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];

    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * dux0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * duy0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hxl0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hyl0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * dux0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * duy0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hxr0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hyr0 * Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];

    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * dux1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * duy1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hxl1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hyl1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * dux1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * duy1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hxr1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -hyr1 * Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];

    // rows 20..25 absolute z bounds gradients
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j]; // zb0
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j]; // zb1
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi0z_[j]; // zl0
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi0z_[j]; // zr0
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi1z_[j]; // zl1
    for (int j = 0; j < n_coeff_z_; ++j) jac_g[idx++] = Phi1z_[j]; // zr1
}

void BoundaryConstraint10Var::getBounds(double* g_lb, double* g_ub, int offset) {
    (void)offset;
    g_lb[0] = -base_xy_tol_; g_ub[0] = base_xy_tol_;
    g_lb[1] = -base_xy_tol_; g_ub[1] = base_xy_tol_;
    g_lb[2] = -base_z_tol_;  g_ub[2] = base_z_tol_;
    g_lb[3] = -base_yaw_tol_; g_ub[3] = base_yaw_tol_;
    g_lb[4] = -base_xy_tol_; g_ub[4] = base_xy_tol_;
    g_lb[5] = -base_xy_tol_; g_ub[5] = base_xy_tol_;
    g_lb[6] = -base_z_tol_;  g_ub[6] = base_z_tol_;
    g_lb[7] = -base_yaw_tol_; g_ub[7] = base_yaw_tol_;
    for (int i = 8; i < 20; ++i) {
        const bool z_row = (i == 10 || i == 13 || i == 16 || i == 19);
        if (z_row) {
            g_lb[i] = -wheel_z_tol_;
            g_ub[i] = wheel_z_tol_;
        } else {
            g_lb[i] = -wheel_xy_tol_;
            g_ub[i] = wheel_xy_tol_;
        }
    }
    g_lb[20] = base_z_min_;  g_ub[20] = base_z_max_;
    g_lb[21] = base_z_min_;  g_ub[21] = base_z_max_;
    g_lb[22] = wheel_z_min_; g_ub[22] = wheel_z_max_;
    g_lb[23] = wheel_z_min_; g_ub[23] = wheel_z_max_;
    g_lb[24] = wheel_z_min_; g_ub[24] = wheel_z_max_;
    g_lb[25] = wheel_z_min_; g_ub[25] = wheel_z_max_;
}

int BoundaryConstraint10Var::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol) return 42 * n_coeff_ + 8 * n_coeff_z_;
    int idx = 0;
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 0; jCol[idx++] = idx_cbx_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 1; jCol[idx++] = idx_cby_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 2; jCol[idx++] = idx_cbz_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 3; jCol[idx++] = idx_cpsi_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 4; jCol[idx++] = idx_cbx_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 5; jCol[idx++] = idx_cby_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 6; jCol[idx++] = idx_cbz_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 7; jCol[idx++] = idx_cpsi_ + j; }

    auto add_row = [&](int row, int c0, int w0, int c1, int w1, int c2, int w2) {
        for (int j = 0; j < w0; ++j) { iRow[idx] = offset + row; jCol[idx++] = c0 + j; }
        for (int j = 0; j < w1; ++j) { iRow[idx] = offset + row; jCol[idx++] = c1 + j; }
        for (int j = 0; j < w2; ++j) { iRow[idx] = offset + row; jCol[idx++] = c2 + j; }
    };
    add_row(8, idx_clx_, n_coeff_, idx_cbx_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(9, idx_cly_, n_coeff_, idx_cby_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(10, idx_clx_, n_coeff_, idx_cly_, n_coeff_, idx_clz_, n_coeff_z_);
    add_row(11, idx_crx_, n_coeff_, idx_cbx_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(12, idx_cry_, n_coeff_, idx_cby_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(13, idx_crx_, n_coeff_, idx_cry_, n_coeff_, idx_crz_, n_coeff_z_);
    add_row(14, idx_clx_, n_coeff_, idx_cbx_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(15, idx_cly_, n_coeff_, idx_cby_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(16, idx_clx_, n_coeff_, idx_cly_, n_coeff_, idx_clz_, n_coeff_z_);
    add_row(17, idx_crx_, n_coeff_, idx_cbx_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(18, idx_cry_, n_coeff_, idx_cby_, n_coeff_, idx_cpsi_, n_coeff_);
    add_row(19, idx_crx_, n_coeff_, idx_cry_, n_coeff_, idx_crz_, n_coeff_z_);
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 20; jCol[idx++] = idx_cbz_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 21; jCol[idx++] = idx_cbz_ + j; }
    for (int j = 0; j < n_coeff_z_; ++j) { iRow[idx] = offset + 22; jCol[idx++] = idx_clz_ + j; }
    for (int j = 0; j < n_coeff_z_; ++j) { iRow[idx] = offset + 23; jCol[idx++] = idx_crz_ + j; }
    for (int j = 0; j < n_coeff_z_; ++j) { iRow[idx] = offset + 24; jCol[idx++] = idx_clz_ + j; }
    for (int j = 0; j < n_coeff_z_; ++j) { iRow[idx] = offset + 25; jCol[idx++] = idx_crz_ + j; }
    return idx;
}
