#include "constraints/boundary_constraint_8var.hpp"
#include "chebyshev_basis.hpp"
#include <cmath>

BoundaryConstraint8Var::BoundaryConstraint8Var(int n_coeff,
                                               double x0, double y0, double psi0,
                                               double x1, double y1, double psi1,
                                               double clear_z,
                                               const RBFTerrain& terrain,
                                               std::shared_ptr<BasisFunction> basis,
                                               double wheels_nominal,
                                               double base_xy_tol,
                                               double base_z_tol,
                                               double base_yaw_tol,
                                               double wheel_xy_tol)
    : n_coeff_(n_coeff),
      x0_(x0), y0_(y0), psi0_(psi0),
      x1_(x1), y1_(y1), psi1_(psi1),
      clear_z_(clear_z),
      wheels_nominal_(wheels_nominal),
      base_xy_tol_(std::abs(base_xy_tol)),
      base_z_tol_(std::abs(base_z_tol)),
      base_yaw_tol_(std::abs(base_yaw_tol)),
      wheel_xy_tol_(std::abs(wheel_xy_tol)),
      terrain_(terrain),
      Phi0_(new double[n_coeff_]),
      Phi1_(new double[n_coeff_]) {
    basis_ = basis ? basis : std::make_shared<ChebyshevBasis>();
    double s0 = 0.0, s1 = 1.0;
    basis_->computeBasisMatrix(&s0, 1, n_coeff_, Phi0_);
    basis_->computeBasisMatrix(&s1, 1, n_coeff_, Phi1_);
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_crx_ = 6 * n_coeff_;
    idx_cry_ = 7 * n_coeff_;
}

BoundaryConstraint8Var::~BoundaryConstraint8Var() {
    delete[] Phi0_;
    delete[] Phi1_;
}

void BoundaryConstraint8Var::eval_g(const double* x, double* g) {
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    const double* cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;

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
    double xr0 = basis_->evalPoint(0.0, crx, n_coeff_);
    double yr0 = basis_->evalPoint(0.0, cry, n_coeff_);
    double xl1 = basis_->evalPoint(1.0, clx, n_coeff_);
    double yl1 = basis_->evalPoint(1.0, cly, n_coeff_);
    double xr1 = basis_->evalPoint(1.0, crx, n_coeff_);
    double yr1 = basis_->evalPoint(1.0, cry, n_coeff_);
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
    g[10] = xr0 - xr0_target;
    g[11] = yr0 - yr0_target;
    g[12] = xl1 - xl1_target;
    g[13] = yl1 - yl1_target;
    g[14] = xr1 - xr1_target;
    g[15] = yr1 - yr1_target;
}

void BoundaryConstraint8Var::eval_jac_g(const double* x, double* jac_g, int offset) {
    (void)offset;
    const double* cpsi = x + idx_cpsi_;
    double yaw0 = basis_->evalPoint(0.0, cpsi, n_coeff_);
    double yaw1 = basis_->evalPoint(1.0, cpsi, n_coeff_);
    double dux0 = -std::cos(yaw0), duy0 = -std::sin(yaw0);
    double dux1 = -std::cos(yaw1), duy1 = -std::sin(yaw1);
    int idx = 0;
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    // Rows 8..15 must follow the same nonzero ordering as getJacobianStructure:
    // [all wxy coeffs][all bxy coeffs][all cpsi coeffs] for each row.
    // row 8: g8 (xl0)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * dux0 * Phi0_[j];
    // row 9: g9 (yl0)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * duy0 * Phi0_[j];
    // row10: g10 (xr0)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * dux0 * Phi0_[j];
    // row11: g11 (yr0)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi0_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * duy0 * Phi0_[j];
    // row12: g12 (xl1)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * dux1 * Phi1_[j];
    // row13: g13 (yl1)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -0.5 * wheels_nominal_ * duy1 * Phi1_[j];
    // row14: g14 (xr1)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * dux1 * Phi1_[j];
    // row15: g15 (yr1)
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = -Phi1_[j];
    for (int j = 0; j < n_coeff_; ++j) jac_g[idx++] = 0.5 * wheels_nominal_ * duy1 * Phi1_[j];
}

void BoundaryConstraint8Var::getBounds(double* g_lb, double* g_ub, int offset) {
    (void)offset;
    // Base start/end x,y
    g_lb[0] = -base_xy_tol_; g_ub[0] = base_xy_tol_;
    g_lb[1] = -base_xy_tol_; g_ub[1] = base_xy_tol_;
    // Base start z
    g_lb[2] = -base_z_tol_; g_ub[2] = base_z_tol_;
    // Base start yaw
    g_lb[3] = -base_yaw_tol_; g_ub[3] = base_yaw_tol_;
    // Base end x,y
    g_lb[4] = -base_xy_tol_; g_ub[4] = base_xy_tol_;
    g_lb[5] = -base_xy_tol_; g_ub[5] = base_xy_tol_;
    // Base end z
    g_lb[6] = -base_z_tol_; g_ub[6] = base_z_tol_;
    // Base end yaw
    g_lb[7] = -base_yaw_tol_; g_ub[7] = base_yaw_tol_;

    // Wheel start/end x,y relative geometry
    for (int i = 8; i < 16; ++i) {
        g_lb[i] = -wheel_xy_tol_;
        g_ub[i] = wheel_xy_tol_;
    }
}

int BoundaryConstraint8Var::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (!iRow || !jCol) return 32 * n_coeff_;
    int idx = 0;
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 0; jCol[idx++] = idx_cbx_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 1; jCol[idx++] = idx_cby_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 2; jCol[idx++] = idx_cbz_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 3; jCol[idx++] = idx_cpsi_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 4; jCol[idx++] = idx_cbx_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 5; jCol[idx++] = idx_cby_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 6; jCol[idx++] = idx_cbz_ + j; }
    for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + 7; jCol[idx++] = idx_cpsi_ + j; }
    int rows[] = {8, 9, 10, 11, 12, 13, 14, 15};
    int wxy[] = {idx_clx_, idx_cly_, idx_crx_, idx_cry_, idx_clx_, idx_cly_, idx_crx_, idx_cry_};
    int bxy[] = {idx_cbx_, idx_cby_, idx_cbx_, idx_cby_, idx_cbx_, idx_cby_, idx_cbx_, idx_cby_};
    for (int r = 0; r < 8; ++r) {
        for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + rows[r]; jCol[idx++] = wxy[r] + j; }
        for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + rows[r]; jCol[idx++] = bxy[r] + j; }
        for (int j = 0; j < n_coeff_; ++j) { iRow[idx] = offset + rows[r]; jCol[idx++] = idx_cpsi_ + j; }
    }
    return idx;
}
