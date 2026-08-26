#include "constraints/boundary_constraint.hpp"
#include "polynomial.hpp"
#include "polynomial_basis.hpp"
#include "chebyshev_basis.hpp"
#include <cstring>
#include <cmath>
#include <memory>

BoundaryConstraint::BoundaryConstraint(int n_coeff,
                                       double x0, double y0, double psi0,
                                       double x1, double y1, double psi1,
                                       double clear_z,
                                       const RBFTerrain& terrain,
                                       std::shared_ptr<BasisFunction> basis,
                                       double wheels_nominal)
    : n_coeff_(n_coeff),
      x0_(x0), y0_(y0), psi0_(psi0),
      x1_(x1), y1_(y1), psi1_(psi1),
      clear_z_(clear_z),
      wheels_nominal_(wheels_nominal),
      terrain_(terrain),
      Phi0_(new double[n_coeff_]),
      Phi1_(new double[n_coeff_]) {
    
    // 如果未提供基函数，默认使用切比雪夫基
    if (basis == nullptr) {
        basis_ = std::make_shared<ChebyshevBasis>();
    } else {
        basis_ = basis;
    }
    
    // 预计算基矩阵（s=0和s=1）
    double s0 = 0.0;
    double s1 = 1.0;
    basis_->computeBasisMatrix(&s0, 1, n_coeff_, Phi0_);
    basis_->computeBasisMatrix(&s1, 1, n_coeff_, Phi1_);
    
    // 变量索引（根据TrajOptNLP中的布局）
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = 7 * n_coeff_;
    idx_cry_ = 8 * n_coeff_;
    idx_crz_ = 9 * n_coeff_;
}

BoundaryConstraint::~BoundaryConstraint() {
    delete[] Phi0_;
    delete[] Phi1_;
}

void BoundaryConstraint::eval_g(const double* x, double* g) {
    // 提取系数
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
    
    // 评估起点 (s=0)
    double xb0 = basis_->evalPoint(0.0, cbx, n_coeff_);
    double yb0 = basis_->evalPoint(0.0, cby, n_coeff_);
    double zb0 = basis_->evalPoint(0.0, cbz, n_coeff_);
    double psi0 = basis_->evalPoint(0.0, cpsi, n_coeff_);
    
    // 评估终点 (s=1)
    double xb1 = basis_->evalPoint(1.0, cbx, n_coeff_);
    double yb1 = basis_->evalPoint(1.0, cby, n_coeff_);
    double zb1 = basis_->evalPoint(1.0, cbz, n_coeff_);
    double psi1 = basis_->evalPoint(1.0, cpsi, n_coeff_);
    
    // 评估轮子起点
    double xl0 = basis_->evalPoint(0.0, clx, n_coeff_);
    double yl0 = basis_->evalPoint(0.0, cly, n_coeff_);
    double zl0 = basis_->evalPoint(0.0, clz, n_coeff_);
    double xr0 = basis_->evalPoint(0.0, crx, n_coeff_);
    double yr0 = basis_->evalPoint(0.0, cry, n_coeff_);
    double zr0 = basis_->evalPoint(0.0, crz, n_coeff_);
    
    // 评估轮子终点
    double xl1 = basis_->evalPoint(1.0, clx, n_coeff_);
    double yl1 = basis_->evalPoint(1.0, cly, n_coeff_);
    double zl1 = basis_->evalPoint(1.0, clz, n_coeff_);
    double xr1 = basis_->evalPoint(1.0, crx, n_coeff_);
    double yr1 = basis_->evalPoint(1.0, cry, n_coeff_);
    double zr1 = basis_->evalPoint(1.0, crz, n_coeff_);
    
    // 计算起点终点高度（使用RBF地形）
    double h0 = terrain_.height(x0_, y0_);
    double h1 = terrain_.height(x1_, y1_);
    double z0_target = h0 + clear_z_;
    double z1_target = h1 + clear_z_;
    
    // 计算轮子位置（基于base和yaw）
    double ux0 = -std::sin(psi0);
    double uy0 = std::cos(psi0);
    double ux1 = -std::sin(psi1);
    double uy1 = std::cos(psi1);
    
    double xl0_target = xb0 + (wheels_nominal_ / 2.0) * ux0;
    double yl0_target = yb0 + (wheels_nominal_ / 2.0) * uy0;
    double xr0_target = xb0 - (wheels_nominal_ / 2.0) * ux0;
    double yr0_target = yb0 - (wheels_nominal_ / 2.0) * uy0;
    
    double xl1_target = xb1 + (wheels_nominal_ / 2.0) * ux1;
    double yl1_target = yb1 + (wheels_nominal_ / 2.0) * uy1;
    double xr1_target = xb1 - (wheels_nominal_ / 2.0) * ux1;
    double yr1_target = yb1 - (wheels_nominal_ / 2.0) * uy1;
    
    // 轮子z坐标应该等于地形高度（触地）
    double zl0_target = terrain_.height(xl0_target, yl0_target);
    double zr0_target = terrain_.height(xr0_target, yr0_target);
    double zl1_target = terrain_.height(xl1_target, yl1_target);
    double zr1_target = terrain_.height(xr1_target, yr1_target);
    
    // Base约束（前8个）
    g[0] = xb0 - x0_;           // x(0) = x0
    g[1] = yb0 - y0_;           // y(0) = y0
    g[2] = zb0 - z0_target;     // z(0) = h(x0,y0) + clear_z
    g[3] = psi0 - psi0_;        // psi(0) = psi0
    g[4] = xb1 - x1_;           // x(1) = x1
    g[5] = yb1 - y1_;           // y(1) = y1
    g[6] = zb1 - z1_target;     // z(1) = h(x1,y1) + clear_z
    g[7] = psi1 - psi1_;        // psi(1) = psi1
    
    // 轮子起点约束
    g[8] = xl0 - xl0_target;   // 左轮x(0) = xb0 + (w/2)*ux0
    g[9] = yl0 - yl0_target;   // 左轮y(0) = yb0 + (w/2)*uy0
    g[10] = zl0 - zl0_target;  // 左轮z(0) = h(xl0, yl0) (触地)
    g[11] = xr0 - xr0_target;  // 右轮x(0) = xb0 - (w/2)*ux0
    g[12] = yr0 - yr0_target;  // 右轮y(0) = yb0 - (w/2)*uy0
    g[13] = zr0 - zr0_target;  // 右轮z(0) = h(xr0, yr0) (触地)
    
    // 轮子终点约束
    g[14] = xl1 - xl1_target;  // 左轮x(1) = xb1 + (w/2)*ux1
    g[15] = yl1 - yl1_target;  // 左轮y(1) = yb1 + (w/2)*uy1
    g[16] = zl1 - zl1_target;  // 左轮z(1) = h(xl1, yl1) (触地)
    g[17] = xr1 - xr1_target;  // 右轮x(1) = xb1 - (w/2)*ux1
    g[18] = yr1 - yr1_target;  // 右轮y(1) = yb1 - (w/2)*uy1
    g[19] = zr1 - zr1_target;  // 右轮z(1) = h(xr1, yr1) (触地)
}

void BoundaryConstraint::eval_jac_g(const double* x, double* jac_g, int offset) {
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    const double* cpsi = x + idx_cpsi_;
    
    // 评估起点和终点的psi值
    double psi0 = basis_->evalPoint(0.0, cpsi, n_coeff_);
    double psi1 = basis_->evalPoint(1.0, cpsi, n_coeff_);
    
    // 计算ux, uy及其对psi的导数
    double ux0 = -std::sin(psi0);
    double uy0 = std::cos(psi0);
    double dux0_dpsi = -std::cos(psi0);
    double duy0_dpsi = -std::sin(psi0);
    
    double ux1 = -std::sin(psi1);
    double uy1 = std::cos(psi1);
    double dux1_dpsi = -std::cos(psi1);
    double duy1_dpsi = -std::sin(psi1);
    
    // 评估起点和终点的base位置
    double xb0 = basis_->evalPoint(0.0, cbx, n_coeff_);
    double yb0 = basis_->evalPoint(0.0, cby, n_coeff_);
    double xb1 = basis_->evalPoint(1.0, cbx, n_coeff_);
    double yb1 = basis_->evalPoint(1.0, cby, n_coeff_);
    
    // 计算轮子目标位置
    double xl0_target = xb0 + (wheels_nominal_ / 2.0) * ux0;
    double yl0_target = yb0 + (wheels_nominal_ / 2.0) * uy0;
    double xr0_target = xb0 - (wheels_nominal_ / 2.0) * ux0;
    double yr0_target = yb0 - (wheels_nominal_ / 2.0) * uy0;
    
    double xl1_target = xb1 + (wheels_nominal_ / 2.0) * ux1;
    double yl1_target = yb1 + (wheels_nominal_ / 2.0) * uy1;
    double xr1_target = xb1 - (wheels_nominal_ / 2.0) * ux1;
    double yr1_target = yb1 - (wheels_nominal_ / 2.0) * uy1;
    
    // 计算地形梯度（用于轮子z约束）
    double hxl0_x, hxl0_y, hxr0_x, hxr0_y;
    double hxl1_x, hxl1_y, hxr1_x, hxr1_y;
    terrain_.gradient(xl0_target, yl0_target, hxl0_x, hxl0_y);
    terrain_.gradient(xr0_target, yr0_target, hxr0_x, hxr0_y);
    terrain_.gradient(xl1_target, yl1_target, hxl1_x, hxl1_y);
    terrain_.gradient(xr1_target, yr1_target, hxr1_x, hxr1_y);
    
    // 按照getJacobianStructure的顺序填充雅可比矩阵
    int idx = 0;
    (void)offset;
    
    // Base约束（前8个约束）
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < n_coeff_; ++j) {
            if (i == 0) jac_g[idx] = Phi0_[j];
            else if (i == 1) jac_g[idx] = Phi0_[j];
            else if (i == 2) jac_g[idx] = Phi0_[j];
            else jac_g[idx] = Phi0_[j];
            idx++;
        }
    }
    for (int i = 4; i < 8; ++i) {
        for (int j = 0; j < n_coeff_; ++j) {
            jac_g[idx] = Phi1_[j];
            idx++;
        }
    }
    
    // 轮子起点约束（g[8-13]）
    // g[8]: xl0约束
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi0_[j];  // ∂g[8]/∂clx
        jac_g[idx++] = -Phi0_[j];  // ∂g[8]/∂cbx
        jac_g[idx++] = -(wheels_nominal_ / 2.0) * dux0_dpsi * Phi0_[j];  // ∂g[8]/∂cpsi
    }
    // g[9]: yl0约束
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi0_[j];  // ∂g[9]/∂cly
        jac_g[idx++] = -Phi0_[j];  // ∂g[9]/∂cby
        jac_g[idx++] = -(wheels_nominal_ / 2.0) * duy0_dpsi * Phi0_[j];  // ∂g[9]/∂cpsi
    }
    // g[10]: zl0约束
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi0_[j];  // ∂g[10]/∂clz
        jac_g[idx++] = -hxl0_x * Phi0_[j];  // ∂g[10]/∂cbx
        jac_g[idx++] = -hxl0_y * Phi0_[j];  // ∂g[10]/∂cby
        jac_g[idx++] = -(wheels_nominal_ / 2.0) * (hxl0_x * dux0_dpsi + hxl0_y * duy0_dpsi) * Phi0_[j];  // ∂g[10]/∂cpsi
    }
    // g[11]: xr0约束
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi0_[j];  // ∂g[11]/∂crx
        jac_g[idx++] = -Phi0_[j];  // ∂g[11]/∂cbx
        jac_g[idx++] = (wheels_nominal_ / 2.0) * dux0_dpsi * Phi0_[j];  // ∂g[11]/∂cpsi
    }
    // g[12]: yr0约束
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi0_[j];  // ∂g[12]/∂cry
        jac_g[idx++] = -Phi0_[j];  // ∂g[12]/∂cby
        jac_g[idx++] = (wheels_nominal_ / 2.0) * duy0_dpsi * Phi0_[j];  // ∂g[12]/∂cpsi
    }
    // g[13]: zr0约束
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi0_[j];  // ∂g[13]/∂crz
        jac_g[idx++] = -hxr0_x * Phi0_[j];  // ∂g[13]/∂cbx
        jac_g[idx++] = -hxr0_y * Phi0_[j];  // ∂g[13]/∂cby
        jac_g[idx++] = (wheels_nominal_ / 2.0) * (hxr0_x * dux0_dpsi + hxr0_y * duy0_dpsi) * Phi0_[j];  // ∂g[13]/∂cpsi
    }
    
    // 轮子终点约束（g[14-19]），类似起点但使用Phi1
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi1_[j];  // g[14]/∂clx
        jac_g[idx++] = -Phi1_[j];  // g[14]/∂cbx
        jac_g[idx++] = -(wheels_nominal_ / 2.0) * dux1_dpsi * Phi1_[j];  // g[14]/∂cpsi
    }
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi1_[j];  // g[15]/∂cly
        jac_g[idx++] = -Phi1_[j];  // g[15]/∂cby
        jac_g[idx++] = -(wheels_nominal_ / 2.0) * duy1_dpsi * Phi1_[j];  // g[15]/∂cpsi
    }
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi1_[j];  // g[16]/∂clz
        jac_g[idx++] = -hxl1_x * Phi1_[j];  // g[16]/∂cbx
        jac_g[idx++] = -hxl1_y * Phi1_[j];  // g[16]/∂cby
        jac_g[idx++] = -(wheels_nominal_ / 2.0) * (hxl1_x * dux1_dpsi + hxl1_y * duy1_dpsi) * Phi1_[j];  // g[16]/∂cpsi
    }
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi1_[j];  // g[17]/∂crx
        jac_g[idx++] = -Phi1_[j];  // g[17]/∂cbx
        jac_g[idx++] = (wheels_nominal_ / 2.0) * dux1_dpsi * Phi1_[j];  // g[17]/∂cpsi
    }
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi1_[j];  // g[18]/∂cry
        jac_g[idx++] = -Phi1_[j];  // g[18]/∂cby
        jac_g[idx++] = (wheels_nominal_ / 2.0) * duy1_dpsi * Phi1_[j];  // g[18]/∂cpsi
    }
    for (int j = 0; j < n_coeff_; ++j) {
        jac_g[idx++] = Phi1_[j];  // g[19]/∂crz
        jac_g[idx++] = -hxr1_x * Phi1_[j];  // g[19]/∂cbx
        jac_g[idx++] = -hxr1_y * Phi1_[j];  // g[19]/∂cby
        jac_g[idx++] = (wheels_nominal_ / 2.0) * (hxr1_x * dux1_dpsi + hxr1_y * duy1_dpsi) * Phi1_[j];  // g[19]/∂cpsi
    }
}

void BoundaryConstraint::getBounds(double* g_lb, double* g_ub, int offset) {
    // 等式约束：g = 0（调用方传入的 g_lb 已偏移，此处只写 [0..19]）
    (void)offset;
    for (int i = 0; i < 20; ++i) {
        g_lb[i] = 0.0;
        g_ub[i] = 0.0;
    }
}

int BoundaryConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        // 只返回非零元素数量
        // Base约束：8个约束 × 1个变量组 × n_coeff = 8 * n_coeff
        // 轮子起点约束：g[8-13]
        //   g[8]: 3个变量组 (clx, cbx, cpsi) × n_coeff
        //   g[9]: 3个变量组 (cly, cby, cpsi) × n_coeff
        //   g[10]: 4个变量组 (clz, cbx, cby, cpsi) × n_coeff
        //   g[11]: 3个变量组 (crx, cbx, cpsi) × n_coeff
        //   g[12]: 3个变量组 (cry, cby, cpsi) × n_coeff
        //   g[13]: 4个变量组 (crz, cbx, cby, cpsi) × n_coeff
        //   起点小计：(3+3+4+3+3+4) × n_coeff = 20 * n_coeff
        // 轮子终点约束：g[14-19]，同样 20 * n_coeff
        // 总共：8 * n_coeff + 20 * n_coeff + 20 * n_coeff = 48 * n_coeff
        return 48 * n_coeff_;
    }
    
    int idx = 0;
    
    // Base约束（前8个约束，每个约束只依赖一个变量组）
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = i + offset;
            if (i == 0) jCol[idx] = idx_cbx_ + j;
            else if (i == 1) jCol[idx] = idx_cby_ + j;
            else if (i == 2) jCol[idx] = idx_cbz_ + j;
            else jCol[idx] = idx_cpsi_ + j;
            idx++;
        }
    }
    
    for (int i = 4; i < 8; ++i) {
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = i + offset;
            if (i == 4) jCol[idx] = idx_cbx_ + j;
            else if (i == 5) jCol[idx] = idx_cby_ + j;
            else if (i == 6) jCol[idx] = idx_cbz_ + j;
            else jCol[idx] = idx_cpsi_ + j;
            idx++;
        }
    }
    
    // 轮子起点约束（g[8-13]）
    // g[8]: xl0约束，依赖clx, cbx, cpsi
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 8 + offset; jCol[idx] = idx_clx_ + j; idx++;
        iRow[idx] = 8 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 8 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    // g[9]: yl0约束，依赖cly, cby, cpsi
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 9 + offset; jCol[idx] = idx_cly_ + j; idx++;
        iRow[idx] = 9 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 9 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    // g[10]: zl0约束，依赖clz, cbx, cby, cpsi
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 10 + offset; jCol[idx] = idx_clz_ + j; idx++;
        iRow[idx] = 10 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 10 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 10 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    // g[11-13]: 右轮起点约束，类似
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 11 + offset; jCol[idx] = idx_crx_ + j; idx++;
        iRow[idx] = 11 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 11 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 12 + offset; jCol[idx] = idx_cry_ + j; idx++;
        iRow[idx] = 12 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 12 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 13 + offset; jCol[idx] = idx_crz_ + j; idx++;
        iRow[idx] = 13 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 13 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 13 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    
    // 轮子终点约束（g[14-19]），类似起点但使用Phi1
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 14 + offset; jCol[idx] = idx_clx_ + j; idx++;
        iRow[idx] = 14 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 14 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 15 + offset; jCol[idx] = idx_cly_ + j; idx++;
        iRow[idx] = 15 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 15 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 16 + offset; jCol[idx] = idx_clz_ + j; idx++;
        iRow[idx] = 16 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 16 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 16 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 17 + offset; jCol[idx] = idx_crx_ + j; idx++;
        iRow[idx] = 17 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 17 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 18 + offset; jCol[idx] = idx_cry_ + j; idx++;
        iRow[idx] = 18 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 18 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    for (int j = 0; j < n_coeff_; ++j) {
        iRow[idx] = 19 + offset; jCol[idx] = idx_crz_ + j; idx++;
        iRow[idx] = 19 + offset; jCol[idx] = idx_cbx_ + j; idx++;
        iRow[idx] = 19 + offset; jCol[idx] = idx_cby_ + j; idx++;
        iRow[idx] = 19 + offset; jCol[idx] = idx_cpsi_ + j; idx++;
    }
    
    return idx;
}

