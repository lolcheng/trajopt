#include "objectives/torque_balance_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

TorqueBalanceObjective::TorqueBalanceObjective(int n_coeff, int n_sample, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      weight_(weight < 0 ? 1.0 / n_sample : weight),  // 降低默认权重：10 -> 1
      enabled_(true) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = 7 * n_coeff_;
    idx_cry_ = 8 * n_coeff_;
    idx_crz_ = 9 * n_coeff_;
    
    int n_poly_vars = 10 * n_coeff_;
    idx_fLx_ = n_poly_vars;
    idx_fLy_ = idx_fLx_ + n_sample_;
    idx_fLz_ = idx_fLy_ + n_sample_;
    idx_fRx_ = idx_fLz_ + n_sample_;
    idx_fRy_ = idx_fRx_ + n_sample_;
    idx_fRz_ = idx_fRy_ + n_sample_;
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void TorqueBalanceObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void TorqueBalanceObjective::crossProduct(
    double ax, double ay, double az,
    double bx, double by, double bz,
    double& cx, double& cy, double& cz) const {
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}

double TorqueBalanceObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 提取接触力
    const double* fLx = x + idx_fLx_;
    const double* fLy = x + idx_fLy_;
    const double* fLz = x + idx_fLz_;
    const double* fRx = x + idx_fRx_;
    const double* fRy = x + idx_fRy_;
    const double* fRz = x + idx_fRz_;
    
    // 评估轨迹
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Zb(n_sample_);
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cbz, Zb.data());
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(clz, Zl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    evalTrajectory(crz, Zr.data());
    
    // 计算目标值: sum(||rL × fL + rR × fR||²)
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        // 计算相对位置向量
        double rLx = Xl[i] - Xb[i];
        double rLy = Yl[i] - Yb[i];
        double rLz = Zl[i] - Zb[i];
        double rRx = Xr[i] - Xb[i];
        double rRy = Yr[i] - Yb[i];
        double rRz = Zr[i] - Zb[i];
        
        // 计算力矩
        double tauLx, tauLy, tauLz;
        double tauRx, tauRy, tauRz;
        crossProduct(rLx, rLy, rLz, fLx[i], fLy[i], fLz[i], tauLx, tauLy, tauLz);
        crossProduct(rRx, rRy, rRz, fRx[i], fRy[i], fRz[i], tauRx, tauRy, tauRz);
        
        // 合矩
        double tau_x = tauLx + tauRx;
        double tau_y = tauLy + tauRy;
        double tau_z = tauLz + tauRz;
        
        obj += tau_x * tau_x + tau_y * tau_y + tau_z * tau_z;
    }
    
    return weight_ * obj;
}

void TorqueBalanceObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 提取接触力
    const double* fLx = x + idx_fLx_;
    const double* fLy = x + idx_fLy_;
    const double* fLz = x + idx_fLz_;
    const double* fRx = x + idx_fRx_;
    const double* fRy = x + idx_fRy_;
    const double* fRz = x + idx_fRz_;
    
    // 评估轨迹
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Zb(n_sample_);
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cbz, Zb.data());
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(clz, Zl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    evalTrajectory(crz, Zr.data());
    
    double* grad_cbx = grad_f + idx_cbx_;
    double* grad_cby = grad_f + idx_cby_;
    double* grad_cbz = grad_f + idx_cbz_;
    double* grad_clx = grad_f + idx_clx_;
    double* grad_cly = grad_f + idx_cly_;
    double* grad_clz = grad_f + idx_clz_;
    double* grad_crx = grad_f + idx_crx_;
    double* grad_cry = grad_f + idx_cry_;
    double* grad_crz = grad_f + idx_crz_;
    double* grad_fLx = grad_f + idx_fLx_;
    double* grad_fLy = grad_f + idx_fLy_;
    double* grad_fLz = grad_f + idx_fLz_;
    double* grad_fRx = grad_f + idx_fRx_;
    double* grad_fRy = grad_f + idx_fRy_;
    double* grad_fRz = grad_f + idx_fRz_;
    
    // 计算梯度
    // d/dcoeff ||rL × fL + rR × fR||² = 2 * (rL × fL + rR × fR) · (drL/dcoeff × fL + rL × dfL/dcoeff)
    for (int i = 0; i < n_sample_; ++i) {
        // 计算相对位置向量
        double rLx = Xl[i] - Xb[i];
        double rLy = Yl[i] - Yb[i];
        double rLz = Zl[i] - Zb[i];
        double rRx = Xr[i] - Xb[i];
        double rRy = Yr[i] - Yb[i];
        double rRz = Zr[i] - Zb[i];
        
        // 计算力矩
        double tauLx, tauLy, tauLz;
        double tauRx, tauRy, tauRz;
        crossProduct(rLx, rLy, rLz, fLx[i], fLy[i], fLz[i], tauLx, tauLy, tauLz);
        crossProduct(rRx, rRy, rRz, fRx[i], fRy[i], fRz[i], tauRx, tauRy, tauRz);
        
        // 合矩
        double tau_x = tauLx + tauRx;
        double tau_y = tauLy + tauRy;
        double tau_z = tauLz + tauRz;
        
        double factor = 2.0 * weight_;
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对接触力的梯度
        // d(tau)/dfL = rL × (单位向量)
        double dtau_x_dfLx = 0.0, dtau_x_dfLy = -rLz, dtau_x_dfLz = rLy;
        double dtau_y_dfLx = rLz, dtau_y_dfLy = 0.0, dtau_y_dfLz = -rLx;
        double dtau_z_dfLx = -rLy, dtau_z_dfLy = rLx, dtau_z_dfLz = 0.0;
        
        grad_fLx[i] += factor * (tau_x * dtau_x_dfLx + tau_y * dtau_y_dfLx + tau_z * dtau_z_dfLx);
        grad_fLy[i] += factor * (tau_x * dtau_x_dfLy + tau_y * dtau_y_dfLy + tau_z * dtau_z_dfLy);
        grad_fLz[i] += factor * (tau_x * dtau_x_dfLz + tau_y * dtau_y_dfLz + tau_z * dtau_z_dfLz);
        
        double dtau_x_dfRx = 0.0, dtau_x_dfRy = -rRz, dtau_x_dfRz = rRy;
        double dtau_y_dfRx = rRz, dtau_y_dfRy = 0.0, dtau_y_dfRz = -rRx;
        double dtau_z_dfRx = -rRy, dtau_z_dfRy = rRx, dtau_z_dfRz = 0.0;
        
        grad_fRx[i] += factor * (tau_x * dtau_x_dfRx + tau_y * dtau_y_dfRx + tau_z * dtau_z_dfRx);
        grad_fRy[i] += factor * (tau_x * dtau_x_dfRy + tau_y * dtau_y_dfRy + tau_z * dtau_z_dfRy);
        grad_fRz[i] += factor * (tau_x * dtau_x_dfRz + tau_y * dtau_y_dfRz + tau_z * dtau_z_dfRz);
        
        // 对位置的梯度（通过rL和rR）
        // d(tau)/drL = fL × (单位向量)
        double dtau_x_drLx = 0.0, dtau_x_drLy = fLz[i], dtau_x_drLz = -fLy[i];
        double dtau_y_drLx = -fLz[i], dtau_y_drLy = 0.0, dtau_y_drLz = fLx[i];
        double dtau_z_drLx = fLy[i], dtau_z_drLy = -fLx[i], dtau_z_drLz = 0.0;
        
        // drL/dXl = 1, drL/dXb = -1
        for (int j = 0; j < n_coeff_; ++j) {
            grad_clx[j] += factor * (tau_x * dtau_x_drLx + tau_y * dtau_y_drLx + tau_z * dtau_z_drLx) * Phi_row[j];
            grad_cly[j] += factor * (tau_x * dtau_x_drLy + tau_y * dtau_y_drLy + tau_z * dtau_z_drLy) * Phi_row[j];
            grad_clz[j] += factor * (tau_x * dtau_x_drLz + tau_y * dtau_y_drLz + tau_z * dtau_z_drLz) * Phi_row[j];
            grad_cbx[j] -= factor * (tau_x * dtau_x_drLx + tau_y * dtau_y_drLx + tau_z * dtau_z_drLx) * Phi_row[j];
            grad_cby[j] -= factor * (tau_x * dtau_x_drLy + tau_y * dtau_y_drLy + tau_z * dtau_z_drLy) * Phi_row[j];
            grad_cbz[j] -= factor * (tau_x * dtau_x_drLz + tau_y * dtau_y_drLz + tau_z * dtau_z_drLz) * Phi_row[j];
        }
        
        double dtau_x_drRx = 0.0, dtau_x_drRy = fRz[i], dtau_x_drRz = -fRy[i];
        double dtau_y_drRx = -fRz[i], dtau_y_drRy = 0.0, dtau_y_drRz = fRx[i];
        double dtau_z_drRx = fRy[i], dtau_z_drRy = -fRx[i], dtau_z_drRz = 0.0;
        
        for (int j = 0; j < n_coeff_; ++j) {
            grad_crx[j] += factor * (tau_x * dtau_x_drRx + tau_y * dtau_y_drRx + tau_z * dtau_z_drRx) * Phi_row[j];
            grad_cry[j] += factor * (tau_x * dtau_x_drRy + tau_y * dtau_y_drRy + tau_z * dtau_z_drRy) * Phi_row[j];
            grad_crz[j] += factor * (tau_x * dtau_x_drRz + tau_y * dtau_y_drRz + tau_z * dtau_z_drRz) * Phi_row[j];
            grad_cbx[j] -= factor * (tau_x * dtau_x_drRx + tau_y * dtau_y_drRx + tau_z * dtau_z_drRx) * Phi_row[j];
            grad_cby[j] -= factor * (tau_x * dtau_x_drRy + tau_y * dtau_y_drRy + tau_z * dtau_z_drRy) * Phi_row[j];
            grad_cbz[j] -= factor * (tau_x * dtau_x_drRz + tau_y * dtau_y_drRz + tau_z * dtau_z_drRz) * Phi_row[j];
        }
    }
}

