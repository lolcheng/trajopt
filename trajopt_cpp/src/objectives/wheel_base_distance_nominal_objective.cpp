#include "objectives/wheel_base_distance_nominal_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

WheelBaseDistanceNominalObjective::WheelBaseDistanceNominalObjective(
    int n_coeff, int n_sample,
    double nominal_dist, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      nominal_dist_sq_(nominal_dist * nominal_dist),
      weight_(weight < 0 ? 1.0 / n_sample : weight),
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
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void WheelBaseDistanceNominalObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double WheelBaseDistanceNominalObjective::computeDistSq(
    double xb, double yb, double zb,
    double xw, double yw, double zw) const {
    double dx = xb - xw;
    double dy = yb - yw;
    double dz = zb - zw;
    return dx * dx + dy * dy + dz * dz;
}

double WheelBaseDistanceNominalObjective::eval_f(const double* x) {
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
    
    // 计算目标值: sum((dist² - nominal²)²)
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        double dist2_l = computeDistSq(Xb[i], Yb[i], Zb[i], Xl[i], Yl[i], Zl[i]);
        double dist2_r = computeDistSq(Xb[i], Yb[i], Zb[i], Xr[i], Yr[i], Zr[i]);
        
        double diff_l = dist2_l - nominal_dist_sq_;
        double diff_r = dist2_r - nominal_dist_sq_;
        obj += diff_l * diff_l + diff_r * diff_r;
    }
    
    return weight_ * obj;
}

void WheelBaseDistanceNominalObjective::eval_grad_f(const double* x, double* grad_f) {
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
    
    // 计算梯度
    // d/dcoeff (dist² - nominal²)² = 2 * (dist² - nominal²) * d(dist²)/dcoeff
    for (int i = 0; i < n_sample_; ++i) {
        double dist2_l = computeDistSq(Xb[i], Yb[i], Zb[i], Xl[i], Yl[i], Zl[i]);
        double dist2_r = computeDistSq(Xb[i], Yb[i], Zb[i], Xr[i], Yr[i], Zr[i]);
        
        double diff_l = dist2_l - nominal_dist_sq_;
        double diff_r = dist2_r - nominal_dist_sq_;
        
        double factor_l = 2.0 * weight_ * diff_l;
        double factor_r = 2.0 * weight_ * diff_r;
        
        double dx_l = Xb[i] - Xl[i];
        double dy_l = Yb[i] - Yl[i];
        double dz_l = Zb[i] - Zl[i];
        double dx_r = Xb[i] - Xr[i];
        double dy_r = Yb[i] - Yr[i];
        double dz_r = Zb[i] - Zr[i];
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 左轮距离的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cbx[j] += factor_l * 2.0 * dx_l * Phi_row[j];
            grad_cby[j] += factor_l * 2.0 * dy_l * Phi_row[j];
            grad_cbz[j] += factor_l * 2.0 * dz_l * Phi_row[j];
            grad_clx[j] -= factor_l * 2.0 * dx_l * Phi_row[j];
            grad_cly[j] -= factor_l * 2.0 * dy_l * Phi_row[j];
            grad_clz[j] -= factor_l * 2.0 * dz_l * Phi_row[j];
        }
        
        // 右轮距离的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cbx[j] += factor_r * 2.0 * dx_r * Phi_row[j];
            grad_cby[j] += factor_r * 2.0 * dy_r * Phi_row[j];
            grad_cbz[j] += factor_r * 2.0 * dz_r * Phi_row[j];
            grad_crx[j] -= factor_r * 2.0 * dx_r * Phi_row[j];
            grad_cry[j] -= factor_r * 2.0 * dy_r * Phi_row[j];
            grad_crz[j] -= factor_r * 2.0 * dz_r * Phi_row[j];
        }
    }
}

