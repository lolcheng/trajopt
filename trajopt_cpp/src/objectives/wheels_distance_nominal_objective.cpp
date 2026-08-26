#include "objectives/wheels_distance_nominal_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

WheelsDistanceNominalObjective::WheelsDistanceNominalObjective(
    int n_coeff, int n_sample,
    double nominal_dist, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      nominal_dist_sq_(nominal_dist * nominal_dist),
      weight_(weight < 0 ? 1.0 / n_sample : weight),
      enabled_(true) {
    
    // 变量索引（根据TrajOptNLP中的布局）
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

void WheelsDistanceNominalObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double WheelsDistanceNominalObjective::computeDistSq(
    double xl, double yl, double zl,
    double xr, double yr, double zr) const {
    double dx = xr - xl;
    double dy = yr - yl;
    double dz = zr - zl;
    return dx * dx + dy * dy + dz * dz;
}

double WheelsDistanceNominalObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取系数
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轨迹
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(clz, Zl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    evalTrajectory(crz, Zr.data());
    
    // 计算目标值: sum((dist² - nominal²)²)
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        double dist2 = computeDistSq(Xl[i], Yl[i], Zl[i], Xr[i], Yr[i], Zr[i]);
        double diff = dist2 - nominal_dist_sq_;
        obj += diff * diff;
    }
    
    return weight_ * obj;
}

void WheelsDistanceNominalObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取系数
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轨迹
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(clz, Zl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    evalTrajectory(crz, Zr.data());
    
    double* grad_clx = grad_f + idx_clx_;
    double* grad_cly = grad_f + idx_cly_;
    double* grad_clz = grad_f + idx_clz_;
    double* grad_crx = grad_f + idx_crx_;
    double* grad_cry = grad_f + idx_cry_;
    double* grad_crz = grad_f + idx_crz_;
    
    // 计算梯度
    // d/dcoeff (dist² - nominal²)² = 2 * (dist² - nominal²) * d(dist²)/dcoeff
    for (int i = 0; i < n_sample_; ++i) {
        double dist2 = computeDistSq(Xl[i], Yl[i], Zl[i], Xr[i], Yr[i], Zr[i]);
        double diff = dist2 - nominal_dist_sq_;
        double factor = 2.0 * weight_ * diff;
        
        double dx = Xr[i] - Xl[i];
        double dy = Yr[i] - Yl[i];
        double dz = Zr[i] - Zl[i];
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对左轮系数的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_clx[j] -= factor * 2.0 * dx * Phi_row[j];
            grad_cly[j] -= factor * 2.0 * dy * Phi_row[j];
            grad_clz[j] -= factor * 2.0 * dz * Phi_row[j];
        }
        
        // 对右轮系数的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_crx[j] += factor * 2.0 * dx * Phi_row[j];
            grad_cry[j] += factor * 2.0 * dy * Phi_row[j];
            grad_crz[j] += factor * 2.0 * dz * Phi_row[j];
        }
    }
}

