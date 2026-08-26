#include "objectives/wheel_terrain_distance_objective.hpp"
#include "polynomial.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

WheelTerrainDistanceObjective::WheelTerrainDistanceObjective(
    int n_coeff, int n_sample,
    const RBFTerrain& terrain,
    double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      terrain_(terrain),
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

void WheelTerrainDistanceObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double WheelTerrainDistanceObjective::eval_f(const double* x) {
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
    
    // 计算地形高度
    std::vector<double> Hl(n_sample_), Hr(n_sample_);
    terrain_.height(Xl.data(), Yl.data(), n_sample_, Hl.data());
    terrain_.height(Xr.data(), Yr.data(), n_sample_, Hr.data());
    
    // 计算目标值: sum((Zl - Hl)² + (Zr - Hr)²)
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        double dist_l = Zl[i] - Hl[i];
        double dist_r = Zr[i] - Hr[i];
        obj += dist_l * dist_l + dist_r * dist_r;
    }
    
    return weight_ * obj;
}

void WheelTerrainDistanceObjective::evalTrajectoryGrad(
    const double* traj_diff, const double* Phi_row, double* grad_coeff) const {
    // grad_coeff += traj_diff * Phi_row
    // 注意：这个方法实际上没有在eval_grad_f中使用，可以删除
    // 但为了接口完整性保留
    for (int j = 0; j < n_coeff_; ++j) {
        grad_coeff[j] += traj_diff[0] * Phi_row[j];
    }
}

void WheelTerrainDistanceObjective::eval_grad_f(const double* x, double* grad_f) {
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
    
    // 计算地形高度和梯度
    std::vector<double> Hl(n_sample_), Hr(n_sample_);
    std::vector<double> Hlx(n_sample_), Hly(n_sample_);
    std::vector<double> Hrx(n_sample_), Hry(n_sample_);
    
    terrain_.height(Xl.data(), Yl.data(), n_sample_, Hl.data());
    terrain_.height(Xr.data(), Yr.data(), n_sample_, Hr.data());
    terrain_.gradient(Xl.data(), Yl.data(), n_sample_, Hlx.data(), Hly.data());
    terrain_.gradient(Xr.data(), Yr.data(), n_sample_, Hrx.data(), Hry.data());
    
    // 计算梯度
    // 对于左轮: d/dcoeff (Zl - Hl)² = 2 * (Zl - Hl) * (dZl/dcoeff - dHl/dcoeff)
    // dHl/dcoeff = dHl/dXl * dXl/dcoeff + dHl/dYl * dYl/dcoeff
    
    double* grad_clx = grad_f + idx_clx_;
    double* grad_cly = grad_f + idx_cly_;
    double* grad_clz = grad_f + idx_clz_;
    double* grad_crx = grad_f + idx_crx_;
    double* grad_cry = grad_f + idx_cry_;
    double* grad_crz = grad_f + idx_crz_;
    
    for (int i = 0; i < n_sample_; ++i) {
        double dist_l = Zl[i] - Hl[i];
        double dist_r = Zr[i] - Hr[i];
        
        double factor_l = 2.0 * weight_ * dist_l;
        double factor_r = 2.0 * weight_ * dist_r;
        
        // 左轮z坐标的梯度
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            grad_clz[j] += factor_l * Phi_row[j];
        }
        
        // 左轮x坐标的梯度（通过Hl对Xl的依赖）
        for (int j = 0; j < n_coeff_; ++j) {
            grad_clx[j] -= factor_l * Hlx[i] * Phi_row[j];
        }
        
        // 左轮y坐标的梯度（通过Hl对Yl的依赖）
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cly[j] -= factor_l * Hly[i] * Phi_row[j];
        }
        
        // 右轮z坐标的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_crz[j] += factor_r * Phi_row[j];
        }
        
        // 右轮x坐标的梯度（通过Hr对Xr的依赖）
        for (int j = 0; j < n_coeff_; ++j) {
            grad_crx[j] -= factor_r * Hrx[i] * Phi_row[j];
        }
        
        // 右轮y坐标的梯度（通过Hr对Yr的依赖）
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cry[j] -= factor_r * Hry[i] * Phi_row[j];
        }
    }
}
