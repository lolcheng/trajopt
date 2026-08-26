#include "objectives/wheel_airborne_penalty_objective.hpp"
#include "polynomial.hpp"
#include <cmath>
#include <algorithm>

WheelAirbornePenaltyObjective::WheelAirbornePenaltyObjective(
    int n_coeff, int n_sample,
    std::shared_ptr<RBFTerrain> terrain,
    double weight,
    double threshold)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      terrain_(terrain),
      weight_(weight < 0 ? 10.0 / n_sample : weight),  // 默认权重
      threshold_(threshold),  // 离地判据阈值（默认5cm）
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

void WheelAirbornePenaltyObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double WheelAirbornePenaltyObjective::smoothPos(double x) const {
    // smooth_pos(x) = max(0, x) 的平滑版本
    // 使用: x > 0 ? x : 0
    return x > 0.0 ? x : 0.0;
}

double WheelAirbornePenaltyObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取轮子轨迹系数
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轮子轨迹
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
    terrain_->height(Xl.data(), Yl.data(), n_sample_, Hl.data());
    terrain_->height(Xr.data(), Yr.data(), n_sample_, Hr.data());
    
    // 计算目标值: sum(max(0, Zl - Hl - threshold)² + max(0, Zr - Hr - threshold)²)
    // 只有当离地距离超过阈值时才惩罚
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        double airborne_l = smoothPos(Zl[i] - Hl[i] - threshold_);  // 左轮离地距离（超过阈值部分）
        double airborne_r = smoothPos(Zr[i] - Hr[i] - threshold_);  // 右轮离地距离（超过阈值部分）
        obj += airborne_l * airborne_l + airborne_r * airborne_r;
    }
    
    return weight_ * obj;
}

void WheelAirbornePenaltyObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取轮子轨迹系数
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轮子轨迹
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
    
    terrain_->height(Xl.data(), Yl.data(), n_sample_, Hl.data());
    terrain_->height(Xr.data(), Yr.data(), n_sample_, Hr.data());
    terrain_->gradient(Xl.data(), Yl.data(), n_sample_, Hlx.data(), Hly.data());
    terrain_->gradient(Xr.data(), Yr.data(), n_sample_, Hrx.data(), Hry.data());
    
    // 计算梯度
    double* grad_clx = grad_f + idx_clx_;
    double* grad_cly = grad_f + idx_cly_;
    double* grad_clz = grad_f + idx_clz_;
    double* grad_crx = grad_f + idx_crx_;
    double* grad_cry = grad_f + idx_cry_;
    double* grad_crz = grad_f + idx_crz_;
    
    for (int i = 0; i < n_sample_; ++i) {
        double airborne_l = smoothPos(Zl[i] - Hl[i]);  // 左轮离地距离
        double airborne_r = smoothPos(Zr[i] - Hr[i]);  // 右轮离地距离
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对左轮的梯度
        if (airborne_l > 0.0) {
            // ∂(max(0, Zl - Hl)²)/∂clz = 2 * max(0, Zl - Hl) * Phi_row
            // ∂(max(0, Zl - Hl)²)/∂clx = -2 * max(0, Zl - Hl) * Hlx * Phi_row
            // ∂(max(0, Zl - Hl)²)/∂cly = -2 * max(0, Zl - Hl) * Hly * Phi_row
            double factor = 2.0 * weight_ * airborne_l;
            for (int j = 0; j < n_coeff_; ++j) {
                grad_clz[j] += factor * Phi_row[j];
                grad_clx[j] -= factor * Hlx[i] * Phi_row[j];
                grad_cly[j] -= factor * Hly[i] * Phi_row[j];
            }
        }
        
        // 对右轮的梯度
        if (airborne_r > 0.0) {
            double factor = 2.0 * weight_ * airborne_r;
            for (int j = 0; j < n_coeff_; ++j) {
                grad_crz[j] += factor * Phi_row[j];
                grad_crx[j] -= factor * Hrx[i] * Phi_row[j];
                grad_cry[j] -= factor * Hry[i] * Phi_row[j];
            }
        }
    }
}

