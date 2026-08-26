#include "objectives/smoothing_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

SmoothingObjective::SmoothingObjective(int n_coeff, int n_sample, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      weight_(weight < 0 ? 1.0 / n_sample : weight),
      enabled_(true) {
    
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
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void SmoothingObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double SmoothingObjective::smoothPenalty(const double* traj) const {
    // 二阶差分: (v3 - 2*v2 + v1)²
    double penalty = 0.0;
    for (int i = 0; i < n_sample_ - 2; ++i) {
        double diff = traj[i + 2] - 2.0 * traj[i + 1] + traj[i];
        penalty += diff * diff;
    }
    return penalty;
}

// smoothPenaltyGrad方法已移除，直接在eval_grad_f中实现

double SmoothingObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取系数
    const double* coeffs[] = {
        x + idx_cbx_, x + idx_cby_, x + idx_cbz_, x + idx_cpsi_,
        x + idx_clx_, x + idx_cly_, x + idx_clz_,
        x + idx_crx_, x + idx_cry_, x + idx_crz_
    };
    
    double obj = 0.0;
    
    // 评估所有轨迹并计算平滑惩罚
    std::vector<double> traj(n_sample_);
    for (int i = 0; i < 10; ++i) {
        evalTrajectory(coeffs[i], traj.data());
        obj += smoothPenalty(traj.data());
    }
    
    return weight_ * obj;
}

void SmoothingObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取系数
    const double* coeffs[] = {
        x + idx_cbx_, x + idx_cby_, x + idx_cbz_, x + idx_cpsi_,
        x + idx_clx_, x + idx_cly_, x + idx_clz_,
        x + idx_crx_, x + idx_cry_, x + idx_crz_
    };
    
    double* grad_coeffs[] = {
        grad_f + idx_cbx_, grad_f + idx_cby_, grad_f + idx_cbz_, grad_f + idx_cpsi_,
        grad_f + idx_clx_, grad_f + idx_cly_, grad_f + idx_clz_,
        grad_f + idx_crx_, grad_f + idx_cry_, grad_f + idx_crz_
    };
    
    // 评估所有轨迹
    std::vector<std::vector<double>> trajs(10, std::vector<double>(n_sample_));
    for (int i = 0; i < 10; ++i) {
        evalTrajectory(coeffs[i], trajs[i].data());
    }
    
    // 计算梯度
    for (int i = 0; i < 10; ++i) {
        // 对每个轨迹计算平滑惩罚的梯度
        for (int j = 0; j < n_sample_; ++j) {
            double grad_traj_j = 0.0;
            
            // 如果j >= 2，贡献来自 (traj[j] - 2*traj[j-1] + traj[j-2])
            if (j >= 2) {
                double diff = trajs[i][j] - 2.0 * trajs[i][j - 1] + trajs[i][j - 2];
                grad_traj_j += 2.0 * diff;
            }
            
            // 如果j >= 1 && j < n_sample - 1，贡献来自 (traj[j+1] - 2*traj[j] + traj[j-1])
            if (j >= 1 && j < n_sample_ - 1) {
                double diff = trajs[i][j + 1] - 2.0 * trajs[i][j] + trajs[i][j - 1];
                grad_traj_j += 2.0 * diff * (-2.0);
            }
            
            // 如果j < n_sample - 2，贡献来自 (traj[j+2] - 2*traj[j+1] + traj[j])
            if (j < n_sample_ - 2) {
                double diff = trajs[i][j + 2] - 2.0 * trajs[i][j + 1] + trajs[i][j];
                grad_traj_j += 2.0 * diff;
            }
            
            // 累加到系数梯度
            const double* Phi_row = Phi_.data() + j * n_coeff_;
            for (int k = 0; k < n_coeff_; ++k) {
                grad_coeffs[i][k] += weight_ * grad_traj_j * Phi_row[k];
            }
        }
    }
}

