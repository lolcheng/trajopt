#include "objectives/yaw_path_alignment_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

YawPathAlignmentObjective::YawPathAlignmentObjective(int n_coeff, int n_sample, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      weight_(weight < 0 ? 50.0 / n_sample : weight),
      enabled_(true) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

void YawPathAlignmentObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void YawPathAlignmentObjective::computePathTangent(
    const std::vector<double>& Xb, const std::vector<double>& Yb,
    std::vector<double>& tx, std::vector<double>& ty) const {
    const double eps_norm = 1e-6;
    
    for (int i = 0; i < n_sample_ - 1; ++i) {
        double dx = Xb[i + 1] - Xb[i];
        double dy = Yb[i + 1] - Yb[i];
        double norm_t_sq = dx * dx + dy * dy + eps_norm;
        double norm_t = std::sqrt(norm_t_sq);
        // 防止除零
        if (norm_t < 1e-10) {
            norm_t = 1e-10;
        }
        tx[i] = dx / norm_t;
        ty[i] = dy / norm_t;
    }
    // 最后一个点使用前一个点的切线
    if (n_sample_ > 1) {
        tx[n_sample_ - 1] = tx[n_sample_ - 2];
        ty[n_sample_ - 1] = ty[n_sample_ - 2];
    }
}

double YawPathAlignmentObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cpsi = x + idx_cpsi_;
    
    // 评估轨迹
    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    std::vector<double> Psi(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cpsi, Psi.data());
    
    // 计算路径切线方向
    std::vector<double> tx(n_sample_), ty(n_sample_);
    computePathTangent(Xb, Yb, tx, ty);
    
    // 计算目标值: sum((1 - cos(alignment))²)
    // alignment = arccos(tx * cos(psi) + ty * sin(psi))
    // cos(alignment) = tx * cos(psi) + ty * sin(psi)
    double obj = 0.0;
    for (int i = 0; i < n_sample_ - 1; ++i) {
        double psi_mid = 0.5 * (Psi[i] + Psi[i + 1]);
        double cos_psi = std::cos(psi_mid);
        double sin_psi = std::sin(psi_mid);
        
        double cos_align = tx[i] * cos_psi + ty[i] * sin_psi;
        double yaw_misalign = 1.0 - cos_align;
        obj += yaw_misalign * yaw_misalign;
    }
    
    return weight_ * obj;
}

void YawPathAlignmentObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cpsi = x + idx_cpsi_;
    
    // 评估轨迹
    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    std::vector<double> Psi(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cpsi, Psi.data());
    
    // 计算路径切线方向
    std::vector<double> tx(n_sample_), ty(n_sample_);
    computePathTangent(Xb, Yb, tx, ty);
    
    double* grad_cbx = grad_f + idx_cbx_;
    double* grad_cby = grad_f + idx_cby_;
    double* grad_cpsi = grad_f + idx_cpsi_;
    
    // 计算梯度
    // d/dcoeff (1 - cos(alignment))² = 2 * (1 - cos(alignment)) * (-d(cos(alignment))/dcoeff)
    // cos(alignment) = tx * cos(psi) + ty * sin(psi)
    // d(cos(alignment))/dpsi = -tx * sin(psi) + ty * cos(psi)
    
    for (int i = 0; i < n_sample_ - 1; ++i) {
        double psi_mid = 0.5 * (Psi[i] + Psi[i + 1]);
        double cos_psi = std::cos(psi_mid);
        double sin_psi = std::sin(psi_mid);
        
        double cos_align = tx[i] * cos_psi + ty[i] * sin_psi;
        double yaw_misalign = 1.0 - cos_align;
        double factor = -2.0 * weight_ * yaw_misalign;
        
        // 对psi的梯度
        double dcos_align_dpsi = -tx[i] * sin_psi + ty[i] * cos_psi;
        const double* Phi_row_i = Phi_.data() + i * n_coeff_;
        const double* Phi_row_i1 = Phi_.data() + (i + 1) * n_coeff_;
        
        for (int j = 0; j < n_coeff_; ++j) {
            // psi_mid = 0.5 * (Psi[i] + Psi[i+1])
            // d(psi_mid)/dPsi[i] = 0.5
            grad_cpsi[j] += factor * dcos_align_dpsi * 0.5 * Phi_row_i[j];
            grad_cpsi[j] += factor * dcos_align_dpsi * 0.5 * Phi_row_i1[j];
        }
        
        // 对Xb, Yb的梯度（通过tx, ty）
        // tx = dx / norm_t, ty = dy / norm_t
        // dx = Xb[i+1] - Xb[i], dy = Yb[i+1] - Yb[i]
        // norm_t = sqrt(dx² + dy² + eps)
        double dx = Xb[i + 1] - Xb[i];
        double dy = Yb[i + 1] - Yb[i];
        double norm_t_sq = dx * dx + dy * dy + 1e-6;
        double norm_t = std::sqrt(norm_t_sq);
        // 防止除零
        if (norm_t < 1e-10) {
            norm_t = 1e-10;
            norm_t_sq = norm_t * norm_t;
        }
        double inv_norm_t = 1.0 / norm_t;
        double inv_norm_t3 = inv_norm_t / norm_t_sq;
        
        // tx = dx / norm_t, ty = dy / norm_t
        // d(tx)/dXb[i] = -1/norm_t + dx²/norm_t³
        // d(tx)/dXb[i+1] = 1/norm_t - dx²/norm_t³
        // d(ty)/dYb[i] = -1/norm_t + dy²/norm_t³
        // d(ty)/dYb[i+1] = 1/norm_t - dy²/norm_t³
        double dcos_align_dtx = cos_psi;
        double dcos_align_dty = sin_psi;
        
        double dtx_dxb_i = -inv_norm_t + dx * dx * inv_norm_t3;
        double dtx_dxb_i1 = inv_norm_t - dx * dx * inv_norm_t3;
        double dty_dyb_i = -inv_norm_t + dy * dy * inv_norm_t3;
        double dty_dyb_i1 = inv_norm_t - dy * dy * inv_norm_t3;
        
        // 注意：tx也依赖于Yb（通过norm_t），ty也依赖于Xb（通过norm_t）
        double dtx_dyb_i = dx * dy * inv_norm_t3;
        double dtx_dyb_i1 = -dx * dy * inv_norm_t3;
        double dty_dxb_i = dx * dy * inv_norm_t3;
        double dty_dxb_i1 = -dx * dy * inv_norm_t3;
        
        for (int j = 0; j < n_coeff_; ++j) {
            // 对cbx的梯度（通过tx和ty）
            grad_cbx[j] += factor * (dcos_align_dtx * dtx_dxb_i + dcos_align_dty * dty_dxb_i) * Phi_row_i[j];
            grad_cbx[j] += factor * (dcos_align_dtx * dtx_dxb_i1 + dcos_align_dty * dty_dxb_i1) * Phi_row_i1[j];
            // 对cby的梯度（通过tx和ty）
            grad_cby[j] += factor * (dcos_align_dtx * dtx_dyb_i + dcos_align_dty * dty_dyb_i) * Phi_row_i[j];
            grad_cby[j] += factor * (dcos_align_dtx * dtx_dyb_i1 + dcos_align_dty * dty_dyb_i1) * Phi_row_i1[j];
        }
    }
}

