#include "objectives/path_length_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

PathLengthObjective::PathLengthObjective(int n_coeff, int n_sample, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      weight_(weight < 0 ? 1.0 / n_sample : weight),  // 默认权重
      enabled_(true) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
    
    // 预计算导数基矩阵（用于计算 dx/ds）
    // 注意：这里不能写死为幂基导数 [0,1,2s,...]，否则在切比雪夫基下会错误。
    // 采用通用方式：对每个单位系数向量 e_j，调用 Polynomial::eval_derivative
    // 得到第 j 列的导数基值，从而兼容任意 BasisFunction。
    Phi_dot_.resize(n_sample_ * n_coeff_);
    std::vector<double> unit_coeff(n_coeff_, 0.0);
    std::vector<double> col_deriv(n_sample_, 0.0);
    for (int j = 0; j < n_coeff_; ++j) {
        std::fill(unit_coeff.begin(), unit_coeff.end(), 0.0);
        unit_coeff[j] = 1.0;
        Polynomial::eval_derivative(s_samples_.data(), n_sample_, unit_coeff.data(), n_coeff_, col_deriv.data());
        for (int i = 0; i < n_sample_; ++i) {
            Phi_dot_[i * n_coeff_ + j] = col_deriv[i];
        }
    }
}

void PathLengthObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void PathLengthObjective::evalTrajectoryDerivative(const double* coeff, double* traj_dot) const {
    // 计算轨迹导数：traj_dot = Phi_dot * coeff
    for (int i = 0; i < n_sample_; ++i) {
        const double* Phi_dot_row = Phi_dot_.data() + i * n_coeff_;
        double sum = 0.0;
        for (int j = 0; j < n_coeff_; ++j) {
            sum += Phi_dot_row[j] * coeff[j];
        }
        traj_dot[i] = sum;
    }
}

double PathLengthObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取base轨迹系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    
    // 评估轨迹和导数
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Zb(n_sample_);
    std::vector<double> Xb_dot(n_sample_), Yb_dot(n_sample_), Zb_dot(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cbz, Zb.data());
    evalTrajectoryDerivative(cbx, Xb_dot.data());
    evalTrajectoryDerivative(cby, Yb_dot.data());
    evalTrajectoryDerivative(cbz, Zb_dot.data());
    
    // 计算轨迹长度（使用梯形法则）
    // L = ∫ sqrt((dx/ds)² + (dy/ds)² + (dz/ds)²) ds
    // 在采样点上使用梯形法则积分
    double length = 0.0;
    double ds = 1.0 / (n_sample_ - 1);  // 采样点间距
    
    for (int i = 0; i < n_sample_ - 1; ++i) {
        // 计算当前点和下一个点的速度大小
        double v1 = std::sqrt(Xb_dot[i] * Xb_dot[i] + 
                              Yb_dot[i] * Yb_dot[i] + 
                              Zb_dot[i] * Zb_dot[i]);
        double v2 = std::sqrt(Xb_dot[i + 1] * Xb_dot[i + 1] + 
                              Yb_dot[i + 1] * Yb_dot[i + 1] + 
                              Zb_dot[i + 1] * Zb_dot[i + 1]);
        
        // 梯形法则：∫ f(s) ds ≈ (f(s1) + f(s2)) / 2 * ds
        length += 0.5 * (v1 + v2) * ds;
    }
    
    return weight_ * length;
}

void PathLengthObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取base轨迹系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    
    // 评估轨迹和导数
    std::vector<double> Xb(n_sample_), Yb(n_sample_), Zb(n_sample_);
    std::vector<double> Xb_dot(n_sample_), Yb_dot(n_sample_), Zb_dot(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cbz, Zb.data());
    evalTrajectoryDerivative(cbx, Xb_dot.data());
    evalTrajectoryDerivative(cby, Yb_dot.data());
    evalTrajectoryDerivative(cbz, Zb_dot.data());
    
    // 计算速度大小
    std::vector<double> v(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        v[i] = std::sqrt(Xb_dot[i] * Xb_dot[i] + 
                        Yb_dot[i] * Yb_dot[i] + 
                        Zb_dot[i] * Zb_dot[i]);
    }
    
    // 计算梯度
    double* grad_cbx = grad_f + idx_cbx_;
    double* grad_cby = grad_f + idx_cby_;
    double* grad_cbz = grad_f + idx_cbz_;
    
    double ds = 1.0 / (n_sample_ - 1);
    
    for (int i = 0; i < n_sample_; ++i) {
        const double* Phi_dot_row = Phi_dot_.data() + i * n_coeff_;
        
        // 计算 ∂L/∂coeff
        // L = ∫ sqrt((dx/ds)² + (dy/ds)² + (dz/ds)²) ds
        // ∂L/∂cbx[j] = ∫ (dx/ds / v) * ∂(dx/ds)/∂cbx[j] ds
        //            = ∫ (dx/ds / v) * Phi_dot[i][j] ds
        
        double factor = 0.0;
        if (i == 0) {
            // 第一个点：只贡献到第一个梯形
            if (v[i] > 1e-10) {  // 避免除零
                factor = 0.5 * ds * weight_ / v[i];
            }
        } else if (i == n_sample_ - 1) {
            // 最后一个点：只贡献到最后一个梯形
            if (v[i] > 1e-10) {
                factor = 0.5 * ds * weight_ / v[i];
            }
        } else {
            // 中间点：贡献到两个梯形
            if (v[i] > 1e-10) {
                factor = ds * weight_ / v[i];
            }
        }
        
        // 对cbx的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cbx[j] += factor * Xb_dot[i] * Phi_dot_row[j];
        }
        
        // 对cby的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cby[j] += factor * Yb_dot[i] * Phi_dot_row[j];
        }
        
        // 对cbz的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            grad_cbz[j] += factor * Zb_dot[i] * Phi_dot_row[j];
        }
    }
}

