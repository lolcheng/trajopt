#include "objectives/force_balance_objective.hpp"
#include <cmath>

ForceBalanceObjective::ForceBalanceObjective(
    int n_coeff, int n_sample,
    double mass, double gravity, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      mg_(mass * gravity),
      weight_(weight < 0 ? 1.0 / n_sample : weight),  // 降低默认权重：10 -> 1
      enabled_(true) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    int n_poly_vars = 10 * n_coeff_;
    idx_fLx_ = n_poly_vars;
    idx_fLy_ = idx_fLx_ + n_sample_;
    idx_fLz_ = idx_fLy_ + n_sample_;
    idx_fRx_ = idx_fLz_ + n_sample_;
    idx_fRy_ = idx_fRx_ + n_sample_;
    idx_fRz_ = idx_fRy_ + n_sample_;
}

double ForceBalanceObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    // 提取接触力
    const double* fLx = x + idx_fLx_;
    const double* fLy = x + idx_fLy_;
    const double* fLz = x + idx_fLz_;
    const double* fRx = x + idx_fRx_;
    const double* fRy = x + idx_fRy_;
    const double* fRz = x + idx_fRz_;
    
    // 计算目标值: sum(||fL + fR + [0, 0, -mg]||²)
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        double Fx_res = fLx[i] + fRx[i];
        double Fy_res = fLy[i] + fRy[i];
        double Fz_res = fLz[i] + fRz[i] - mg_;
        
        obj += Fx_res * Fx_res + Fy_res * Fy_res + Fz_res * Fz_res;
    }
    
    return weight_ * obj;
}

void ForceBalanceObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 提取接触力
    const double* fLx = x + idx_fLx_;
    const double* fLy = x + idx_fLy_;
    const double* fLz = x + idx_fLz_;
    const double* fRx = x + idx_fRx_;
    const double* fRy = x + idx_fRy_;
    const double* fRz = x + idx_fRz_;
    
    double* grad_fLx = grad_f + idx_fLx_;
    double* grad_fLy = grad_f + idx_fLy_;
    double* grad_fLz = grad_f + idx_fLz_;
    double* grad_fRx = grad_f + idx_fRx_;
    double* grad_fRy = grad_f + idx_fRy_;
    double* grad_fRz = grad_f + idx_fRz_;
    
    // 计算梯度
    // d/dcoeff ||fL + fR + [0, 0, -mg]||² = 2 * (fL + fR + [0, 0, -mg])
    for (int i = 0; i < n_sample_; ++i) {
        double Fx_res = fLx[i] + fRx[i];
        double Fy_res = fLy[i] + fRy[i];
        double Fz_res = fLz[i] + fRz[i] - mg_;
        
        double factor = 2.0 * weight_;
        
        grad_fLx[i] += factor * Fx_res;
        grad_fLy[i] += factor * Fy_res;
        grad_fLz[i] += factor * Fz_res;
        grad_fRx[i] += factor * Fx_res;
        grad_fRy[i] += factor * Fy_res;
        grad_fRz[i] += factor * Fz_res;
    }
}
