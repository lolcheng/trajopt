#include "objectives/regularization_objective.hpp"

RegularizationObjective::RegularizationObjective(int n_coeff, double weight)
    : n_coeff_(n_coeff), weight_(weight), enabled_(true) {
    
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

double RegularizationObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
    double obj = 0.0;
    
    // 累加所有系数的平方
    const double* coeffs[] = {
        x + idx_cbx_, x + idx_cby_, x + idx_cbz_, x + idx_cpsi_,
        x + idx_clx_, x + idx_cly_, x + idx_clz_,
        x + idx_crx_, x + idx_cry_, x + idx_crz_
    };
    
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < n_coeff_; ++j) {
            double c = coeffs[i][j];
            obj += c * c;
        }
    }
    
    return weight_ * obj;
}

void RegularizationObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
    // 梯度：2 * weight * coeff
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
    
    double factor = 2.0 * weight_;
    
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < n_coeff_; ++j) {
            grad_coeffs[i][j] += factor * coeffs[i][j];
        }
    }
}

