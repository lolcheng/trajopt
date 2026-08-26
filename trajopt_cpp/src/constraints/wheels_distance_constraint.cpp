#include "constraints/wheels_distance_constraint.hpp"
#include "polynomial.hpp"
#include <cmath>

WheelsDistanceConstraint::WheelsDistanceConstraint(
    int n_coeff, int n_sample,
    double min_dist, double max_dist)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      min_dist_sq_(min_dist * min_dist),
      max_dist_sq_(max_dist * max_dist) {
    
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

WheelsDistanceConstraint::~WheelsDistanceConstraint() = default;

int WheelsDistanceConstraint::getNumConstraints() const {
    return n_sample_;
}

void WheelsDistanceConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double WheelsDistanceConstraint::computeDistSq(
    double xl, double yl, double zl,
    double xr, double yr, double zr) const {
    double dx = xr - xl;
    double dy = yr - yl;
    double dz = zr - zl;
    return dx * dx + dy * dy + dz * dz;
}

void WheelsDistanceConstraint::eval_g(const double* x, double* g) {
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
    
    // 计算距离平方（约束值 = dist²）
    for (int i = 0; i < n_sample_; ++i) {
        g[i] = computeDistSq(Xl[i], Yl[i], Zl[i], Xr[i], Yr[i], Zr[i]);
    }
}

void WheelsDistanceConstraint::eval_jac_g(const double* x, double* values, int offset) {
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
    
    // 雅可比矩阵：∂(dist²)/∂coeff
    // dist² = (Xr - Xl)² + (Yr - Yl)² + (Zr - Zl)²
    // ∂(dist²)/∂clx = -2 * (Xr - Xl) * ∂Xl/∂clx
    // ∂(dist²)/∂crx = 2 * (Xr - Xl) * ∂Xr/∂crx
    
    int idx = 0;
    (void)offset;
    
    for (int i = 0; i < n_sample_; ++i) {
        double dx = Xr[i] - Xl[i];
        double dy = Yr[i] - Yl[i];
        double dz = Zr[i] - Zl[i];
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对clx的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -2.0 * dx * Phi_row[j];
        }
        
        // 对cly的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -2.0 * dy * Phi_row[j];
        }
        
        // 对clz的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -2.0 * dz * Phi_row[j];
        }
        
        // 对crx的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dx * Phi_row[j];
        }
        
        // 对cry的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dy * Phi_row[j];
        }
        
        // 对crz的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dz * Phi_row[j];
        }
    }
}

void WheelsDistanceConstraint::getBounds(double* g_l, double* g_u, int offset) {
    // 约束：min_dist_sq ≤ dist² ≤ max_dist_sq
    (void)offset;
    for (int i = 0; i < n_sample_; ++i) {
        g_l[i] = min_dist_sq_;
        g_u[i] = max_dist_sq_;
    }
}

int WheelsDistanceConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        // 只返回非零元素数量
        // 每个约束依赖6个变量组（clx, cly, clz, crx, cry, crz），每个变量组n_coeff个系数
        return n_sample_ * 6 * n_coeff_;
    }
    
    int idx = 0;
    
    for (int i = 0; i < n_sample_; ++i) {
        int row = i + offset;
        
        // clx, cly, clz, crx, cry, crz
        int col_bases[] = {idx_clx_, idx_cly_, idx_clz_, idx_crx_, idx_cry_, idx_crz_};
        for (int var_idx = 0; var_idx < 6; ++var_idx) {
            for (int j = 0; j < n_coeff_; ++j) {
                iRow[idx] = row;
                jCol[idx] = col_bases[var_idx] + j;
                idx++;
            }
        }
    }
    
    return idx;
}

