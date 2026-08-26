#include "constraints/wheel_base_distance_constraint.hpp"
#include "polynomial.hpp"
#include <cmath>

WheelBaseDistanceConstraint::WheelBaseDistanceConstraint(
    int n_coeff, int n_sample,
    double min_dist, double max_dist)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      min_dist_sq_(min_dist * min_dist),
      max_dist_sq_(max_dist * max_dist) {
    
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

WheelBaseDistanceConstraint::~WheelBaseDistanceConstraint() = default;

int WheelBaseDistanceConstraint::getNumConstraints() const {
    // 左轮和右轮各n_sample个约束
    return 2 * n_sample_;
}

void WheelBaseDistanceConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double WheelBaseDistanceConstraint::computeDistSq(
    double xb, double yb, double zb,
    double xw, double yw, double zw) const {
    double dx = xb - xw;
    double dy = yb - yw;
    double dz = zb - zw;
    return dx * dx + dy * dy + dz * dz;
}

void WheelBaseDistanceConstraint::eval_g(const double* x, double* g) {
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
    
    // 计算距离平方（约束值 = dist²）
    // 前n_sample个：左轮距离平方
    // 后n_sample个：右轮距离平方
    for (int i = 0; i < n_sample_; ++i) {
        g[i] = computeDistSq(Xb[i], Yb[i], Zb[i], Xl[i], Yl[i], Zl[i]);
        g[i + n_sample_] = computeDistSq(Xb[i], Yb[i], Zb[i], Xr[i], Yr[i], Zr[i]);
    }
}

void WheelBaseDistanceConstraint::eval_jac_g(const double* x, double* values, int offset) {
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
    
    // 雅可比矩阵：∂(dist²)/∂coeff
    // dist² = (Xb - Xw)² + (Yb - Yw)² + (Zb - Zw)²
    // ∂(dist²)/∂cbx = 2 * (Xb - Xw) * ∂Xb/∂cbx
    // ∂(dist²)/∂clx = -2 * (Xb - Xw) * ∂Xl/∂clx
    
    int idx = 0;
    
    // 左轮约束的雅可比（前n_sample个约束）
    for (int i = 0; i < n_sample_; ++i) {
        double dx = Xb[i] - Xl[i];
        double dy = Yb[i] - Yl[i];
        double dz = Zb[i] - Zl[i];
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对cbx的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dx * Phi_row[j];
        }
        
        // 对cby的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dy * Phi_row[j];
        }
        
        // 对cbz的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dz * Phi_row[j];
        }
        
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
    }
    
    // 右轮约束的雅可比（后n_sample个约束）
    for (int i = 0; i < n_sample_; ++i) {
        double dx = Xb[i] - Xr[i];
        double dy = Yb[i] - Yr[i];
        double dz = Zb[i] - Zr[i];
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对cbx的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dx * Phi_row[j];
        }
        
        // 对cby的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dy * Phi_row[j];
        }
        
        // 对cbz的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 2.0 * dz * Phi_row[j];
        }
        
        // 对crx的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -2.0 * dx * Phi_row[j];
        }
        
        // 对cry的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -2.0 * dy * Phi_row[j];
        }
        
        // 对crz的梯度
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -2.0 * dz * Phi_row[j];
        }
    }
}

void WheelBaseDistanceConstraint::getBounds(double* g_l, double* g_u, int offset) {
    // 约束：min_dist_sq ≤ dist² ≤ max_dist_sq
    // 注意：g_l和g_u指针已经在调用时偏移了offset，所以这里不需要再使用offset
    (void)offset;  // offset参数保留用于接口一致性，但不使用
    for (int i = 0; i < 2 * n_sample_; ++i) {
        g_l[i] = min_dist_sq_;
        g_u[i] = max_dist_sq_;
    }
}

int WheelBaseDistanceConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        // 只返回非零元素数量
        return 2 * n_sample_ * 6 * n_coeff_;
    }
    
    int idx = 0;
    
    // 左轮约束的雅可比结构
    for (int i = 0; i < n_sample_; ++i) {
        int row = i + offset;
        
        // cbx, cby, cbz, clx, cly, clz
        for (int var_idx = 0; var_idx < 6; ++var_idx) {
            int col_base = (var_idx < 3) ? var_idx * n_coeff_ : (var_idx + 1) * n_coeff_;
            for (int j = 0; j < n_coeff_; ++j) {
                iRow[idx] = row;
                jCol[idx] = col_base + j;
                idx++;
            }
        }
    }
    
    // 右轮约束的雅可比结构
    for (int i = 0; i < n_sample_; ++i) {
        int row = i + n_sample_ + offset;
        
        // cbx, cby, cbz, crx, cry, crz
        for (int var_idx = 0; var_idx < 6; ++var_idx) {
            int col_base = (var_idx < 3) ? var_idx * n_coeff_ : (var_idx + 4) * n_coeff_;
            for (int j = 0; j < n_coeff_; ++j) {
                iRow[idx] = row;
                jCol[idx] = col_base + j;
                idx++;
            }
        }
    }
    
    return idx;
}

