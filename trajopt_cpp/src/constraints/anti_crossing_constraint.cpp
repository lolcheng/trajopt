#include "constraints/anti_crossing_constraint.hpp"
#include "polynomial.hpp"
#include <cmath>

AntiCrossingConstraint::AntiCrossingConstraint(int n_coeff, int n_sample)
    : n_coeff_(n_coeff), n_sample_(n_sample) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_crx_ = 7 * n_coeff_;
    idx_cry_ = 8 * n_coeff_;
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

AntiCrossingConstraint::~AntiCrossingConstraint() = default;

int AntiCrossingConstraint::getNumConstraints() const {
    return n_sample_;
}

void AntiCrossingConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void AntiCrossingConstraint::eval_g(const double* x, double* g) {
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    
    // 评估轨迹
    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    std::vector<double> Xl(n_sample_), Yl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_);
    std::vector<double> Psi(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cpsi, Psi.data());
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    
    // 计算约束值：g = yL_base - yR_base
    // yL_base = (Xl - Xb) · u_lat_x + (Yl - Yb) · u_lat_y
    // yR_base = (Xr - Xb) · u_lat_x + (Yr - Yb) · u_lat_y
    // u_lat = [-sin(psi), cos(psi)]
    for (int i = 0; i < n_sample_; ++i) {
        double u_lat_x = -std::sin(Psi[i]);
        double u_lat_y = std::cos(Psi[i]);
        
        double vLx = Xl[i] - Xb[i];
        double vLy = Yl[i] - Yb[i];
        double vRx = Xr[i] - Xb[i];
        double vRy = Yr[i] - Yb[i];
        
        double yL_base = vLx * u_lat_x + vLy * u_lat_y;
        double yR_base = vRx * u_lat_x + vRy * u_lat_y;
        
        g[i] = yL_base - yR_base;  // 约束：g ≥ 0
    }
}

void AntiCrossingConstraint::eval_jac_g(const double* x, double* values, int offset) {
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    
    // 评估轨迹
    std::vector<double> Xb(n_sample_), Yb(n_sample_);
    std::vector<double> Xl(n_sample_), Yl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_);
    std::vector<double> Psi(n_sample_);
    
    evalTrajectory(cbx, Xb.data());
    evalTrajectory(cby, Yb.data());
    evalTrajectory(cpsi, Psi.data());
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    
    // 雅可比矩阵：∂g/∂coeff
    // g = yL_base - yR_base
    // yL_base = (Xl - Xb) · u_lat_x + (Yl - Yb) · u_lat_y
    // yR_base = (Xr - Xb) · u_lat_x + (Yr - Yb) · u_lat_y
    // u_lat = [-sin(psi), cos(psi)]
    
    int idx = 0;
    (void)offset;
    
    for (int i = 0; i < n_sample_; ++i) {
        double psi = Psi[i];
        double cos_psi = std::cos(psi);
        double sin_psi = std::sin(psi);
        double u_lat_x = -sin_psi;
        double u_lat_y = cos_psi;
        
        double vLx = Xl[i] - Xb[i];
        double vLy = Yl[i] - Yb[i];
        double vRx = Xr[i] - Xb[i];
        double vRy = Yr[i] - Yb[i];
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // ∂g/∂cbx = -u_lat_x (通过yL_base和yR_base)
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -u_lat_x * Phi_row[j];
        }
        
        // ∂g/∂cby = -u_lat_y
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -u_lat_y * Phi_row[j];
        }
        
        // ∂g/∂cpsi = (vLx * (-cos_psi) + vLy * (-sin_psi)) - (vRx * (-cos_psi) + vRy * (-sin_psi))
        //          = -(vLx - vRx) * cos_psi - (vLy - vRy) * sin_psi
        double dg_dpsi = -(vLx - vRx) * cos_psi - (vLy - vRy) * sin_psi;
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = dg_dpsi * Phi_row[j];
        }
        
        // ∂g/∂clx = u_lat_x
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = u_lat_x * Phi_row[j];
        }
        
        // ∂g/∂cly = u_lat_y
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = u_lat_y * Phi_row[j];
        }
        
        // ∂g/∂crx = -u_lat_x
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -u_lat_x * Phi_row[j];
        }
        
        // ∂g/∂cry = -u_lat_y
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -u_lat_y * Phi_row[j];
        }
    }
}

void AntiCrossingConstraint::getBounds(double* g_l, double* g_u, int offset) {
    // 约束：g = yL_base - yR_base ≥ 0
    (void)offset;
    for (int i = 0; i < n_sample_; ++i) {
        g_l[i] = 0.0;
        g_u[i] = 1e20;  // 无上界
    }
}

int AntiCrossingConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        // 只返回非零元素数量
        // 每个约束依赖7个变量组（cbx, cby, cpsi, clx, cly, crx, cry），每个变量组n_coeff个系数
        return n_sample_ * 7 * n_coeff_;
    }
    
    int idx = 0;
    
    for (int i = 0; i < n_sample_; ++i) {
        int row = i + offset;
        
        // cbx, cby, cpsi, clx, cly, crx, cry
        int col_bases[] = {idx_cbx_, idx_cby_, idx_cpsi_, idx_clx_, idx_cly_, idx_crx_, idx_cry_};
        for (int var_idx = 0; var_idx < 7; ++var_idx) {
            for (int j = 0; j < n_coeff_; ++j) {
                iRow[idx] = row;
                jCol[idx] = col_bases[var_idx] + j;
                idx++;
            }
        }
    }
    
    return idx;
}

