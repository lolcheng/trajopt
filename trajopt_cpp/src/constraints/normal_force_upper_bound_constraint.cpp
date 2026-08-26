#include "constraints/normal_force_upper_bound_constraint.hpp"
#include "polynomial.hpp"
#include <cmath>

NormalForceUpperBoundConstraint::NormalForceUpperBoundConstraint(
    int n_coeff, int n_sample,
    std::shared_ptr<RBFTerrain> terrain,
    double mass, double gravity, double max_accel)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      fN_max_(mass * (gravity + max_accel)),
      terrain_(terrain) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    int n_poly_vars = 10 * n_coeff_;
    idx_fLx_ = n_poly_vars;
    idx_fLy_ = idx_fLx_ + n_sample_;
    idx_fLz_ = idx_fLy_ + n_sample_;
    idx_fRx_ = idx_fLz_ + n_sample_;
    idx_fRy_ = idx_fRx_ + n_sample_;
    idx_fRz_ = idx_fRy_ + n_sample_;
    
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

NormalForceUpperBoundConstraint::~NormalForceUpperBoundConstraint() = default;

int NormalForceUpperBoundConstraint::getNumConstraints() const {
    return 2 * n_sample_;
}

void NormalForceUpperBoundConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double NormalForceUpperBoundConstraint::computeNormalForce(
    double fx, double fy, double fz,
    double nx, double ny, double nz) const {
    return fx * nx + fy * ny + fz * nz;
}

void NormalForceUpperBoundConstraint::eval_g(const double* x, double* g) {
    // 提取接触力
    const double* fLx = x + idx_fLx_;
    const double* fLy = x + idx_fLy_;
    const double* fLz = x + idx_fLz_;
    const double* fRx = x + idx_fRx_;
    const double* fRy = x + idx_fRy_;
    const double* fRz = x + idx_fRz_;
    
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
    
    // 计算约束值：g = fN（法向力）
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形法向量
        double nx_l, ny_l, nz_l;
        double nx_r, ny_r, nz_r;
        terrain_->normal(Xl[i], Yl[i], nx_l, ny_l, nz_l);
        terrain_->normal(Xr[i], Yr[i], nx_r, ny_r, nz_r);
        
        // 计算法向力
        g[i] = computeNormalForce(fLx[i], fLy[i], fLz[i], nx_l, ny_l, nz_l);
        g[i + n_sample_] = computeNormalForce(fRx[i], fRy[i], fRz[i], nx_r, ny_r, nz_r);
    }
}

void NormalForceUpperBoundConstraint::eval_jac_g(const double* x, double* values, int offset) {
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
    
    // 雅可比矩阵：∂fN/∂coeff
    // fN = f · n = fx*nx + fy*ny + fz*nz
    // ∂fN/∂fx = nx, ∂fN/∂fy = ny, ∂fN/∂fz = nz
    // ∂fN/∂Xl = fx*∂nx/∂Xl + fy*∂ny/∂Xl + fz*∂nz/∂Xl (通过地形法向量)
    
    int idx = 0;
    (void)offset;
    
    // 左轮约束（g[0] 到 g[n_sample-1]）
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形法向量
        double nx_l, ny_l, nz_l;
        terrain_->normal(Xl[i], Yl[i], nx_l, ny_l, nz_l);
        
        // 对接触力的梯度
        values[idx++] = nx_l;  // ∂fN/∂fLx
        values[idx++] = ny_l;  // ∂fN/∂fLy
        values[idx++] = nz_l;  // ∂fN/∂fLz
        
        // 对轮子位置的梯度（通过法向量，这里简化处理，只考虑法向量对位置的依赖）
        // 注意：法向量的梯度计算较复杂，这里先简化，假设法向量对位置变化不敏感
        // 实际应用中，可以通过数值微分或解析计算法向量的梯度
        // 暂时设为0（如果地形变化不大，这个近似是合理的）
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 0.0;  // ∂fN/∂clx (简化)
            values[idx++] = 0.0;  // ∂fN/∂cly (简化)
            values[idx++] = 0.0;  // ∂fN/∂clz (简化)
        }
    }
    
    // 右轮约束（g[n_sample] 到 g[2*n_sample-1]）
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形法向量
        double nx_r, ny_r, nz_r;
        terrain_->normal(Xr[i], Yr[i], nx_r, ny_r, nz_r);
        
        // 对接触力的梯度
        values[idx++] = nx_r;  // ∂fN/∂fRx
        values[idx++] = ny_r;  // ∂fN/∂fRy
        values[idx++] = nz_r;  // ∂fN/∂fRz
        
        // 对轮子位置的梯度（简化处理）
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = 0.0;  // ∂fN/∂crx (简化)
            values[idx++] = 0.0;  // ∂fN/∂cry (简化)
            values[idx++] = 0.0;  // ∂fN/∂crz (简化)
        }
    }
}

void NormalForceUpperBoundConstraint::getBounds(double* g_l, double* g_u, int offset) {
    // 约束：fN ≤ fN_max
    (void)offset;
    for (int i = 0; i < 2 * n_sample_; ++i) {
        g_l[i] = -1e20;  // 无下界
        g_u[i] = fN_max_;
    }
}

int NormalForceUpperBoundConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        // 只返回非零元素数量
        // 每个约束依赖3个接触力变量 + 3个轮子位置变量组（每个n_coeff个系数）
        return 2 * n_sample_ * (3 + 3 * n_coeff_);
    }
    
    int idx = 0;
    
    // 左轮约束
    for (int i = 0; i < n_sample_; ++i) {
        int row = i + offset;
        
        // 接触力变量
        iRow[idx] = row; jCol[idx] = idx_fLx_ + i; idx++;
        iRow[idx] = row; jCol[idx] = idx_fLy_ + i; idx++;
        iRow[idx] = row; jCol[idx] = idx_fLz_ + i; idx++;
        
        // 轮子位置变量（简化，实际应该考虑法向量对位置的依赖）
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row; jCol[idx] = idx_clx_ + j; idx++;
            iRow[idx] = row; jCol[idx] = idx_cly_ + j; idx++;
            iRow[idx] = row; jCol[idx] = idx_clz_ + j; idx++;
        }
    }
    
    // 右轮约束
    for (int i = 0; i < n_sample_; ++i) {
        int row = i + n_sample_ + offset;
        
        // 接触力变量
        iRow[idx] = row; jCol[idx] = idx_fRx_ + i; idx++;
        iRow[idx] = row; jCol[idx] = idx_fRy_ + i; idx++;
        iRow[idx] = row; jCol[idx] = idx_fRz_ + i; idx++;
        
        // 轮子位置变量
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[idx] = row; jCol[idx] = idx_crx_ + j; idx++;
            iRow[idx] = row; jCol[idx] = idx_cry_ + j; idx++;
            iRow[idx] = row; jCol[idx] = idx_crz_ + j; idx++;
        }
    }
    
    return idx;
}

