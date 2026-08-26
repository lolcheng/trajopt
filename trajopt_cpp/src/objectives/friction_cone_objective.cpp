#include "objectives/friction_cone_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

FrictionConeObjective::FrictionConeObjective(
    int n_coeff, int n_sample,
    std::shared_ptr<RBFTerrain> terrain,
    double mu_friction, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      mu_friction_(mu_friction),
      mu2_(mu_friction * mu_friction),
      weight_(weight < 0 ? 10.0 / n_sample : weight),  // 进一步降低默认权重：1000 -> 100 -> 10
      enabled_(true), terrain_(terrain) {
    
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

void FrictionConeObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double FrictionConeObjective::smoothPos(double x) const {
    // smooth_pos(x) = max(0, x) 的平滑版本
    // 使用: x > 0 ? x : 0 (简化，实际可以用更平滑的函数)
    return x > 0.0 ? x : 0.0;
}

double FrictionConeObjective::computeNormalForce(
    double fx, double fy, double fz,
    double nx, double ny, double nz) const {
    return fx * nx + fy * ny + fz * nz;
}

double FrictionConeObjective::computeTangentialForceSq(
    double fx, double fy, double fz,
    double f_n, double nx, double ny, double nz) const {
    // 切向力: f_t = f - f_n * n
    double f_tx = fx - f_n * nx;
    double f_ty = fy - f_n * ny;
    double f_tz = fz - f_n * nz;
    return f_tx * f_tx + f_ty * f_ty + f_tz * f_tz;
}

double FrictionConeObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    
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
    
    // 计算目标值
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形法向量
        double nx_l, ny_l, nz_l;
        double nx_r, ny_r, nz_r;
        terrain_->normal(Xl[i], Yl[i], nx_l, ny_l, nz_l);
        terrain_->normal(Xr[i], Yr[i], nx_r, ny_r, nz_r);
        
        // 计算法向力
        double fLn = computeNormalForce(fLx[i], fLy[i], fLz[i], nx_l, ny_l, nz_l);
        double fRn = computeNormalForce(fRx[i], fRy[i], fRz[i], nx_r, ny_r, nz_r);
        
        // 计算切向力平方
        double fLt_sq = computeTangentialForceSq(fLx[i], fLy[i], fLz[i], fLn, nx_l, ny_l, nz_l);
        double fRt_sq = computeTangentialForceSq(fRx[i], fRy[i], fRz[i], fRn, nx_r, ny_r, nz_r);
        
        // 摩擦锥违反: max(0, ||f_t||² - μ² × f_n²)
        double viol_fric_L = smoothPos(fLt_sq - mu2_ * fLn * fLn);
        double viol_fric_R = smoothPos(fRt_sq - mu2_ * fRn * fRn);
        
        // 拉力违反: max(0, -f_n)
        double viol_normal_L = smoothPos(-fLn);
        double viol_normal_R = smoothPos(-fRn);
        
        obj += viol_fric_L * viol_fric_L + viol_fric_R * viol_fric_R +
               viol_normal_L * viol_normal_L + viol_normal_R * viol_normal_R;
    }
    
    return weight_ * obj;
}

void FrictionConeObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    
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
    
    double* grad_fLx = grad_f + idx_fLx_;
    double* grad_fLy = grad_f + idx_fLy_;
    double* grad_fLz = grad_f + idx_fLz_;
    double* grad_fRx = grad_f + idx_fRx_;
    double* grad_fRy = grad_f + idx_fRy_;
    double* grad_fRz = grad_f + idx_fRz_;
    // 轮子位置的梯度暂时设为0（通过法向量的梯度较复杂，简化处理）
    (void)grad_f;  // 避免未使用变量警告
    (void)idx_clx_; (void)idx_cly_; (void)idx_clz_;
    (void)idx_crx_; (void)idx_cry_; (void)idx_crz_;
    
    // 计算梯度
    // d/dcoeff max(0, ||f_t||² - μ² × f_n²)² = 2 * max(0, ||f_t||² - μ² × f_n²) * d(||f_t||² - μ² × f_n²)/dcoeff
    // d/dcoeff max(0, -f_n)² = 2 * max(0, -f_n) * d(-f_n)/dcoeff
    
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形法向量
        double nx_l, ny_l, nz_l;
        double nx_r, ny_r, nz_r;
        terrain_->normal(Xl[i], Yl[i], nx_l, ny_l, nz_l);
        terrain_->normal(Xr[i], Yr[i], nx_r, ny_r, nz_r);
        
        // 计算法向力
        double fLn = computeNormalForce(fLx[i], fLy[i], fLz[i], nx_l, ny_l, nz_l);
        double fRn = computeNormalForce(fRx[i], fRy[i], fRz[i], nx_r, ny_r, nz_r);
        
        // 计算切向力平方
        double fLt_sq = computeTangentialForceSq(fLx[i], fLy[i], fLz[i], fLn, nx_l, ny_l, nz_l);
        double fRt_sq = computeTangentialForceSq(fRx[i], fRy[i], fRz[i], fRn, nx_r, ny_r, nz_r);
        
        // 摩擦锥违反和拉力违反
        double viol_fric_L_val = fLt_sq - mu2_ * fLn * fLn;
        double viol_fric_R_val = fRt_sq - mu2_ * fRn * fRn;
        double viol_fric_L = smoothPos(viol_fric_L_val);
        double viol_fric_R = smoothPos(viol_fric_R_val);
        double viol_normal_L = smoothPos(-fLn);
        double viol_normal_R = smoothPos(-fRn);
        
        double factor_fric_L = 2.0 * weight_ * viol_fric_L;
        double factor_fric_R = 2.0 * weight_ * viol_fric_R;
        double factor_normal_L = 2.0 * weight_ * viol_normal_L;
        double factor_normal_R = 2.0 * weight_ * viol_normal_R;
        
        // 计算切向力
        double fLtx = fLx[i] - fLn * nx_l;
        double fLty = fLy[i] - fLn * ny_l;
        double fLtz = fLz[i] - fLn * nz_l;
        double fRtx = fRx[i] - fRn * nx_r;
        double fRty = fRy[i] - fRn * ny_r;
        double fRtz = fRz[i] - fRn * nz_r;
        
        // 对接触力的梯度
        // d(||f_t||²)/df = 2 * f_t
        // d(||f_t||²)/df_n = -2 * f_t · n = -2 * (f - f_n * n) · n = -2 * (f · n - f_n) = 0
        // d(||f_t||² - μ² × f_n²)/df = 2 * f_t - 2 * μ² * f_n * n
        // d(-f_n)/df = -n
        
        if (viol_fric_L_val > 0.0) {
            grad_fLx[i] += factor_fric_L * (2.0 * fLtx - 2.0 * mu2_ * fLn * nx_l);
            grad_fLy[i] += factor_fric_L * (2.0 * fLty - 2.0 * mu2_ * fLn * ny_l);
            grad_fLz[i] += factor_fric_L * (2.0 * fLtz - 2.0 * mu2_ * fLn * nz_l);
        }
        
        if (viol_fric_R_val > 0.0) {
            grad_fRx[i] += factor_fric_R * (2.0 * fRtx - 2.0 * mu2_ * fRn * nx_r);
            grad_fRy[i] += factor_fric_R * (2.0 * fRty - 2.0 * mu2_ * fRn * ny_r);
            grad_fRz[i] += factor_fric_R * (2.0 * fRtz - 2.0 * mu2_ * fRn * nz_r);
        }
        
        if (viol_normal_L > 0.0) {
            grad_fLx[i] += factor_normal_L * (-nx_l);
            grad_fLy[i] += factor_normal_L * (-ny_l);
            grad_fLz[i] += factor_normal_L * (-nz_l);
        }
        
        if (viol_normal_R > 0.0) {
            grad_fRx[i] += factor_normal_R * (-nx_r);
            grad_fRy[i] += factor_normal_R * (-ny_r);
            grad_fRz[i] += factor_normal_R * (-nz_r);
        }
        
        // 对轮子位置的梯度（通过法向量，这里简化处理，设为0）
        // 实际应用中，可以通过数值微分或解析计算法向量的梯度
        // 暂时设为0（如果地形变化不大，这个近似是合理的）
    }
}

