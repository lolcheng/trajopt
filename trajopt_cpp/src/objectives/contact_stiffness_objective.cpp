#include "objectives/contact_stiffness_objective.hpp"
#include "polynomial.hpp"
#include <cmath>

ContactStiffnessObjective::ContactStiffnessObjective(
    int n_coeff, int n_sample,
    std::shared_ptr<RBFTerrain> terrain,
    double k_contact, double weight)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      k_contact_(k_contact),
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

void ContactStiffnessObjective::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

double ContactStiffnessObjective::smoothPos(double x) const {
    // smooth_pos(x) = max(0, x) 的平滑版本
    // 使用: x > 0 ? x : 0 (简化，实际可以用更平滑的函数)
    return x > 0.0 ? x : 0.0;
}

double ContactStiffnessObjective::computePenetration(double z_wheel, double h_terrain) const {
    // pen = smooth_pos(-(z_wheel - h_terrain)) = smooth_pos(h_terrain - z_wheel)
    return smoothPos(h_terrain - z_wheel);
}

double ContactStiffnessObjective::computeNormalForce(
    double fx, double fy, double fz,
    double nx, double ny, double nz) const {
    return fx * nx + fy * ny + fz * nz;
}

double ContactStiffnessObjective::eval_f(const double* x) {
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
    
    // 计算目标值: sum((fLn - K_CONTACT × pen_l)²)
    double obj = 0.0;
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形高度和法向量
        double h_l = terrain_->height(Xl[i], Yl[i]);
        double h_r = terrain_->height(Xr[i], Yr[i]);
        double nx_l, ny_l, nz_l;
        double nx_r, ny_r, nz_r;
        terrain_->normal(Xl[i], Yl[i], nx_l, ny_l, nz_l);
        terrain_->normal(Xr[i], Yr[i], nx_r, ny_r, nz_r);
        
        // 计算穿透量
        double pen_l = computePenetration(Zl[i], h_l);
        double pen_r = computePenetration(Zr[i], h_r);
        
        // 计算法向力
        double fLn = computeNormalForce(fLx[i], fLy[i], fLz[i], nx_l, ny_l, nz_l);
        double fRn = computeNormalForce(fRx[i], fRy[i], fRz[i], nx_r, ny_r, nz_r);
        
        // 弹簧模型预测的法向力
        double fLn_spring = k_contact_ * pen_l;
        double fRn_spring = k_contact_ * pen_r;
        
        double diff_l = fLn - fLn_spring;
        double diff_r = fRn - fRn_spring;
        obj += diff_l * diff_l + diff_r * diff_r;
    }
    
    return weight_ * obj;
}

void ContactStiffnessObjective::eval_grad_f(const double* x, double* grad_f) {
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
    double* grad_clx = grad_f + idx_clx_;
    double* grad_cly = grad_f + idx_cly_;
    double* grad_clz = grad_f + idx_clz_;
    double* grad_crx = grad_f + idx_crx_;
    double* grad_cry = grad_f + idx_cry_;
    double* grad_crz = grad_f + idx_crz_;
    
    // 计算梯度
    // d/dcoeff (fLn - K_CONTACT × pen_l)² = 2 * (fLn - K_CONTACT × pen_l) * (dfLn/dcoeff - K_CONTACT × dpen_l/dcoeff)
    for (int i = 0; i < n_sample_; ++i) {
        // 计算地形高度和法向量
        double h_l = terrain_->height(Xl[i], Yl[i]);
        double h_r = terrain_->height(Xr[i], Yr[i]);
        double nx_l, ny_l, nz_l;
        double nx_r, ny_r, nz_r;
        terrain_->normal(Xl[i], Yl[i], nx_l, ny_l, nz_l);
        terrain_->normal(Xr[i], Yr[i], nx_r, ny_r, nz_r);
        
        // 计算穿透量
        double pen_l = computePenetration(Zl[i], h_l);
        double pen_r = computePenetration(Zr[i], h_r);
        
        // 计算法向力
        double fLn = computeNormalForce(fLx[i], fLy[i], fLz[i], nx_l, ny_l, nz_l);
        double fRn = computeNormalForce(fRx[i], fRy[i], fRz[i], nx_r, ny_r, nz_r);
        
        // 弹簧模型预测的法向力
        double fLn_spring = k_contact_ * pen_l;
        double fRn_spring = k_contact_ * pen_r;
        
        double diff_l = fLn - fLn_spring;
        double diff_r = fRn - fRn_spring;
        double factor_l = 2.0 * weight_ * diff_l;
        double factor_r = 2.0 * weight_ * diff_r;
        
        // 对接触力的梯度
        grad_fLx[i] += factor_l * nx_l;
        grad_fLy[i] += factor_l * ny_l;
        grad_fLz[i] += factor_l * nz_l;
        grad_fRx[i] += factor_r * nx_r;
        grad_fRy[i] += factor_r * ny_r;
        grad_fRz[i] += factor_r * nz_r;
        
        // 对轮子位置的梯度（通过穿透量和法向量）
        // dpen_l/dZl = -1 (如果pen_l > 0), 0 (否则)
        double dpen_l_dzl = (pen_l > 0.0) ? -1.0 : 0.0;
        double dpen_r_dzr = (pen_r > 0.0) ? -1.0 : 0.0;
        
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        for (int j = 0; j < n_coeff_; ++j) {
            grad_clz[j] += factor_l * (-k_contact_ * dpen_l_dzl) * Phi_row[j];
            grad_crz[j] += factor_r * (-k_contact_ * dpen_r_dzr) * Phi_row[j];
        }
    }
}
