#ifndef CONTACT_STIFFNESS_OBJECTIVE_HPP
#define CONTACT_STIFFNESS_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include "polynomial.hpp"
#include <vector>
#include <memory>

/**
 * 接触刚度目标项
 * 惩罚接触力与穿透量弹簧模型不符: (fLn - K_CONTACT × pen_l)²
 * 参数: K_CONTACT = 500.0
 * 权重: W_CONTACT_STIFF / n_sample = 1000.0 / n_sample
 */
class ContactStiffnessObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param terrain RBF地形
     * @param k_contact 接触刚度（默认 500.0）
     * @param weight 权重（默认 1000.0 / n_sample）
     */
    ContactStiffnessObjective(int n_coeff, int n_sample,
                             std::shared_ptr<RBFTerrain> terrain,
                             double k_contact = 500.0,
                             double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double k_contact_;
    double weight_;
    bool enabled_;
    std::shared_ptr<RBFTerrain> terrain_;
    
    // 变量索引
    int idx_fLx_, idx_fLy_, idx_fLz_;
    int idx_fRx_, idx_fRy_, idx_fRz_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：平滑正函数（smooth_pos）
    double smoothPos(double x) const;
    
    // 辅助方法：计算穿透量
    double computePenetration(double z_wheel, double h_terrain) const;
    
    // 辅助方法：计算法向力
    double computeNormalForce(double fx, double fy, double fz,
                              double nx, double ny, double nz) const;
};

#endif // CONTACT_STIFFNESS_OBJECTIVE_HPP

