#ifndef PATH_LENGTH_OBJECTIVE_HPP
#define PATH_LENGTH_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 轨迹长度惩罚项
 * 惩罚轨迹的总长度，让轨迹尽可能短
 * 
 * 使用数值积分计算轨迹长度：
 * L = ∫ sqrt((dx/ds)² + (dy/ds)² + (dz/ds)²) ds
 * 
 * 对于多项式轨迹，在采样点上计算导数，然后使用梯形法则积分
 */
class PathLengthObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param weight 权重（默认值会根据n_sample自动设置）
     */
    PathLengthObjective(int n_coeff, int n_sample, double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double weight_;
    bool enabled_;
    
    // 变量索引
    int idx_cbx_, idx_cby_, idx_cbz_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 预计算的导数基矩阵（用于计算 dx/ds, dy/ds, dz/ds）
    std::vector<double> Phi_dot_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：评估轨迹导数
    void evalTrajectoryDerivative(const double* coeff, double* traj_dot) const;
};

#endif // PATH_LENGTH_OBJECTIVE_HPP

