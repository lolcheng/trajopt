#ifndef YAW_PATH_ALIGNMENT_OBJECTIVE_HPP
#define YAW_PATH_ALIGNMENT_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "polynomial.hpp"
#include <vector>

/**
 * 偏航-路径对齐目标项
 * 惩罚偏航角与路径切线方向不对齐: (1 - cos(alignment))²
 * 权重: W_YAW_ALIGN / n_sample = 50.0 / n_sample
 */
class YawPathAlignmentObjective : public ObjectiveBase {
public:
    /**
     * 构造函数
     * @param n_coeff 多项式系数数量
     * @param n_sample 采样点数量
     * @param weight 权重（默认 50.0 / n_sample）
     */
    YawPathAlignmentObjective(int n_coeff, int n_sample, double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    
    double eval_f(const double* x) override;
    
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    double weight_;
    bool enabled_;
    
    // 变量索引
    int idx_cbx_, idx_cby_, idx_cpsi_;
    
    // 预计算的采样点和基矩阵
    std::vector<double> s_samples_;
    std::vector<double> Phi_;  // (n_sample × n_coeff)
    
    // 辅助方法：评估轨迹
    void evalTrajectory(const double* coeff, double* traj) const;
    
    // 辅助方法：计算路径切线方向
    void computePathTangent(const std::vector<double>& Xb, const std::vector<double>& Yb,
                           std::vector<double>& tx, std::vector<double>& ty) const;
};

#endif // YAW_PATH_ALIGNMENT_OBJECTIVE_HPP

