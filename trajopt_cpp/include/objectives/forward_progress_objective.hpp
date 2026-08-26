#ifndef FORWARD_PROGRESS_OBJECTIVE_HPP
#define FORWARD_PROGRESS_OBJECTIVE_HPP

#include "objective_base.hpp"
#include <vector>

/**
 * 前向进度惩罚项（软约束）
 * 仅惩罚沿起点->终点方向的“回头”段，减少蛇形和反复转弯。
 */
class ForwardProgressObjective : public ObjectiveBase {
public:
    ForwardProgressObjective(int n_coeff, int n_sample,
                             double xb0, double yb0,
                             double xb1, double yb1,
                             double weight = -1.0);

    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;
    bool isEnabled() const override { return enabled_; }

private:
    int n_coeff_;
    int n_sample_;
    double weight_;
    bool enabled_;

    // 目标方向单位向量（起点 -> 终点）
    double dir_x_;
    double dir_y_;

    int idx_cbx_;
    int idx_cby_;

    std::vector<double> s_samples_;
    std::vector<double> Phi_;

    void evalTrajectory(const double* coeff, double* traj) const;

    static double smoothPos(double z, double eps = 1e-6);
    static double smoothPosGrad(double z, double eps = 1e-6);
};

#endif // FORWARD_PROGRESS_OBJECTIVE_HPP


