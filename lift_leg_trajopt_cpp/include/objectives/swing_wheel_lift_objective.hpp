#ifndef SWING_WHEEL_LIFT_OBJECTIVE_HPP
#define SWING_WHEEL_LIFT_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include <memory>
#include <vector>

class SwingWheelLiftObjective : public ObjectiveBase {
public:
    SwingWheelLiftObjective(int n_coeff, int n_sample,
                            int n_seg_z,
                            std::shared_ptr<RBFTerrain> terrain,
                            bool swing_left,
                            double lift_peak,
                            double stance_clearance_target,
                            double weight);

    bool isEnabled() const override { return enabled_; }
    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    double profile(double s) const;

    int n_coeff_, n_seg_z_, n_coeff_z_, n_sample_;
    bool swing_left_;
    double lift_peak_;
    double stance_clearance_target_;
    double weight_;
    bool enabled_;
    std::shared_ptr<RBFTerrain> terrain_;
    int idx_clx_, idx_cly_, idx_clz_, idx_crx_, idx_cry_, idx_crz_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    std::vector<double> Phi_z_;
};

#endif
