#ifndef BASE_HEIGHT_STABILITY_OBJECTIVE_HPP
#define BASE_HEIGHT_STABILITY_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include <memory>
#include <vector>

class BaseHeightStabilityObjective : public ObjectiveBase {
public:
    BaseHeightStabilityObjective(int n_coeff, int n_sample,
                                 std::shared_ptr<RBFTerrain> terrain,
                                 double clear_z_nominal, double weight);
    bool isEnabled() const override { return enabled_; }
    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_, n_sample_;
    double clear_z_nominal_, weight_;
    bool enabled_;
    std::shared_ptr<RBFTerrain> terrain_;
    int idx_cbx_, idx_cby_, idx_cbz_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
};

#endif
