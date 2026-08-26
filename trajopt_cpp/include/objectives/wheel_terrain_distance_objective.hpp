#ifndef WHEEL_TERRAIN_DISTANCE_OBJECTIVE_HPP
#define WHEEL_TERRAIN_DISTANCE_OBJECTIVE_HPP

#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include "polynomial.hpp"
#include <vector>

class WheelTerrainDistanceObjective : public ObjectiveBase {
public:
    WheelTerrainDistanceObjective(int n_coeff, int n_sample,
                                  const RBFTerrain& terrain,
                                  double weight = -1.0);
    
    bool isEnabled() const override { return enabled_; }
    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;

private:
    int n_coeff_;
    int n_sample_;
    const RBFTerrain& terrain_;
    double weight_;
    bool enabled_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    void evalTrajectory(const double* coeff, double* traj) const;
    void evalTrajectoryGrad(const double* traj_diff, const double* Phi_row,
                           double* grad_coeff) const;
};

#endif
