#ifndef WHEELS_DISTANCE_TERRAIN_CONSTRAINT_HPP
#define WHEELS_DISTANCE_TERRAIN_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include <vector>

class WheelsDistanceTerrainConstraint : public ConstraintBase {
public:
    WheelsDistanceTerrainConstraint(int n_coeff, int n_sample,
                                    int n_seg_z,
                                    double min_dist, double max_dist);
    ~WheelsDistanceTerrainConstraint() override = default;

    int getNumConstraints() const override { return n_sample_; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_, n_seg_z_, n_coeff_z_, n_sample_;
    double min_dist_sq_, max_dist_sq_;
    int idx_clx_, idx_cly_, idx_clz_, idx_crx_, idx_cry_, idx_crz_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    std::vector<double> Phi_z_;
};

#endif
