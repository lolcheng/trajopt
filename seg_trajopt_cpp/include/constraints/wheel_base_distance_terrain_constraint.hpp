#ifndef WHEEL_BASE_DISTANCE_TERRAIN_CONSTRAINT_HPP
#define WHEEL_BASE_DISTANCE_TERRAIN_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "rbf_terrain.hpp"
#include <memory>
#include <vector>

class WheelBaseDistanceTerrainConstraint : public ConstraintBase {
public:
    WheelBaseDistanceTerrainConstraint(int n_coeff, int n_sample,
                                       std::shared_ptr<RBFTerrain> terrain,
                                       double min_dist, double max_dist);
    ~WheelBaseDistanceTerrainConstraint() override = default;

    int getNumConstraints() const override { return 2 * n_sample_; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_, n_sample_;
    double min_dist_sq_, max_dist_sq_;
    std::shared_ptr<RBFTerrain> terrain_;
    int idx_cbx_, idx_cby_, idx_cbz_, idx_clx_, idx_cly_, idx_crx_, idx_cry_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
};

#endif
