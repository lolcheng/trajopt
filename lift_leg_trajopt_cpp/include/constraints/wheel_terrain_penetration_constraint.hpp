#ifndef WHEEL_TERRAIN_PENETRATION_CONSTRAINT_HPP
#define WHEEL_TERRAIN_PENETRATION_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "rbf_terrain.hpp"
#include <memory>
#include <vector>

class WheelTerrainPenetrationConstraint : public ConstraintBase {
public:
    WheelTerrainPenetrationConstraint(int n_coeff, int n_sample,
                                      int n_seg_z,
                                      std::shared_ptr<RBFTerrain> terrain,
                                      double max_penetration);
    ~WheelTerrainPenetrationConstraint() override = default;

    int getNumConstraints() const override { return 2 * n_sample_; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_, n_seg_z_, n_coeff_z_, n_sample_;
    double max_penetration_;
    std::shared_ptr<RBFTerrain> terrain_;
    int idx_clx_, idx_cly_, idx_clz_, idx_crx_, idx_cry_, idx_crz_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    std::vector<double> Phi_z_;
};

#endif
