#ifndef PENETRATION_BASED_FORCE_CONSTRAINT_HPP
#define PENETRATION_BASED_FORCE_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "rbf_terrain.hpp"
#include <memory>
#include <vector>

class PenetrationBasedForceConstraint : public ConstraintBase {
public:
    PenetrationBasedForceConstraint(int n_coeff, int n_sample,
                                    int n_seg_z,
                                    std::shared_ptr<RBFTerrain> terrain,
                                    double stiffness,
                                    const std::vector<int>& left_contact_mask,
                                    const std::vector<int>& right_contact_mask,
                                    double stance_force_min,
                                    double stance_force_max,
                                    double swing_force_max,
                                    double smooth_eps = 1e-6);
    ~PenetrationBasedForceConstraint() override = default;

    int getNumConstraints() const override { return 2 * n_sample_; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_l, double* g_u, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    double smoothPositive(double d) const;
    double smoothPositiveGrad(double d) const;
    double estimateForce(double penetration) const;

    int n_coeff_, n_seg_z_, n_coeff_z_, n_sample_;
    double stiffness_;
    double stance_force_min_, stance_force_max_, swing_force_max_;
    double smooth_eps_;
    std::shared_ptr<RBFTerrain> terrain_;
    std::vector<int> left_contact_mask_;
    std::vector<int> right_contact_mask_;
    int idx_clx_, idx_cly_, idx_clz_, idx_crx_, idx_cry_, idx_crz_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    std::vector<double> Phi_z_;
};

#endif
