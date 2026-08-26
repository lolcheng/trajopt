#ifndef BOUNDARY_CONSTRAINT_8VAR_HPP
#define BOUNDARY_CONSTRAINT_8VAR_HPP

#include "constraint_base.hpp"
#include "rbf_terrain.hpp"
#include "basis_function.hpp"
#include <memory>

class BoundaryConstraint8Var : public ConstraintBase {
public:
    BoundaryConstraint8Var(int n_coeff,
                           double x0, double y0, double psi0,
                           double x1, double y1, double psi1,
                           double clear_z,
                           const RBFTerrain& terrain,
                           std::shared_ptr<BasisFunction> basis = nullptr,
                           double wheels_nominal = 0.5,
                           double base_xy_tol = 0.0,
                           double base_z_tol = 0.0,
                           double base_yaw_tol = 0.0,
                           double wheel_xy_tol = 0.0);
    ~BoundaryConstraint8Var() override;

    int getNumConstraints() const override { return 16; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* jac_g, int offset = 0) override;
    void getBounds(double* g_lb, double* g_ub, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    double x0_, y0_, psi0_;
    double x1_, y1_, psi1_;
    double clear_z_;
    double wheels_nominal_;
    double base_xy_tol_;
    double base_z_tol_;
    double base_yaw_tol_;
    double wheel_xy_tol_;
    const RBFTerrain& terrain_;
    std::shared_ptr<BasisFunction> basis_;
    double* Phi0_;
    double* Phi1_;
    int idx_cbx_, idx_cby_, idx_cbz_, idx_cpsi_;
    int idx_clx_, idx_cly_, idx_crx_, idx_cry_;
};

#endif
