#ifndef SEG_TRAJOPT_NLP_HPP
#define SEG_TRAJOPT_NLP_HPP

#include <stddef.h>
#include <vector>
#include <memory>
#include "IpTNLP.hpp"
#include "IpSmartPtr.hpp"
#include "constraint_base.hpp"
#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include "basis_function.hpp"

/** 分段轨迹优化 NLP：决策变量仅为切比雪夫多项式系数（8 组，无接触力）。 */
class SegTrajOptNLP : public Ipopt::TNLP {
public:
    SegTrajOptNLP(int n_coeff, int n_sample,
                  std::shared_ptr<BasisFunction> basis = nullptr);

    virtual ~SegTrajOptNLP();

    void addConstraint(std::shared_ptr<ConstraintBase> constraint);
    void addObjective(std::shared_ptr<ObjectiveBase> objective);
    void setTerrain(std::shared_ptr<RBFTerrain> terrain);
    std::shared_ptr<RBFTerrain> getTerrain() const { return terrain_; }

    void extractTrajectory(const Ipopt::Number* x,
                           std::vector<double>& Xb, std::vector<double>& Yb,
                           std::vector<double>& Zb, std::vector<double>& Psib,
                           std::vector<double>& Xl, std::vector<double>& Yl, std::vector<double>& Zl,
                           std::vector<double>& Xr, std::vector<double>& Yr, std::vector<double>& Zr) const;

    void extractCoefficients(const Ipopt::Number* x,
                              std::vector<double>& cbx, std::vector<double>& cby,
                              std::vector<double>& cbz, std::vector<double>& cpsi,
                              std::vector<double>& clx, std::vector<double>& cly, std::vector<double>& clz,
                              std::vector<double>& crx, std::vector<double>& cry, std::vector<double>& crz) const;

    bool hasSolution() const { return solution_valid_; }
    const std::vector<Ipopt::Number>& getSolution() const { return solution_; }
    Ipopt::SolverReturn getSolverStatus() const { return solver_status_; }
    int getNumVars() const { return n_vars_; }
    int getNumConstraints() const { return n_constraints_; }
    Ipopt::Index get_n_constraints() const { return static_cast<Ipopt::Index>(n_constraints_); }

    void setInitialGuess(const Ipopt::Number* x);
    const Ipopt::Number* getInitialGuess() const { return initial_guess_.empty() ? nullptr : initial_guess_.data(); }

    int get_n_coeff() const { return n_coeff_; }
    int get_n_sample() const { return n_sample_; }
    std::shared_ptr<BasisFunction> getBasis() const { return basis_; }

    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                      Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style) override;
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;
    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;
    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Index m, Ipopt::Number* g) override;
    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                    Ipopt::Index* jCol, Ipopt::Number* values) override;
    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                Ipopt::Index* jCol, Ipopt::Number* values) override;
    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                          const Ipopt::Number* x, const Ipopt::Number* z_L,
                          const Ipopt::Number* z_U, Ipopt::Index m,
                          const Ipopt::Number* g, const Ipopt::Number* lambda,
                          Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
                          Ipopt::IpoptCalculatedQuantities* ip_cq) override;

private:
    int n_coeff_;
    int n_sample_;
    int n_vars_;
    int n_constraints_;

    std::vector<std::shared_ptr<ConstraintBase>> constraints_;
    std::vector<std::shared_ptr<ObjectiveBase>> objectives_;
    std::shared_ptr<RBFTerrain> terrain_;
    std::shared_ptr<BasisFunction> basis_;

    int idx_cbx_, idx_cby_, idx_cbz_, idx_cpsi_;
    int idx_clx_, idx_cly_;
    int idx_crx_, idx_cry_;

    void updateDimensions();
    void computeJacobianStructure();
    std::vector<Ipopt::Index> jac_iRow_;
    std::vector<Ipopt::Index> jac_jCol_;
    bool jac_structure_computed_;

    std::vector<Ipopt::Number> solution_;
    bool solution_valid_;
    Ipopt::SolverReturn solver_status_;
    std::vector<Ipopt::Number> initial_guess_;
};

#endif
