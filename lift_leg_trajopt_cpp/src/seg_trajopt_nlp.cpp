#include "seg_trajopt_nlp.hpp"
#include "polynomial.hpp"
#include "chebyshev_basis.hpp"
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>

SegTrajOptNLP::SegTrajOptNLP(int n_coeff, int n_sample,
                             std::shared_ptr<BasisFunction> basis,
                             int n_seg_z)
    : n_coeff_(n_coeff), n_seg_z_(std::max(1, n_seg_z)),
      n_coeff_z_(n_coeff_ * n_seg_z_), n_sample_(n_sample),
      n_vars_(0), n_constraints_(0),
      jac_structure_computed_(false),
      solution_valid_(false),
      solver_status_(static_cast<Ipopt::SolverReturn>(0)) {

    if (basis == nullptr)
        basis_ = std::make_shared<ChebyshevBasis>();
    else
        basis_ = basis;

    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = idx_clz_ + n_coeff_z_;
    idx_cry_ = idx_crx_ + n_coeff_;
    idx_crz_ = idx_cry_ + n_coeff_;

    n_vars_ = idx_crz_ + n_coeff_z_;
}
void SegTrajOptNLP::buildPiecewiseZBasis(std::vector<double>& Phi_z) const {
    Phi_z.assign(n_sample_ * n_coeff_z_, 0.0);
    for (int i = 0; i < n_sample_; ++i) {
        const double s = (n_sample_ > 1) ? static_cast<double>(i) / (n_sample_ - 1) : 0.0;
        int seg = std::min(static_cast<int>(std::floor(s * n_seg_z_)), n_seg_z_ - 1);
        const double local = std::clamp(s * n_seg_z_ - seg, 0.0, 1.0);
        std::vector<double> phi_local(n_coeff_);
        basis_->computeBasisMatrix(&local, 1, n_coeff_, phi_local.data());
        const int row_off = i * n_coeff_z_;
        const int col_off = seg * n_coeff_;
        for (int k = 0; k < n_coeff_; ++k) {
            Phi_z[row_off + col_off + k] = phi_local[k];
        }
    }
}


SegTrajOptNLP::~SegTrajOptNLP() = default;

void SegTrajOptNLP::addConstraint(std::shared_ptr<ConstraintBase> constraint) {
    constraints_.push_back(constraint);
    updateDimensions();
    jac_structure_computed_ = false;
}

void SegTrajOptNLP::addObjective(std::shared_ptr<ObjectiveBase> objective) {
    objectives_.push_back(objective);
}

void SegTrajOptNLP::setTerrain(std::shared_ptr<RBFTerrain> terrain) {
    terrain_ = terrain;
}

void SegTrajOptNLP::updateDimensions() {
    n_constraints_ = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled())
            n_constraints_ += c->getNumConstraints();
    }
}

bool SegTrajOptNLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                                 Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                                 Ipopt::TNLP::IndexStyleEnum& index_style) {
    n = n_vars_;
    m = n_constraints_;
    if (!jac_structure_computed_)
        computeJacobianStructure();
    nnz_jac_g = static_cast<Ipopt::Index>(jac_iRow_.size());
    nnz_h_lag = 0;
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
}

bool SegTrajOptNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                    Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) {
    for (Ipopt::Index i = 0; i < n; ++i) {
        x_l[i] = -1e20;
        x_u[i] = 1e20;
    }
    int offset = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            c->getBounds(g_l + offset, g_u + offset, offset);
            offset += c->getNumConstraints();
        }
    }
    return true;
}

bool SegTrajOptNLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                       bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                       Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) {
    (void)z_L;
    (void)z_U;
    (void)init_lambda;
    (void)lambda;
    if (init_x) {
        if (!initial_guess_.empty() && initial_guess_.size() == static_cast<size_t>(n))
            std::copy(initial_guess_.begin(), initial_guess_.end(), x);
        else {
            for (Ipopt::Index i = 0; i < n; ++i)
                x[i] = 0.0;
        }
    }
    return true;
}

void SegTrajOptNLP::setInitialGuess(const Ipopt::Number* x) {
    if (x != nullptr) {
        initial_guess_.resize(n_vars_);
        std::copy(x, x + n_vars_, initial_guess_.begin());
    } else
        initial_guess_.clear();
}

bool SegTrajOptNLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                           Ipopt::Number& obj_value) {
    (void)new_x;
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (!std::isfinite(x[i])) {
            std::cerr << "SegTrajOptNLP: invalid x[" << i << "] = " << x[i] << std::endl;
            return false;
        }
    }
    obj_value = 0.0;
    for (const auto& obj : objectives_) {
        if (obj->isEnabled()) {
            double v = obj->eval_f(x);
            if (!std::isfinite(v)) {
                std::cerr << "SegTrajOptNLP: invalid objective value" << std::endl;
                return false;
            }
            obj_value += v;
        }
    }
    return true;
}

bool SegTrajOptNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                                Ipopt::Number* grad_f) {
    (void)new_x;
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (!std::isfinite(x[i])) return false;
    }
    for (Ipopt::Index i = 0; i < n; ++i)
        grad_f[i] = 0.0;
    for (const auto& obj : objectives_) {
        if (obj->isEnabled())
            obj->eval_grad_f(x, grad_f);
    }
    return true;
}

bool SegTrajOptNLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                           Ipopt::Index m, Ipopt::Number* g) {
    (void)new_x;
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (!std::isfinite(x[i])) return false;
    }
    int offset = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            c->eval_g(x, g + offset);
            offset += c->getNumConstraints();
        }
    }
    return true;
}

bool SegTrajOptNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                               Ipopt::Index m, Ipopt::Index nele_jac,
                               Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) {
    (void)n;
    (void)x;
    (void)new_x;
    (void)m;
    if (values == nullptr) {
        if (!jac_structure_computed_)
            computeJacobianStructure();
        for (Ipopt::Index i = 0; i < nele_jac; ++i) {
            iRow[i] = jac_iRow_[i];
            jCol[i] = jac_jCol_[i];
        }
    } else {
        int constraint_offset = 0;
        int jac_offset = 0;
        for (const auto& c : constraints_) {
            if (c->isEnabled()) {
                int nnz = c->getJacobianStructure(nullptr, nullptr, 0);
                c->eval_jac_g(x, values + jac_offset, constraint_offset);
                jac_offset += nnz;
                constraint_offset += c->getNumConstraints();
            }
        }
    }
    return true;
}

bool SegTrajOptNLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                           Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                           bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                           Ipopt::Index* jCol, Ipopt::Number* values) {
    (void)n;
    (void)x;
    (void)new_x;
    (void)obj_factor;
    (void)m;
    (void)lambda;
    (void)new_lambda;
    (void)nele_hess;
    (void)iRow;
    (void)jCol;
    (void)values;
    return true;
}

void SegTrajOptNLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                                      const Ipopt::Number* x, const Ipopt::Number* z_L,
                                      const Ipopt::Number* z_U, Ipopt::Index m,
                                      const Ipopt::Number* g, const Ipopt::Number* lambda,
                                      Ipopt::Number obj_value,
                                      const Ipopt::IpoptData* ip_data,
                                      Ipopt::IpoptCalculatedQuantities* ip_cq) {
    (void)z_L;
    (void)z_U;
    (void)m;
    (void)g;
    (void)lambda;
    (void)ip_data;
    (void)ip_cq;
    solver_status_ = status;
    std::cout << "\n========== SegTrajOpt Solution ==========\n";
    std::cout << "Status: " << status << "\n";
    std::cout << "Objective: " << obj_value << "\n";
    solution_.resize(n);
    std::memcpy(solution_.data(), x, n * sizeof(Ipopt::Number));
    solution_valid_ = true;
}

void SegTrajOptNLP::computeJacobianStructure() {
    jac_iRow_.clear();
    jac_jCol_.clear();
    int constraint_offset = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            int nnz = c->getJacobianStructure(nullptr, nullptr, 0);
            std::vector<int> iRow(nnz), jCol(nnz);
            int actual = c->getJacobianStructure(iRow.data(), jCol.data(), constraint_offset);
            for (int i = 0; i < actual; ++i) {
                jac_iRow_.push_back(iRow[i]);
                jac_jCol_.push_back(jCol[i]);
            }
            constraint_offset += c->getNumConstraints();
        }
    }
    jac_structure_computed_ = true;
}

void SegTrajOptNLP::extractTrajectory(const Ipopt::Number* x_in,
                                      std::vector<double>& Xb, std::vector<double>& Yb,
                                      std::vector<double>& Zb, std::vector<double>& Psib,
                                      std::vector<double>& Xl, std::vector<double>& Yl, std::vector<double>& Zl,
                                      std::vector<double>& Xr, std::vector<double>& Yr, std::vector<double>& Zr) const {
    const Ipopt::Number* x = x_in;
    if (x == nullptr) {
        if (!solution_valid_ || solution_.empty()) {
            std::cerr << "SegTrajOptNLP: no solution for trajectory extraction" << std::endl;
            return;
        }
        x = solution_.data();
    }
    std::vector<double> s(n_sample_);
    for (int i = 0; i < n_sample_; ++i)
        s[i] = static_cast<double>(i) / (n_sample_ - 1);
    std::vector<double> Phi(n_sample_ * n_coeff_);
    basis_->computeBasisMatrix(s.data(), n_sample_, n_coeff_, Phi.data());

    const double* cbx = x + idx_cbx_, * cby = x + idx_cby_, * cbz = x + idx_cbz_, * cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_, * cly = x + idx_cly_, * clz = x + idx_clz_;
    const double* crx = x + idx_crx_, * cry = x + idx_cry_, * crz = x + idx_crz_;

    Xb.resize(n_sample_); Yb.resize(n_sample_); Zb.resize(n_sample_); Psib.resize(n_sample_);
    Xl.resize(n_sample_); Yl.resize(n_sample_); Zl.resize(n_sample_);
    Xr.resize(n_sample_); Yr.resize(n_sample_); Zr.resize(n_sample_);

    basis_->eval(Phi.data(), n_sample_, cbx, n_coeff_, Xb.data());
    basis_->eval(Phi.data(), n_sample_, cby, n_coeff_, Yb.data());
    basis_->eval(Phi.data(), n_sample_, cbz, n_coeff_, Zb.data());
    basis_->eval(Phi.data(), n_sample_, cpsi, n_coeff_, Psib.data());
    basis_->eval(Phi.data(), n_sample_, clx, n_coeff_, Xl.data());
    basis_->eval(Phi.data(), n_sample_, cly, n_coeff_, Yl.data());
    std::vector<double> Phi_z;
    buildPiecewiseZBasis(Phi_z);
    basis_->eval(Phi_z.data(), n_sample_, clz, n_coeff_z_, Zl.data());
    basis_->eval(Phi.data(), n_sample_, crx, n_coeff_, Xr.data());
    basis_->eval(Phi.data(), n_sample_, cry, n_coeff_, Yr.data());
    basis_->eval(Phi_z.data(), n_sample_, crz, n_coeff_z_, Zr.data());
}

void SegTrajOptNLP::extractCoefficients(const Ipopt::Number* x_in,
                                        std::vector<double>& cbx, std::vector<double>& cby,
                                        std::vector<double>& cbz, std::vector<double>& cpsi,
                                        std::vector<double>& clx, std::vector<double>& cly, std::vector<double>& clz,
                                        std::vector<double>& crx, std::vector<double>& cry, std::vector<double>& crz) const {
    const Ipopt::Number* x = x_in;
    if (x == nullptr) {
        if (!solution_valid_ || solution_.empty()) return;
        x = solution_.data();
    }
    cbx.resize(n_coeff_); cby.resize(n_coeff_); cbz.resize(n_coeff_); cpsi.resize(n_coeff_);
    clx.resize(n_coeff_); cly.resize(n_coeff_); clz.resize(n_coeff_z_);
    crx.resize(n_coeff_); cry.resize(n_coeff_); crz.resize(n_coeff_z_);
    std::copy(x + idx_cbx_, x + idx_cbx_ + n_coeff_, cbx.begin());
    std::copy(x + idx_cby_, x + idx_cby_ + n_coeff_, cby.begin());
    std::copy(x + idx_cbz_, x + idx_cbz_ + n_coeff_, cbz.begin());
    std::copy(x + idx_cpsi_, x + idx_cpsi_ + n_coeff_, cpsi.begin());
    std::copy(x + idx_clx_, x + idx_clx_ + n_coeff_, clx.begin());
    std::copy(x + idx_cly_, x + idx_cly_ + n_coeff_, cly.begin());
    std::copy(x + idx_clz_, x + idx_clz_ + n_coeff_z_, clz.begin());
    std::copy(x + idx_crx_, x + idx_crx_ + n_coeff_, crx.begin());
    std::copy(x + idx_cry_, x + idx_cry_ + n_coeff_, cry.begin());
    std::copy(x + idx_crz_, x + idx_crz_ + n_coeff_z_, crz.begin());
}
