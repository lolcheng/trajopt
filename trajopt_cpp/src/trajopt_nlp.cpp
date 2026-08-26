#include "trajopt_nlp.hpp"
#include "polynomial.hpp"
#include "polynomial_basis.hpp"
#include "chebyshev_basis.hpp"
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cerrno>

TrajOptNLP::TrajOptNLP(int n_coeff, int n_sample, 
                       std::shared_ptr<BasisFunction> basis)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      n_vars_(0), n_constraints_(0),
      jac_structure_computed_(false),
      solution_valid_(false),
      solver_status_(static_cast<Ipopt::SolverReturn>(0)) {  // 初始化为0，会在finalize_solution中更新
    
    // 如果未提供基函数，默认使用切比雪夫基
    if (basis == nullptr) {
        basis_ = std::make_shared<ChebyshevBasis>();
    } else {
        basis_ = basis;
    }
    
    // 计算变量索引
    idx_cbx_ = 0;
    idx_cby_ = n_coeff_;
    idx_cbz_ = 2 * n_coeff_;
    idx_cpsi_ = 3 * n_coeff_;
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = 7 * n_coeff_;
    idx_cry_ = 8 * n_coeff_;
    idx_crz_ = 9 * n_coeff_;
    
    // 接触力变量索引
    int n_poly_vars = 10 * n_coeff_;  // 10个轨迹 × n_coeff
    idx_fLx_ = n_poly_vars;
    idx_fLy_ = idx_fLx_ + n_sample;
    idx_fLz_ = idx_fLy_ + n_sample;
    idx_fRx_ = idx_fLz_ + n_sample;
    idx_fRy_ = idx_fRx_ + n_sample;
    idx_fRz_ = idx_fRy_ + n_sample;
    
    // 总变量数：多项式系数 + 接触力
    n_vars_ = n_poly_vars + 6 * n_sample;  // 10个轨迹 × n_coeff + 6个力 × n_sample
}

TrajOptNLP::~TrajOptNLP() = default;

void TrajOptNLP::addConstraint(std::shared_ptr<ConstraintBase> constraint) {
    constraints_.push_back(constraint);
    updateDimensions();
    jac_structure_computed_ = false;
}

void TrajOptNLP::addObjective(std::shared_ptr<ObjectiveBase> objective) {
    objectives_.push_back(objective);
}

void TrajOptNLP::setTerrain(std::shared_ptr<RBFTerrain> terrain) {
    terrain_ = terrain;
}

void TrajOptNLP::updateDimensions() {
    n_constraints_ = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            n_constraints_ += c->getNumConstraints();
        }
    }
}

bool TrajOptNLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                               Ipopt::Index& nnz_jac_g,
                               Ipopt::Index& nnz_h_lag,
                               Ipopt::TNLP::IndexStyleEnum& index_style) {
    n = n_vars_;
    m = n_constraints_;
    
    // 计算雅可比稀疏结构
    if (!jac_structure_computed_) {
        computeJacobianStructure();
    }
    nnz_jac_g = static_cast<Ipopt::Index>(jac_iRow_.size());
    
    // 海塞矩阵（当前版本设为0，使用有限差分）
    nnz_h_lag = 0;
    
    index_style = Ipopt::TNLP::C_STYLE;
    
    return true;
}

bool TrajOptNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                  Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) {
    // 变量边界
    // 多项式系数：无界
    for (Ipopt::Index i = 0; i < 10 * n_coeff_; ++i) {
        x_l[i] = -1e20;
        x_u[i] = 1e20;
    }
    // 接触力：无界（约束会限制法向力）
    for (Ipopt::Index i = 10 * n_coeff_; i < n; ++i) {
        x_l[i] = -1e20;
        x_u[i] = 1e20;
    }
    
    // 约束边界
    int offset = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            c->getBounds(g_l + offset, g_u + offset, offset);
            offset += c->getNumConstraints();
        }
    }
    
    return true;
}

bool TrajOptNLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) {
    if (init_x) {
        if (!initial_guess_.empty() && initial_guess_.size() == static_cast<size_t>(n)) {
            // 使用设置的初始猜测
            std::copy(initial_guess_.begin(), initial_guess_.end(), x);
        } else {
            // 默认初始猜测：设为0
            for (Ipopt::Index i = 0; i < n; ++i) {
                x[i] = 0.0;
            }
        }
    }
    
    return true;
}

void TrajOptNLP::setInitialGuess(const Ipopt::Number* x) {
    if (x != nullptr) {
        initial_guess_.resize(n_vars_);
        std::copy(x, x + n_vars_, initial_guess_.begin());
    } else {
        initial_guess_.clear();
    }
}

bool TrajOptNLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Number& obj_value) {
    // 检查输入变量是否有NaN/Inf（静默检查，只在出错时输出）
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (!std::isfinite(x[i])) {
            std::cerr << "\n*** ERROR: Invalid number detected in variable x[" << i << "] = " << x[i] << std::endl;
            // 输出变量类型信息
            if (i < 10 * n_coeff_) {
                int var_idx = i / n_coeff_;
                const char* var_names[] = {"cbx", "cby", "cbz", "cpsi", "clx", "cly", "clz", "crx", "cry", "crz"};
                if (var_idx < 10) {
                    std::cerr << "  Variable: " << var_names[var_idx] << "[" << (i % n_coeff_) << "]" << std::endl;
                }
            } else {
                std::cerr << "  This is a contact force variable" << std::endl;
            }
            return false;
        }
    }
    
    obj_value = 0.0;
    
    // 累加所有目标函数项
    int obj_idx = 0;
    for (const auto& obj : objectives_) {
        if (obj->isEnabled()) {
            double obj_contrib = obj->eval_f(x);
            if (!std::isfinite(obj_contrib)) {
                std::cerr << "\n*** ERROR: Invalid number in objective function #" << obj_idx << std::endl;
                std::cerr << "  Contribution value: " << obj_contrib << std::endl;
                // 尝试识别是哪个目标函数
                const char* obj_names[] = {
                    "WheelTerrainDistance", "Regularization", "Smoothing",
                    "WheelBaseDistanceNominal", "WheelsDistanceNominal",
                    "YawPathAlignment", "ContactStiffness", "ForceBalance",
                    "TorqueBalance", "FrictionCone"
                };
                if (obj_idx < 10) {
                    std::cerr << "  Objective: " << obj_names[obj_idx] << std::endl;
                }
                return false;
            }
            obj_value += obj_contrib;
            obj_idx++;
        }
    }
    
    if (!std::isfinite(obj_value)) {
        std::cerr << "\n*** ERROR: Invalid number in total objective value: " << obj_value << std::endl;
        return false;
    }
    
    return true;
}

bool TrajOptNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                              Ipopt::Number* grad_f) {
    // 检查输入变量（静默检查）
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (!std::isfinite(x[i])) {
            std::cerr << "\n*** ERROR: Invalid number in x[" << i << "] in eval_grad_f: " << x[i] << std::endl;
            return false;
        }
    }
    
    // 初始化梯度为0
    for (Ipopt::Index i = 0; i < n; ++i) {
        grad_f[i] = 0.0;
    }
    
    // 累加所有目标函数项的梯度
    int obj_idx = 0;
    for (const auto& obj : objectives_) {
        if (obj->isEnabled()) {
            obj->eval_grad_f(x, grad_f);
            // 检查梯度是否有NaN/Inf
            for (Ipopt::Index i = 0; i < n; ++i) {
                if (!std::isfinite(grad_f[i])) {
                    std::cerr << "\n*** ERROR: Invalid number in grad_f[" << i << "] after objective #" << obj_idx << std::endl;
                    std::cerr << "  Gradient value: " << grad_f[i] << std::endl;
                    const char* obj_names[] = {
                        "WheelTerrainDistance", "Regularization", "Smoothing",
                        "WheelBaseDistanceNominal", "WheelsDistanceNominal",
                        "YawPathAlignment", "ContactStiffness", "ForceBalance",
                        "TorqueBalance", "FrictionCone"
                    };
                    if (obj_idx < 10) {
                        std::cerr << "  Objective: " << obj_names[obj_idx] << std::endl;
                    }
                    return false;
                }
            }
            obj_idx++;
        }
    }
    
    return true;
}

bool TrajOptNLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Index m, Ipopt::Number* g) {
    // 检查输入变量是否有NaN/Inf（静默检查，只在出错时输出）
    for (Ipopt::Index i = 0; i < n; ++i) {
        if (!std::isfinite(x[i])) {
            std::cerr << "\n*** ERROR: Invalid number detected in variable x[" << i << "] = " << x[i] << std::endl;
            return false;
        }
    }
    
    int offset = 0;
    int constraint_idx = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            c->eval_g(x, g + offset);
            // 检查约束值是否有NaN/Inf
            int n_cons = c->getNumConstraints();
            for (int i = 0; i < n_cons; ++i) {
                if (!std::isfinite(g[offset + i])) {
                    std::cerr << "\n*** ERROR: Invalid number in constraint #" << constraint_idx 
                              << ", g[" << (offset + i) << "] = " << g[offset + i] << std::endl;
                    // 尝试识别是哪个约束
                    const char* constraint_names[] = {
                        "BoundaryConstraint", "WheelBaseDistanceConstraint",
                        "WheelsDistanceConstraint", "AntiCrossingConstraint",
                        "NormalForceUpperBoundConstraint"
                    };
                    if (constraint_idx < 5) {
                        std::cerr << "  Constraint: " << constraint_names[constraint_idx] << std::endl;
                    }
                    return false;
                }
            }
            offset += n_cons;
            constraint_idx++;
        }
    }
    return true;
}

bool TrajOptNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                            Ipopt::Index m, Ipopt::Index nele_jac,
                            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) {
    if (values == nullptr) {
        // 返回稀疏结构
        if (!jac_structure_computed_) {
            computeJacobianStructure();
        }
        for (Ipopt::Index i = 0; i < nele_jac; ++i) {
            iRow[i] = jac_iRow_[i];
            jCol[i] = jac_jCol_[i];
        }
    } else {
        // 返回雅可比值
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

bool TrajOptNLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Number obj_factor, Ipopt::Index m,
                        const Ipopt::Number* lambda, bool new_lambda,
                        Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values) {
    // 当前版本不使用海塞矩阵
    return true;
}

void TrajOptNLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                                   const Ipopt::Number* x, const Ipopt::Number* z_L,
                                   const Ipopt::Number* z_U, Ipopt::Index m,
                                   const Ipopt::Number* g, const Ipopt::Number* lambda,
                                   Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data,
                                   Ipopt::IpoptCalculatedQuantities* ip_cq) {
    // 保存求解器状态码（SolverReturn，与ApplicationReturnStatus可能不同）
    solver_status_ = status;
    
    std::cout << "\n========== Solution ==========\n";
    std::cout << "Status: " << status << "\n";
    std::cout << "Objective value: " << obj_value << "\n";
    std::cout << "Variables:\n";
    for (Ipopt::Index i = 0; i < n && i < 20; ++i) {
        std::cout << "  x[" << i << "] = " << x[i] << "\n";
    }
    if (n > 20) {
        std::cout << "  ... (showing first 20)\n";
    }
    
    // 保存解供后续使用（即使失败也保存，用于可视化）
    solution_.resize(n);
    std::memcpy(solution_.data(), x, n * sizeof(Ipopt::Number));
    // 总是标记为有效，即使失败也允许可视化
    solution_valid_ = true;
}

void TrajOptNLP::computeJacobianStructure() {
    jac_iRow_.clear();
    jac_jCol_.clear();
    
    int constraint_offset = 0;
    for (const auto& c : constraints_) {
        if (c->isEnabled()) {
            // 先获取非零元素数量
            int nnz = c->getJacobianStructure(nullptr, nullptr, 0);
            std::vector<int> iRow(nnz), jCol(nnz);
            int actual_nnz = c->getJacobianStructure(iRow.data(), jCol.data(), constraint_offset);
            
            for (int i = 0; i < actual_nnz; ++i) {
                jac_iRow_.push_back(iRow[i]);
                jac_jCol_.push_back(jCol[i]);
            }
            
            constraint_offset += c->getNumConstraints();
        }
    }
    
    jac_structure_computed_ = true;
}

void TrajOptNLP::extractTrajectory(const Ipopt::Number* x_in,
                                   std::vector<double>& Xb, std::vector<double>& Yb,
                                   std::vector<double>& Zb, std::vector<double>& Psib,
                                   std::vector<double>& Xl, std::vector<double>& Yl,
                                   std::vector<double>& Zl,
                                   std::vector<double>& Xr, std::vector<double>& Yr,
                                   std::vector<double>& Zr) const {
    // 如果x_in为nullptr，使用保存的解
    const Ipopt::Number* x = x_in;
    if (x == nullptr) {
        if (!solution_valid_ || solution_.empty()) {
            std::cerr << "Error: No solution available for trajectory extraction!" << std::endl;
            return;
        }
        x = solution_.data();
    }
    // 创建采样点 s ∈ [0, 1]
    std::vector<double> s(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 创建基矩阵
    std::vector<double> Phi(n_sample_ * n_coeff_);
    basis_->computeBasisMatrix(s.data(), n_sample_, n_coeff_, Phi.data());
    
    // 提取系数
    const double* cbx = x + idx_cbx_;
    const double* cby = x + idx_cby_;
    const double* cbz = x + idx_cbz_;
    const double* cpsi = x + idx_cpsi_;
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轨迹
    Xb.resize(n_sample_);
    Yb.resize(n_sample_);
    Zb.resize(n_sample_);
    Psib.resize(n_sample_);
    Xl.resize(n_sample_);
    Yl.resize(n_sample_);
    Zl.resize(n_sample_);
    Xr.resize(n_sample_);
    Yr.resize(n_sample_);
    Zr.resize(n_sample_);
    
    basis_->eval(Phi.data(), n_sample_, cbx, n_coeff_, Xb.data());
    basis_->eval(Phi.data(), n_sample_, cby, n_coeff_, Yb.data());
    basis_->eval(Phi.data(), n_sample_, cbz, n_coeff_, Zb.data());
    basis_->eval(Phi.data(), n_sample_, cpsi, n_coeff_, Psib.data());
    basis_->eval(Phi.data(), n_sample_, clx, n_coeff_, Xl.data());
    basis_->eval(Phi.data(), n_sample_, cly, n_coeff_, Yl.data());
    basis_->eval(Phi.data(), n_sample_, clz, n_coeff_, Zl.data());
    basis_->eval(Phi.data(), n_sample_, crx, n_coeff_, Xr.data());
    basis_->eval(Phi.data(), n_sample_, cry, n_coeff_, Yr.data());
    basis_->eval(Phi.data(), n_sample_, crz, n_coeff_, Zr.data());
}

void TrajOptNLP::extractCoefficients(const Ipopt::Number* x_in,
                                      std::vector<double>& cbx, std::vector<double>& cby,
                                      std::vector<double>& cbz, std::vector<double>& cpsi,
                                      std::vector<double>& clx, std::vector<double>& cly,
                                      std::vector<double>& clz,
                                      std::vector<double>& crx, std::vector<double>& cry,
                                      std::vector<double>& crz) const {
    // 如果x_in为nullptr，使用保存的解
    const Ipopt::Number* x = x_in;
    if (x == nullptr) {
        if (!solution_valid_ || solution_.empty()) {
            std::cerr << "Error: No solution available for coefficient extraction!" << std::endl;
            return;
        }
        x = solution_.data();
    }
    
    // 提取系数
    cbx.resize(n_coeff_);
    cby.resize(n_coeff_);
    cbz.resize(n_coeff_);
    cpsi.resize(n_coeff_);
    clx.resize(n_coeff_);
    cly.resize(n_coeff_);
    clz.resize(n_coeff_);
    crx.resize(n_coeff_);
    cry.resize(n_coeff_);
    crz.resize(n_coeff_);
    
    std::copy(x + idx_cbx_, x + idx_cbx_ + n_coeff_, cbx.begin());
    std::copy(x + idx_cby_, x + idx_cby_ + n_coeff_, cby.begin());
    std::copy(x + idx_cbz_, x + idx_cbz_ + n_coeff_, cbz.begin());
    std::copy(x + idx_cpsi_, x + idx_cpsi_ + n_coeff_, cpsi.begin());
    std::copy(x + idx_clx_, x + idx_clx_ + n_coeff_, clx.begin());
    std::copy(x + idx_cly_, x + idx_cly_ + n_coeff_, cly.begin());
    std::copy(x + idx_clz_, x + idx_clz_ + n_coeff_, clz.begin());
    std::copy(x + idx_crx_, x + idx_crx_ + n_coeff_, crx.begin());
    std::copy(x + idx_cry_, x + idx_cry_ + n_coeff_, cry.begin());
    std::copy(x + idx_crz_, x + idx_crz_ + n_coeff_, crz.begin());
}

void TrajOptNLP::extractContactForces(const Ipopt::Number* x_in,
                                      std::vector<double>& fLx, std::vector<double>& fLy,
                                      std::vector<double>& fLz,
                                      std::vector<double>& fRx, std::vector<double>& fRy,
                                      std::vector<double>& fRz) const {
    // 如果x_in为nullptr，使用保存的解
    const Ipopt::Number* x = x_in;
    if (x == nullptr) {
        if (!solution_valid_ || solution_.empty()) {
            std::cerr << "Error: No solution available for contact force extraction!" << std::endl;
            return;
        }
        x = solution_.data();
    }
    
    // 提取接触力
    fLx.resize(n_sample_);
    fLy.resize(n_sample_);
    fLz.resize(n_sample_);
    fRx.resize(n_sample_);
    fRy.resize(n_sample_);
    fRz.resize(n_sample_);
    
    std::copy(x + idx_fLx_, x + idx_fLx_ + n_sample_, fLx.begin());
    std::copy(x + idx_fLy_, x + idx_fLy_ + n_sample_, fLy.begin());
    std::copy(x + idx_fLz_, x + idx_fLz_ + n_sample_, fLz.begin());
    std::copy(x + idx_fRx_, x + idx_fRx_ + n_sample_, fRx.begin());
    std::copy(x + idx_fRy_, x + idx_fRy_ + n_sample_, fRy.begin());
    std::copy(x + idx_fRz_, x + idx_fRz_ + n_sample_, fRz.begin());
}

