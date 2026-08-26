#ifndef TRAJOPT_NLP_HPP
#define TRAJOPT_NLP_HPP

#include <stddef.h>  // Must be before Ipopt headers (C-style header)
#include <vector>
#include <memory>
#include "IpTNLP.hpp"
#include "IpSmartPtr.hpp"
#include "constraint_base.hpp"
#include "objective_base.hpp"
#include "rbf_terrain.hpp"
#include "basis_function.hpp"
#include <memory>

/**
 * 轨迹优化NLP问题类
 * 实现Ipopt的TNLP接口
 */
class TrajOptNLP : public Ipopt::TNLP {
public:
    /**
     * 构造函数
     * @param n_coeff 基函数系数数量
     * @param n_sample 采样点数量
     * @param basis 基函数指针（如果为nullptr，则默认使用切比雪夫基）
     */
    TrajOptNLP(int n_coeff, int n_sample, 
               std::shared_ptr<BasisFunction> basis = nullptr);
    
    virtual ~TrajOptNLP();
    
    /**
     * 添加约束
     */
    void addConstraint(std::shared_ptr<ConstraintBase> constraint);
    
    /**
     * 添加目标函数项
     */
    void addObjective(std::shared_ptr<ObjectiveBase> objective);
    
    /**
     * 设置RBF地形
     */
    void setTerrain(std::shared_ptr<RBFTerrain> terrain);
    
    /**
     * 获取地形
     */
    std::shared_ptr<RBFTerrain> getTerrain() const { return terrain_; }
    
    /**
     * 从解中提取轨迹数据
     * @param x 优化变量（解），如果为nullptr则使用保存的解
     * @param Xb 输出的基座x轨迹
     * @param Yb 输出的基座y轨迹
     * @param Zb 输出的基座z轨迹
     * @param Psib 输出的基座偏航角轨迹
     * @param Xl 输出的左轮x轨迹
     * @param Yl 输出的左轮y轨迹
     * @param Zl 输出的左轮z轨迹
     * @param Xr 输出的右轮x轨迹
     * @param Yr 输出的右轮y轨迹
     * @param Zr 输出的右轮z轨迹
     */
    void extractTrajectory(const Ipopt::Number* x,
                          std::vector<double>& Xb, std::vector<double>& Yb,
                          std::vector<double>& Zb, std::vector<double>& Psib,
                          std::vector<double>& Xl, std::vector<double>& Yl,
                          std::vector<double>& Zl,
                          std::vector<double>& Xr, std::vector<double>& Yr,
                          std::vector<double>& Zr) const;
    
    /**
     * 从解中提取轨迹系数
     * @param x 优化变量（解），如果为nullptr则使用保存的解
     * @param cbx 输出的基座x系数
     * @param cby 输出的基座y系数
     * @param cbz 输出的基座z系数
     * @param cpsi 输出的基座偏航角系数
     * @param clx 输出的左轮x系数
     * @param cly 输出的左轮y系数
     * @param clz 输出的左轮z系数
     * @param crx 输出的右轮x系数
     * @param cry 输出的右轮y系数
     * @param crz 输出的右轮z系数
     */
    void extractCoefficients(const Ipopt::Number* x,
                            std::vector<double>& cbx, std::vector<double>& cby,
                            std::vector<double>& cbz, std::vector<double>& cpsi,
                            std::vector<double>& clx, std::vector<double>& cly,
                            std::vector<double>& clz,
                            std::vector<double>& crx, std::vector<double>& cry,
                            std::vector<double>& crz) const;
    
    /**
     * 从解中提取接触力
     * @param x 优化变量（解），如果为nullptr则使用保存的解
     * @param fLx 输出的左轮x方向力
     * @param fLy 输出的左轮y方向力
     * @param fLz 输出的左轮z方向力
     * @param fRx 输出的右轮x方向力
     * @param fRy 输出的右轮y方向力
     * @param fRz 输出的右轮z方向力
     */
    void extractContactForces(const Ipopt::Number* x,
                             std::vector<double>& fLx, std::vector<double>& fLy,
                             std::vector<double>& fLz,
                             std::vector<double>& fRx, std::vector<double>& fRy,
                             std::vector<double>& fRz) const;
    
    /**
     * 检查是否有有效的解
     */
    bool hasSolution() const { return solution_valid_; }
    
    /**
     * 获取保存的解（用于可视化）
     */
    const std::vector<Ipopt::Number>& getSolution() const { return solution_; }
    
    /**
     * 获取求解器返回的状态码（SolverReturn）
     * 注意：这与ApplicationReturnStatus可能不同
     */
    Ipopt::SolverReturn getSolverStatus() const { return solver_status_; }
    
    /**
     * 获取变量总数
     */
    int getNumVars() const { return n_vars_; }
    
    /**
     * 获取约束总数
     */
    int getNumConstraints() const { return n_constraints_; }
    
    /**
     * 获取约束总数（Ipopt接口）
     */
    Ipopt::Index get_n_constraints() const { return static_cast<Ipopt::Index>(n_constraints_); }
    
    /**
     * 设置初始猜测
     * @param x 初始猜测数组（长度应为n_vars_）
     */
    void setInitialGuess(const Ipopt::Number* x);
    
    /**
     * 获取初始猜测
     * @return 初始猜测数组的指针（如果未设置则返回nullptr）
     */
    const Ipopt::Number* getInitialGuess() const { return initial_guess_.empty() ? nullptr : initial_guess_.data(); }
    
    // Ipopt接口实现
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                      Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style) override;
    
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;
    
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                           bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                           Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;
    
    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number& obj_value) override;
    
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                     Ipopt::Number* grad_f) override;
    
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
    std::shared_ptr<BasisFunction> basis_;  // 基函数指针
    
    // 变量索引（多项式系数）
    int idx_cbx_, idx_cby_, idx_cbz_, idx_cpsi_;
    int idx_clx_, idx_cly_, idx_clz_;
    int idx_crx_, idx_cry_, idx_crz_;
    
    // 变量索引（接触力）
    int idx_fLx_, idx_fLy_, idx_fLz_;
    int idx_fRx_, idx_fRy_, idx_fRz_;
    
    // 计算变量总数和约束总数
    void updateDimensions();
    
    // 计算雅可比稀疏结构
    void computeJacobianStructure();
    std::vector<Ipopt::Index> jac_iRow_;
    std::vector<Ipopt::Index> jac_jCol_;
    bool jac_structure_computed_;
    
    // 保存的解
    std::vector<Ipopt::Number> solution_;
    bool solution_valid_;
    
    // 求解器状态码（SolverReturn，与ApplicationReturnStatus可能不同）
    Ipopt::SolverReturn solver_status_;
    
    // 初始猜测
    std::vector<Ipopt::Number> initial_guess_;
};

#endif // TRAJOPT_NLP_HPP

