# C++ 架构设计文档

## 设计原则

1. **模块化**: 每个约束/目标项都是独立的类，可以单独启用/禁用
2. **可扩展**: 方便添加新的约束或目标项
3. **高效**: 使用C++实现，避免符号图构建
4. **清晰**: 接口明确，易于理解和维护

## 目录结构

```
trajopt_cpp/
├── include/
│   ├── trajopt_nlp.hpp          # 主NLP问题类
│   ├── constraint_base.hpp       # 约束基类
│   ├── objective_base.hpp       # 目标项基类
│   ├── rbf_terrain.hpp          # RBF地形类
│   ├── polynomial.hpp            # 多项式工具
│   └── constraints/              # 各约束模块
│       ├── boundary_constraint.hpp
│       ├── wheel_base_constraint.hpp
│       ├── wheel_distance_constraint.hpp
│       ├── anti_cross_constraint.hpp
│       └── friction_constraint.hpp
│   └── objectives/              # 各目标项模块
│       ├── smoothing_objective.hpp
│       ├── regularization_objective.hpp
│       ├── yaw_alignment_objective.hpp
│       └── ...
├── src/
│   ├── trajopt_nlp.cpp
│   ├── rbf_terrain.cpp
│   ├── polynomial.cpp
│   └── constraints/
│       └── ...
└── main.cpp
```

## 核心类设计

### 1. ConstraintBase (约束基类)

```cpp
class ConstraintBase {
public:
    virtual ~ConstraintBase() = default;
    
    // 返回约束数量
    virtual int getNumConstraints() const = 0;
    
    // 评估约束值 g(x)
    virtual void eval_g(const double* x, double* g) = 0;
    
    // 评估约束雅可比 J = ∂g/∂x
    virtual void eval_jac_g(const double* x, double* jac_g) = 0;
    
    // 获取约束边界 [g_lb, g_ub]
    virtual void getBounds(double* g_lb, double* g_ub) = 0;
    
    // 获取雅可比稀疏结构
    virtual void getJacobianStructure(int* iRow, int* jCol) = 0;
};
```

### 2. ObjectiveBase (目标项基类)

```cpp
class ObjectiveBase {
public:
    virtual ~ObjectiveBase() = default;
    
    // 评估目标项值
    virtual double eval_f(const double* x) = 0;
    
    // 评估目标项梯度
    virtual void eval_grad_f(const double* x, double* grad_f) = 0;
    
    // 评估目标项海塞矩阵贡献 (可选)
    virtual void eval_hessian(const double* x, const double* lambda,
                              double* hessian) = 0;
    
    // 是否启用
    virtual bool isEnabled() const = 0;
};
```

### 3. TrajOptNLP (主NLP问题类)

```cpp
class TrajOptNLP : public Ipopt::TNLP {
public:
    // 添加约束
    void addConstraint(std::shared_ptr<ConstraintBase> constraint);
    
    // 添加目标项
    void addObjective(std::shared_ptr<ObjectiveBase> objective);
    
    // 设置RBF地形
    void setTerrain(std::shared_ptr<RBFTerrain> terrain);
    
    // Ipopt接口实现
    bool get_nlp_info(...) override;
    bool get_bounds_info(...) override;
    bool eval_f(...) override;
    bool eval_grad_f(...) override;
    bool eval_g(...) override;
    bool eval_jac_g(...) override;
    bool eval_h(...) override;
    
private:
    std::vector<std::shared_ptr<ConstraintBase>> constraints_;
    std::vector<std::shared_ptr<ObjectiveBase>> objectives_;
    std::shared_ptr<RBFTerrain> terrain_;
    
    // 变量布局信息
    int n_vars_;
    int n_constraints_;
    VariableLayout layout_;
};
```

### 4. RBFTerrain (RBF地形类)

```cpp
class RBFTerrain {
public:
    RBFTerrain(const std::vector<double>& centers_x,
               const std::vector<double>& centers_y,
               const std::vector<double>& weights,
               double sigma);
    
    // 计算地形高度
    double height(double x, double y) const;
    
    // 计算地形高度 (批量)
    void height(const double* x, const double* y, int n, double* h) const;
    
    // 计算梯度
    void gradient(double x, double y, double& hx, double& hy) const;
    
    // 计算梯度 (批量)
    void gradient(const double* x, const double* y, int n,
                  double* hx, double* hy) const;
    
    // 计算表面法向量
    void normal(double x, double y, double& nx, double& ny, double& nz) const;
    
private:
    std::vector<double> centers_x_, centers_y_, weights_;
    double sigma_;
    double sigma2_;
};
```

### 5. Polynomial (多项式工具类)

```cpp
class Polynomial {
public:
    // 构造多项式基矩阵 Phi
    static void vandermonde(const double* s, int n_samples, int n_coeff,
                           double* Phi);
    
    // 评估轨迹: x(s) = Phi * coeff
    static void eval(const double* s, int n_samples, const double* coeff,
                    int n_coeff, double* x);
    
    // 评估轨迹导数: dx/ds = dPhi/ds * coeff
    static void eval_derivative(const double* s, int n_samples,
                               const double* coeff, int n_coeff, double* dx);
};
```

## 变量布局

```cpp
struct VariableLayout {
    // 多项式系数索引
    int idx_cbx, idx_cby, idx_cbz, idx_cpsi;
    int idx_clx, idx_cly, idx_clz;
    int idx_crx, idx_cry, idx_crz;
    
    // 接触力索引 (每个采样点)
    int idx_fLx, idx_fLy, idx_fLz;
    int idx_fRx, idx_fRy, idx_fRz;
    
    // 参数
    int n_coeff;
    int n_sample;
    int n_vars;
};
```

## 约束实现示例

### BoundaryConstraint

```cpp
class BoundaryConstraint : public ConstraintBase {
public:
    BoundaryConstraint(const VariableLayout& layout,
                      double x0, double y0, double z0, double psi0,
                      double x1, double y1, double z1, double psi1,
                      const RBFTerrain& terrain);
    
    int getNumConstraints() const override { return 8; }
    
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* jac_g) override;
    void getBounds(double* g_lb, double* g_ub) override;
    
private:
    VariableLayout layout_;
    double x0_, y0_, z0_, psi0_;
    double x1_, y1_, z1_, psi1_;
    const RBFTerrain& terrain_;
};
```

## 使用流程

```cpp
// 1. 创建NLP问题
auto nlp = std::make_shared<TrajOptNLP>();

// 2. 设置地形
auto terrain = std::make_shared<RBFTerrain>(centers_x, centers_y, weights, sigma);
nlp->setTerrain(terrain);

// 3. 添加约束
nlp->addConstraint(std::make_shared<BoundaryConstraint>(...));
nlp->addConstraint(std::make_shared<WheelBaseConstraint>(...));
// ...

// 4. 添加目标项
nlp->addObjective(std::make_shared<SmoothingObjective>(...));
// ...

// 5. 求解
Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
app->Options()->SetStringValue("linear_solver", "mumps");
app->Initialize();
Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(nlp);
```

## 扩展性

添加新约束只需：
1. 继承 `ConstraintBase`
2. 实现接口方法
3. 在 `TrajOptNLP` 中注册

添加新目标项只需：
1. 继承 `ObjectiveBase`
2. 实现接口方法
3. 在 `TrajOptNLP` 中注册

