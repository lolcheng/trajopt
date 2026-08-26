# 基函数系统说明

## 概述

代码现在支持多种基函数，通过继承 `BasisFunction` 接口类实现。当前实现了两种基函数：

1. **多项式基（PolynomialBasis）**：传统的幂基 `[1, s, s^2, s^3, ..., s^(n-1)]`
2. **切比雪夫基（ChebyshevBasis）**：第一类切比雪夫多项式基

## 什么是切比雪夫基？

### 切比雪夫多项式

切比雪夫多项式是定义在区间 `[-1, 1]` 上的一组正交多项式：

- **T_0(x) = 1**
- **T_1(x) = x**
- **T_{n+1}(x) = 2x * T_n(x) - T_{n-1}(x)**（递推关系）

### 在轨迹优化中的应用

由于我们的参数 `s ∈ [0, 1]`，需要将 `s` 映射到 `[-1, 1]`：
```
x = 2s - 1
```

基函数为：`[T_0(2s-1), T_1(2s-1), T_2(2s-1), ..., T_{n-1}(2s-1)]`

### 切比雪夫基的优势

1. **数值稳定性更好**：
   - 高阶项不会像幂基那样快速增长
   - 避免了 Vandermonde 矩阵的病态问题

2. **近似误差更均匀**：
   - 在区间内误差分布更均匀（等波纹性质）
   - 最大误差最小化（Chebyshev 交替定理）

3. **条件数更小**：
   - 基矩阵的条件数通常比 Vandermonde 矩阵小得多
   - 使得优化问题更容易求解

4. **更好的收敛性**：
   - 对于光滑函数，切比雪夫基通常需要更少的项就能达到相同的精度

## 使用方法

### 在 main.cpp 中选择基函数

```cpp
// 基函数选择
const std::string BASIS_TYPE = "chebyshev";  // 或 "polynomial"

std::shared_ptr<BasisFunction> basis;
if (BASIS_TYPE == "chebyshev") {
    basis = std::make_shared<ChebyshevBasis>();
    std::cout << "Using Chebyshev basis functions" << std::endl;
} else {
    basis = std::make_shared<PolynomialBasis>();
    std::cout << "Using Polynomial basis functions (power basis)" << std::endl;
}

// 创建NLP时传入基函数
Ipopt::SmartPtr<TrajOptNLP> nlp = new TrajOptNLP(N_COEFF, N_SAMPLE, basis);
```

### 向后兼容性

为了保持向后兼容，`Polynomial` 类仍然存在，但现在是一个静态包装类，内部使用 `PolynomialBasis`。所有使用 `Polynomial::` 的代码仍然可以正常工作。

## 基函数接口

所有基函数都继承自 `BasisFunction` 接口，需要实现以下方法：

```cpp
class BasisFunction {
public:
    // 构造基矩阵
    virtual void computeBasisMatrix(const double* s, int n_samples, int n_coeff,
                                    double* Phi) = 0;
    
    // 评估轨迹在单个点的值
    virtual double evalPoint(double s, const double* coeff, int n_coeff) = 0;
    
    // 评估轨迹
    virtual void eval(const double* Phi, int n_samples, const double* coeff,
                     int n_coeff, double* x) = 0;
    
    // 评估轨迹导数
    virtual void evalDerivative(const double* s, int n_samples,
                               const double* coeff, int n_coeff, double* dx) = 0;
    
    // 评估轨迹导数在单个点的值
    virtual double evalDerivativePoint(double s, const double* coeff, int n_coeff) = 0;
};
```

## 添加新的基函数

要添加新的基函数（例如 B样条基、Legendre 基等），只需：

1. 创建新类继承自 `BasisFunction`
2. 实现所有纯虚函数
3. 在 `main.cpp` 中添加选择逻辑

示例：

```cpp
class LegendreBasis : public BasisFunction {
    // 实现所有接口方法
};
```

## 文件结构

```
include/
├── basis_function.hpp          # 基函数接口
├── polynomial_basis.hpp         # 多项式基实现
├── chebyshev_basis.hpp          # 切比雪夫基实现
└── polynomial.hpp               # 向后兼容的包装类

src/
├── polynomial_basis.cpp
├── chebyshev_basis.cpp
└── polynomial.cpp               # 向后兼容的包装实现
```

## 性能对比

| 特性 | 多项式基 | 切比雪夫基 |
|------|---------|-----------|
| 数值稳定性 | 较差（高阶项快速增长） | 较好 |
| 条件数 | 大（Vandermonde 矩阵） | 小 |
| 近似误差分布 | 不均匀 | 均匀（等波纹） |
| 实现复杂度 | 简单 | 中等 |
| 计算开销 | 低 | 略高（递推计算） |

## 建议

- **默认使用切比雪夫基**：对于大多数轨迹优化问题，切比雪夫基能提供更好的数值稳定性和收敛性
- **多项式基保留**：作为向后兼容和简单测试使用
- **根据问题选择**：如果遇到特定问题，可以尝试不同的基函数

