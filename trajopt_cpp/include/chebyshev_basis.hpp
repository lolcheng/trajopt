#ifndef CHEBYSHEV_BASIS_HPP
#define CHEBYSHEV_BASIS_HPP

#include "basis_function.hpp"
#include <cmath>

/**
 * 切比雪夫基函数（第一类切比雪夫多项式）
 * 
 * 切比雪夫多项式定义在区间 [-1, 1] 上：
 *   T_0(x) = 1
 *   T_1(x) = x
 *   T_{n+1}(x) = 2x * T_n(x) - T_{n-1}(x)
 * 
 * 由于我们的参数 s ∈ [0, 1]，需要将 s 映射到 [-1, 1]：
 *   x = 2s - 1
 * 
 * 基函数为: [T_0(2s-1), T_1(2s-1), T_2(2s-1), ..., T_{n-1}(2s-1)]
 * 
 * 切比雪夫基的优势：
 * 1. 数值稳定性更好：高阶项不会像幂基那样快速增长
 * 2. 近似误差更均匀：在区间内误差分布更均匀（等波纹性质）
 * 3. 条件数更小：基矩阵的条件数通常比Vandermonde矩阵小得多
 */
class ChebyshevBasis : public BasisFunction {
public:
    /**
     * 构造切比雪夫基矩阵
     * Phi[i, j] = T_j(2*s[i] - 1)，其中 T_j 是第 j 阶切比雪夫多项式
     */
    void computeBasisMatrix(const double* s, int n_samples, int n_coeff,
                            double* Phi) override;
    
    /**
     * 评估轨迹在单个点的值: x(s) = sum(coeff[j] * T_j(2s-1))
     */
    double evalPoint(double s, const double* coeff, int n_coeff) override;
    
    /**
     * 评估轨迹: x(s) = Phi * coeff
     */
    void eval(const double* Phi, int n_samples, const double* coeff,
             int n_coeff, double* x) override;
    
    /**
     * 评估轨迹导数: dx/ds = sum(coeff[j] * dT_j(2s-1)/ds)
     * 注意: dT_j(2s-1)/ds = 2 * dT_j(x)/dx，其中 x = 2s-1
     */
    void evalDerivative(const double* s, int n_samples,
                       const double* coeff, int n_coeff, double* dx) override;
    
    /**
     * 评估轨迹导数在单个点的值
     */
    double evalDerivativePoint(double s, const double* coeff, int n_coeff) override;

    void fillDerivativeBasisAt(double s, int n_coeff,
                               double* dphi_ds, double* d2phi_ds2) override;

private:
    /**
     * 计算第 n 阶切比雪夫多项式在 x 处的值
     * @param n 多项式阶数
     * @param x 自变量，x ∈ [-1, 1]
     * @return T_n(x)
     */
    static double chebyshevT(int n, double x);
    
    /**
     * 计算第 n 阶切比雪夫多项式的导数在 x 处的值
     * @param n 多项式阶数
     * @param x 自变量，x ∈ [-1, 1]
     * @return dT_n(x)/dx
     */
    static double chebyshevTDerivative(int n, double x);

    /** d^2 T_n(x) / dx^2，x ∈ [-1,1]（用于 fillDerivativeBasisAt） */
    static double chebyshevTSecondDerivative(int n, double x);
};

#endif // CHEBYSHEV_BASIS_HPP

