#ifndef POLYNOMIAL_BASIS_HPP
#define POLYNOMIAL_BASIS_HPP

#include "basis_function.hpp"

/**
 * 多项式基函数（幂基）
 * 基函数为: [1, s, s^2, s^3, ..., s^(n-1)]
 * 
 * 这是传统的Vandermonde基，也称为幂基或单项式基
 */
class PolynomialBasis : public BasisFunction {
public:
    /**
     * 构造Vandermonde矩阵（多项式基矩阵）
     * Phi[i, j] = s[i]^j, 其中 j = 0, 1, ..., n_coeff-1
     */
    void computeBasisMatrix(const double* s, int n_samples, int n_coeff,
                            double* Phi) override;
    
    /**
     * 评估轨迹在单个点的值: x(s) = [1, s, s^2, ..., s^(n_coeff-1)] * coeff
     */
    double evalPoint(double s, const double* coeff, int n_coeff) override;
    
    /**
     * 评估轨迹: x(s) = Phi * coeff
     */
    void eval(const double* Phi, int n_samples, const double* coeff,
             int n_coeff, double* x) override;
    
    /**
     * 评估轨迹导数: dx/ds = [0, 1, 2s, 3s^2, ..., (n_coeff-1)s^(n_coeff-2)] * coeff
     */
    void evalDerivative(const double* s, int n_samples,
                       const double* coeff, int n_coeff, double* dx) override;
    
    /**
     * 评估轨迹导数在单个点的值
     */
    double evalDerivativePoint(double s, const double* coeff, int n_coeff) override;
};

#endif // POLYNOMIAL_BASIS_HPP

