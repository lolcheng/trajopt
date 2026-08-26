#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "basis_function.hpp"
#include "polynomial_basis.hpp"
#include "chebyshev_basis.hpp"
#include <memory>

/**
 * 多项式工具类（向后兼容包装）
 * 用于轨迹的多项式表示和求值
 * 
 * 注意：这个类现在是一个静态包装类，内部使用可配置基函数（默认切比雪夫基）
 * 为了向后兼容，保留了原有的静态接口
 */
class Polynomial {
private:
    static std::shared_ptr<BasisFunction> basis_;  // 当前全局基函数
    
public:
    /**
     * 设置全局基函数（供所有仍使用 Polynomial:: 接口的模块共享）
     * @param basis 基函数；若为空则回退到默认切比雪夫基
     */
    static void set_basis(std::shared_ptr<BasisFunction> basis) {
        if (basis) {
            basis_ = std::move(basis);
        } else {
            basis_ = std::make_shared<ChebyshevBasis>();
        }
    }

    /**
     * 构造Vandermonde矩阵（多项式基矩阵）
     * Phi[i, j] = s[i]^j, 其中 j = 0, 1, ..., n_coeff-1
     * 
     * @param s 采样点数组 (n_samples)
     * @param n_samples 采样点数量
     * @param n_coeff 多项式系数数量（阶数+1）
     * @param Phi 输出的基矩阵 (n_samples × n_coeff)，按行优先存储
     */
    static void vandermonde(const double* s, int n_samples, int n_coeff,
                           double* Phi) {
        basis_->computeBasisMatrix(s, n_samples, n_coeff, Phi);
    }
    
    /**
     * 评估轨迹: x(s) = Phi * coeff
     * 
     * @param Phi 基矩阵 (n_samples × n_coeff)，按行优先存储
     * @param n_samples 采样点数量
     * @param coeff 多项式系数 (n_coeff)
     * @param n_coeff 系数数量
     * @param x 输出的轨迹值 (n_samples)
     */
    static void eval(const double* Phi, int n_samples, const double* coeff,
                    int n_coeff, double* x) {
        basis_->eval(Phi, n_samples, coeff, n_coeff, x);
    }
    
    /**
     * 评估轨迹在单个点的值: x(s) = [1, s, s^2, ..., s^(n_coeff-1)] * coeff
     * 
     * @param s 采样点
     * @param coeff 多项式系数
     * @param n_coeff 系数数量
     * @return 轨迹值
     */
    static double eval_point(double s, const double* coeff, int n_coeff) {
        return basis_->evalPoint(s, coeff, n_coeff);
    }
    
    /**
     * 评估轨迹导数: dx/ds = dPhi/ds * coeff
     * 
     * @param s 采样点数组
     * @param n_samples 采样点数量
     * @param coeff 多项式系数
     * @param n_coeff 系数数量
     * @param dx 输出的导数值 (n_samples)
     */
    static void eval_derivative(const double* s, int n_samples,
                               const double* coeff, int n_coeff, double* dx) {
        basis_->evalDerivative(s, n_samples, coeff, n_coeff, dx);
    }
};

#endif // POLYNOMIAL_HPP

