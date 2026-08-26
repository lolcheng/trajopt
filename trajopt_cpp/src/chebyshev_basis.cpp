#include "chebyshev_basis.hpp"
#include <cmath>
#include <cassert>

double ChebyshevBasis::chebyshevT(int n, double x) {
    if (n == 0) {
        return 1.0;
    } else if (n == 1) {
        return x;
    } else {
        // 使用递推关系: T_{n+1}(x) = 2x * T_n(x) - T_{n-1}(x)
        double T_prev = 1.0;  // T_0(x) = 1
        double T_curr = x;    // T_1(x) = x
        for (int i = 2; i <= n; ++i) {
            double T_next = 2.0 * x * T_curr - T_prev;
            T_prev = T_curr;
            T_curr = T_next;
        }
        return T_curr;
    }
}

double ChebyshevBasis::chebyshevTDerivative(int n, double x) {
    if (n == 0) {
        return 0.0;
    } else if (n == 1) {
        return 1.0;
    } else {
        // 切比雪夫多项式的导数公式:
        // dT_n(x)/dx = n * U_{n-1}(x)，其中 U_n 是第二类切比雪夫多项式
        // 但更简单的方法是使用递推关系:
        // dT_n(x)/dx = n * (T_{n-1}(x) - x*T_n(x)) / (1 - x^2)  (当 x ≠ ±1)
        // 或者使用: dT_n(x)/dx = n * sum_{k=0}^{floor((n-1)/2)} (-1)^k * (n-k-1)! / (k! * (n-2k-1)!) * 2^{n-2k-1} * x^{n-2k-1}
        // 
        // 更实用的方法是使用递推关系计算 U_{n-1}(x):
        // U_0(x) = 1
        // U_1(x) = 2x
        // U_{n+1}(x) = 2x * U_n(x) - U_{n-1}(x)
        // 然后 dT_n(x)/dx = n * U_{n-1}(x)
        
        // 特殊情况：x = ±1
        if (std::abs(x - 1.0) < 1e-10) {
            return n * n;  // dT_n(1)/dx = n^2
        } else if (std::abs(x + 1.0) < 1e-10) {
            return n * n * ((n % 2 == 0) ? 1.0 : -1.0);  // dT_n(-1)/dx = (-1)^{n+1} * n^2
        }
        
        // 计算 U_{n-1}(x)
        if (n == 1) {
            return 1.0;  // dT_1(x)/dx = 1
        }
        
        double U_prev = 1.0;  // U_0(x) = 1
        double U_curr = 2.0 * x;  // U_1(x) = 2x
        for (int i = 2; i < n; ++i) {
            double U_next = 2.0 * x * U_curr - U_prev;
            U_prev = U_curr;
            U_curr = U_next;
        }
        
        return n * U_curr;  // dT_n(x)/dx = n * U_{n-1}(x)
    }
}

void ChebyshevBasis::computeBasisMatrix(const double* s, int n_samples, int n_coeff,
                                       double* Phi) {
    for (int i = 0; i < n_samples; ++i) {
        double s_val = s[i];
        // 将 s ∈ [0, 1] 映射到 x ∈ [-1, 1]
        double x = 2.0 * s_val - 1.0;
        
        // 计算所有阶数的切比雪夫多项式值
        for (int j = 0; j < n_coeff; ++j) {
            Phi[i * n_coeff + j] = chebyshevT(j, x);
        }
    }
}

void ChebyshevBasis::eval(const double* Phi, int n_samples, const double* coeff,
                         int n_coeff, double* x) {
    for (int i = 0; i < n_samples; ++i) {
        x[i] = 0.0;
        for (int j = 0; j < n_coeff; ++j) {
            x[i] += Phi[i * n_coeff + j] * coeff[j];
        }
    }
}

double ChebyshevBasis::evalPoint(double s, const double* coeff, int n_coeff) {
    // 将 s ∈ [0, 1] 映射到 x ∈ [-1, 1]
    double x = 2.0 * s - 1.0;
    
    double result = 0.0;
    for (int j = 0; j < n_coeff; ++j) {
        result += coeff[j] * chebyshevT(j, x);
    }
    return result;
}

void ChebyshevBasis::evalDerivative(const double* s, int n_samples,
                                    const double* coeff, int n_coeff, double* dx) {
    for (int i = 0; i < n_samples; ++i) {
        double s_val = s[i];
        // 将 s ∈ [0, 1] 映射到 x ∈ [-1, 1]
        double x = 2.0 * s_val - 1.0;
        
        dx[i] = 0.0;
        for (int j = 0; j < n_coeff; ++j) {
            // dT_j(2s-1)/ds = 2 * dT_j(x)/dx
            double dT_dx = chebyshevTDerivative(j, x);
            dx[i] += coeff[j] * 2.0 * dT_dx;
        }
    }
}

double ChebyshevBasis::evalDerivativePoint(double s, const double* coeff, int n_coeff) {
    // 将 s ∈ [0, 1] 映射到 x ∈ [-1, 1]
    double x = 2.0 * s - 1.0;
    
    double result = 0.0;
    for (int j = 0; j < n_coeff; ++j) {
        // dT_j(2s-1)/ds = 2 * dT_j(x)/dx
        double dT_dx = chebyshevTDerivative(j, x);
        result += coeff[j] * 2.0 * dT_dx;
    }
    return result;
}

double ChebyshevBasis::chebyshevTSecondDerivative(int n, double x) {
    if (n <= 1)
        return 0.0;
    // 切比雪夫方程: (1-x^2) T_n'' - x T_n' + n^2 T_n = 0
    // => T_n'' = (x T_n' - n^2 T_n) / (1 - x^2)
    if (std::abs(x - 1.0) < 1e-10) {
        return n * n * (static_cast<double>(n * n) - 1.0) / 3.0;
    }
    if (std::abs(x + 1.0) < 1e-10) {
        const double sign = (n % 2 == 0) ? 1.0 : -1.0;
        return sign * n * n * (static_cast<double>(n * n) - 1.0) / 3.0;
    }
    const double T = chebyshevT(n, x);
    const double dT = chebyshevTDerivative(n, x);
    const double one_mx2 = 1.0 - x * x;
    return (x * dT - static_cast<double>(n * n) * T) / one_mx2;
}

void ChebyshevBasis::fillDerivativeBasisAt(double s, int n_coeff,
                                           double* dphi_ds, double* d2phi_ds2) {
    assert(dphi_ds != nullptr && d2phi_ds2 != nullptr);
    const double x = 2.0 * s - 1.0;
    for (int j = 0; j < n_coeff; ++j) {
        const double dT_dx = chebyshevTDerivative(j, x);
        const double d2T_dx2 = chebyshevTSecondDerivative(j, x);
        // phi_j(s) = T_j(x), x = 2s-1 => d/ds = 2 d/dx, d^2/ds^2 = 4 d^2/dx^2
        dphi_ds[j] = 2.0 * dT_dx;
        d2phi_ds2[j] = 4.0 * d2T_dx2;
    }
}

