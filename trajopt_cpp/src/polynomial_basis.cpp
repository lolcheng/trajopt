#include "polynomial_basis.hpp"
#include <cmath>
#include <cassert>

void PolynomialBasis::computeBasisMatrix(const double* s, int n_samples, int n_coeff,
                                         double* Phi) {
    for (int i = 0; i < n_samples; ++i) {
        double s_val = s[i];
        double s_power = 1.0;  // s^0 = 1
        for (int j = 0; j < n_coeff; ++j) {
            Phi[i * n_coeff + j] = s_power;
            s_power *= s_val;  // s^j
        }
    }
}

void PolynomialBasis::eval(const double* Phi, int n_samples, const double* coeff,
                          int n_coeff, double* x) {
    for (int i = 0; i < n_samples; ++i) {
        x[i] = 0.0;
        for (int j = 0; j < n_coeff; ++j) {
            x[i] += Phi[i * n_coeff + j] * coeff[j];
        }
    }
}

double PolynomialBasis::evalPoint(double s, const double* coeff, int n_coeff) {
    double result = 0.0;
    double s_power = 1.0;
    for (int j = 0; j < n_coeff; ++j) {
        result += coeff[j] * s_power;
        s_power *= s;
    }
    return result;
}

void PolynomialBasis::evalDerivative(const double* s, int n_samples,
                                     const double* coeff, int n_coeff, double* dx) {
    for (int i = 0; i < n_samples; ++i) {
        double s_val = s[i];
        dx[i] = 0.0;
        double s_power = 1.0;  // 从s^0开始
        for (int j = 1; j < n_coeff; ++j) {  // 跳过常数项
            dx[i] += j * coeff[j] * s_power;
            s_power *= s_val;
        }
    }
}

double PolynomialBasis::evalDerivativePoint(double s, const double* coeff, int n_coeff) {
    double result = 0.0;
    double s_power = 1.0;  // 从s^0开始
    for (int j = 1; j < n_coeff; ++j) {  // 跳过常数项
        result += j * coeff[j] * s_power;
        s_power *= s;
    }
    return result;
}

