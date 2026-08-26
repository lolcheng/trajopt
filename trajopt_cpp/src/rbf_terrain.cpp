#include "rbf_terrain.hpp"
#include <cmath>
#include <cassert>
#include <algorithm>

RBFTerrain::RBFTerrain(const std::vector<double>& centers_x,
                       const std::vector<double>& centers_y,
                       const std::vector<double>& weights,
                       double sigma)
    : centers_x_(centers_x), centers_y_(centers_y), weights_(weights),
      sigma_(sigma), sigma2_(sigma * sigma) {
    assert(centers_x.size() == centers_y.size());
    assert(centers_x.size() == weights.size());
    assert(sigma > 0.0);
}

double RBFTerrain::height(double x, double y) const {
    double h = 0.0;
    const int n = numCenters();
    
    for (int i = 0; i < n; ++i) {
        double dx = x - centers_x_[i];
        double dy = y - centers_y_[i];
        double r2 = dx * dx + dy * dy;
        h += weights_[i] * std::exp(-0.5 * r2 / sigma2_);
    }
    
    return h;
}

void RBFTerrain::height(const double* x, const double* y, int n, double* h) const {
    for (int j = 0; j < n; ++j) {
        h[j] = height(x[j], y[j]);
    }
}

void RBFTerrain::gradient(double x, double y, double& hx, double& hy) const {
    hx = 0.0;
    hy = 0.0;
    const int n = numCenters();
    
    for (int i = 0; i < n; ++i) {
        double dx = x - centers_x_[i];
        double dy = y - centers_y_[i];
        double r2 = dx * dx + dy * dy;
        double g = weights_[i] * std::exp(-0.5 * r2 / sigma2_);
        hx += g * (-dx / sigma2_);
        hy += g * (-dy / sigma2_);
    }
}

void RBFTerrain::gradient(const double* x, const double* y, int n,
                          double* hx, double* hy) const {
    for (int j = 0; j < n; ++j) {
        gradient(x[j], y[j], hx[j], hy[j]);
    }
}

void RBFTerrain::normal(double x, double y, double& nx, double& ny, double& nz) const {
    double hx, hy;
    gradient(x, y, hx, hy);
    // 法向量为 (-∂H/∂x, -∂H/∂y, 1)，需要归一化
    nx = -hx;
    ny = -hy;
    nz = 1.0;
    
    // 归一化法向量
    double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    const double eps = 1e-8;
    if (norm > eps) {
        nx /= norm;
        ny /= norm;
        nz /= norm;
    } else {
        // 如果法向量太小，使用默认值 (0, 0, 1)
        nx = 0.0;
        ny = 0.0;
        nz = 1.0;
    }
}

void RBFTerrain::normal(const double* x, const double* y, int n,
                       double* nx, double* ny, double* nz) const {
    for (int j = 0; j < n; ++j) {
        normal(x[j], y[j], nx[j], ny[j], nz[j]);
    }
}

void RBFTerrain::getXRange(double& xmin, double& xmax) const {
    if (centers_x_.empty()) {
        xmin = 0.0;
        xmax = 0.0;
        return;
    }
    xmin = *std::min_element(centers_x_.begin(), centers_x_.end());
    xmax = *std::max_element(centers_x_.begin(), centers_x_.end());
}

void RBFTerrain::getYRange(double& ymin, double& ymax) const {
    if (centers_y_.empty()) {
        ymin = 0.0;
        ymax = 0.0;
        return;
    }
    ymin = *std::min_element(centers_y_.begin(), centers_y_.end());
    ymax = *std::max_element(centers_y_.begin(), centers_y_.end());
}

void RBFTerrain::getBounds(double& xmin, double& xmax, double& ymin, double& ymax) const {
    getXRange(xmin, xmax);
    getYRange(ymin, ymax);
}

