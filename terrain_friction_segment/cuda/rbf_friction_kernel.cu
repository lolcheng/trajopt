/**
 * GPU 内核：对网格点 (x,y) 计算 RBF 法向量，并判断摩擦锥约束
 * 约束：竖直支撑力 f=[0,0,mg] 在锥内等价于 n_z >= 1/sqrt(1+mu^2)
 */

#include <cuda_runtime.h>
#include <cmath>

namespace {

__device__ void rbf_normal_at(double x, double y,
                              const double* cx, const double* cy, const double* w,
                              int n_centers, double sigma2,
                              double& nx, double& ny, double& nz) {
    double hx = 0.0, hy = 0.0;
    for (int i = 0; i < n_centers; ++i) {
        double dx = x - cx[i];
        double dy = y - cy[i];
        double r2 = dx * dx + dy * dy;
        double g = w[i] * exp(-0.5 * r2 / sigma2);
        hx += g * (-dx / sigma2);
        hy += g * (-dy / sigma2);
    }
    nx = -hx;
    ny = -hy;
    nz = 1.0;
    double norm = sqrt(nx * nx + ny * ny + nz * nz);
    const double eps = 1e-8;
    if (norm > eps) {
        nx /= norm;
        ny /= norm;
        nz /= norm;
    } else {
        nx = 0.0;
        ny = 0.0;
        nz = 1.0;
    }
}

__global__ void friction_cone_kernel(
    const double* __restrict__ centers_x,
    const double* __restrict__ centers_y,
    const double* __restrict__ weights,
    int n_centers,
    double sigma2,
    double xmin, double xmax, double ymin, double ymax,
    int nx, int ny,
    double nz_min,
    int* __restrict__ labels)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) return;

    double x = (nx > 1) ? (xmin + (xmax - xmin) * (double)i / (double)(nx - 1)) : xmin;
    double y = (ny > 1) ? (ymin + (ymax - ymin) * (double)j / (double)(ny - 1)) : ymin;

    double nx_val, ny_val, nz_val;
    rbf_normal_at(x, y, centers_x, centers_y, weights, n_centers, sigma2,
                  nx_val, ny_val, nz_val);

    int idx = j * nx + i;
    labels[idx] = (nz_val >= nz_min) ? 1 : 0;
}

} // anonymous namespace

extern "C" {

void launch_friction_cone(
    const double* d_cx, const double* d_cy, const double* d_w,
    int n_centers, double sigma2,
    double xmin, double xmax, double ymin, double ymax,
    int nx, int ny, double nz_min,
    int* d_labels)
{
    dim3 block(16, 16);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    friction_cone_kernel<<<grid, block>>>(
        d_cx, d_cy, d_w, n_centers, sigma2,
        xmin, xmax, ymin, ymax, nx, ny, nz_min, d_labels);
}

} // extern "C"
