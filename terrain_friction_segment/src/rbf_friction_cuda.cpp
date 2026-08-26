#include "rbf_friction_cuda.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cuda_runtime.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace rbf_friction_cuda {

bool run(const std::vector<double>& centers_x,
         const std::vector<double>& centers_y,
         const std::vector<double>& weights,
         double sigma,
         const Params& params,
         std::vector<int>& out_labels) {
    const int n_centers = static_cast<int>(centers_x.size());
    if (n_centers != static_cast<int>(centers_y.size()) ||
        n_centers != static_cast<int>(weights.size()) || n_centers == 0) {
        std::cerr << "rbf_friction_cuda: invalid RBF sizes" << std::endl;
        return false;
    }
    const int nx = params.nx;
    const int ny = params.ny;
    const int n = nx * ny;
    const double sigma2 = sigma * sigma;
    const double mu = params.mu;
    const double nz_min = 1.0 / std::sqrt(1.0 + mu * mu);

    double* d_cx = nullptr;
    double* d_cy = nullptr;
    double* d_w = nullptr;
    int* d_labels = nullptr;

    cudaError_t err;
    err = cudaMalloc(&d_cx, n_centers * sizeof(double));
    if (err != cudaSuccess) { std::cerr << "cudaMalloc d_cx: " << cudaGetErrorString(err) << std::endl; return false; }
    err = cudaMalloc(&d_cy, n_centers * sizeof(double));
    if (err != cudaSuccess) { cudaFree(d_cx); std::cerr << "cudaMalloc d_cy: " << cudaGetErrorString(err) << std::endl; return false; }
    err = cudaMalloc(&d_w, n_centers * sizeof(double));
    if (err != cudaSuccess) { cudaFree(d_cx); cudaFree(d_cy); std::cerr << "cudaMalloc d_w: " << cudaGetErrorString(err) << std::endl; return false; }
    err = cudaMalloc(&d_labels, n * sizeof(int));
    if (err != cudaSuccess) {
        cudaFree(d_cx); cudaFree(d_cy); cudaFree(d_w);
        std::cerr << "cudaMalloc d_labels: " << cudaGetErrorString(err) << std::endl;
        return false;
    }

    err = cudaMemcpy(d_cx, centers_x.data(), n_centers * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) goto cleanup;
    err = cudaMemcpy(d_cy, centers_y.data(), n_centers * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) goto cleanup;
    err = cudaMemcpy(d_w, weights.data(), n_centers * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) goto cleanup;

    launch_friction_cone(d_cx, d_cy, d_w, n_centers, sigma2,
                         params.xmin, params.xmax, params.ymin, params.ymax,
                         nx, ny, nz_min, d_labels);

    err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        std::cerr << "kernel: " << cudaGetErrorString(err) << std::endl;
        goto cleanup;
    }

    out_labels.resize(n);
    err = cudaMemcpy(out_labels.data(), d_labels, n * sizeof(int), cudaMemcpyDeviceToHost);

cleanup:
    cudaFree(d_cx);
    cudaFree(d_cy);
    cudaFree(d_w);
    cudaFree(d_labels);

    if (err != cudaSuccess) {
        std::cerr << "cudaMemcpy D2H: " << cudaGetErrorString(err) << std::endl;
        return false;
    }
    return true;
}

bool saveSegmentation(const std::string& path,
                      int nx, int ny,
                      double xmin, double xmax, double ymin, double ymax,
                      const std::vector<int>& labels) {
    if (static_cast<size_t>(nx * ny) != labels.size()) return false;
    std::ofstream out(path);
    if (!out) return false;
    out.precision(16);
    out << nx << " " << ny << " " << xmin << " " << xmax << " " << ymin << " " << ymax << "\n";
    for (int v : labels) out << v << " ";
    out << "\n";
    return !!out;
}

bool saveSegmentationImage(const std::string& path,
                          int nx, int ny,
                          const std::vector<int>& labels,
                          bool open_after_save) {
    if (static_cast<size_t>(nx * ny) != labels.size()) return false;
    std::vector<unsigned char> rgb(static_cast<size_t>(nx) * ny * 3);
    const unsigned char red[] = { 0xc0, 0x39, 0x2b };
    const unsigned char green[] = { 0x27, 0xae, 0x60 };
    for (int j = 0; j < ny; ++j) {
        int img_row = ny - 1 - j;
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;
            const unsigned char* c = labels[idx] ? green : red;
            size_t out_idx = (static_cast<size_t>(img_row) * nx + i) * 3;
            rgb[out_idx] = c[0];
            rgb[out_idx + 1] = c[1];
            rgb[out_idx + 2] = c[2];
        }
    }
    if (!stbi_write_png(path.c_str(), nx, ny, 3, rgb.data(), 0)) {
        std::cerr << "Failed to write PNG: " << path << std::endl;
        return false;
    }
    if (open_after_save) {
        std::string cmd = "xdg-open \"" + path + "\" 2>/dev/null || open \"" + path + "\" 2>/dev/null || true";
        (void)std::system(cmd.c_str());
    }
    return true;
}

static void heightToRgb(double t, unsigned char& r, unsigned char& g, unsigned char& b) {
    if (t <= 0.0) { r = 30; g = 60; b = 180; return; }
    if (t >= 1.0) { r = 200; g = 50; b = 50; return; }
    if (t < 0.5) {
        double s = t * 2.0;
        r = static_cast<unsigned char>(30 + s * (100 - 30));
        g = static_cast<unsigned char>(60 + s * (180 - 60));
        b = static_cast<unsigned char>(180 + s * (60 - 180));
    } else {
        double s = (t - 0.5) * 2.0;
        r = static_cast<unsigned char>(100 + s * (200 - 100));
        g = static_cast<unsigned char>(180 + s * (50 - 180));
        b = static_cast<unsigned char>(60 + s * (50 - 60));
    }
}

bool saveTerrainSegmentationComparison(
    const std::string& path,
    int nx, int ny,
    const std::vector<double>& terrain_height,
    const std::vector<int>& labels,
    bool open_after_save) {
    const size_t n = static_cast<size_t>(nx) * ny;
    if (terrain_height.size() != n || labels.size() != n) return false;
    const int gap = 12;
    const int total_w = 2 * nx + gap;
    std::vector<unsigned char> rgb(static_cast<size_t>(total_w) * ny * 3, 255);
    double hmin = terrain_height[0], hmax = terrain_height[0];
    for (double h : terrain_height) {
        if (h < hmin) hmin = h;
        if (h > hmax) hmax = h;
    }
    double inv_span = (hmax > hmin) ? (1.0 / (hmax - hmin)) : 1.0;
    const unsigned char red[] = { 0xc0, 0x39, 0x2b };
    const unsigned char green[] = { 0x27, 0xae, 0x60 };
    for (int j = 0; j < ny; ++j) {
        int img_row = ny - 1 - j;
        for (int i = 0; i < nx; ++i) {
            size_t idx = static_cast<size_t>(j) * nx + i;
            double t = (terrain_height[idx] - hmin) * inv_span;
            unsigned char r, g, b;
            heightToRgb(t, r, g, b);
            size_t out_left = (static_cast<size_t>(img_row) * total_w + i) * 3;
            rgb[out_left] = r;
            rgb[out_left + 1] = g;
            rgb[out_left + 2] = b;
            const unsigned char* c = labels[idx] ? green : red;
            size_t out_right = (static_cast<size_t>(img_row) * total_w + gap + nx + i) * 3;
            rgb[out_right] = c[0];
            rgb[out_right + 1] = c[1];
            rgb[out_right + 2] = c[2];
        }
    }
    if (!stbi_write_png(path.c_str(), total_w, ny, 3, rgb.data(), 0)) {
        std::cerr << "Failed to write comparison PNG: " << path << std::endl;
        return false;
    }
    if (open_after_save) {
        std::string cmd = "xdg-open \"" + path + "\" 2>/dev/null || open \"" + path + "\" 2>/dev/null || true";
        (void)std::system(cmd.c_str());
    }
    return true;
}

} // namespace rbf_friction_cuda
