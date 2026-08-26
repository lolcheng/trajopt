#include "rrt_star_region.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace rrt_star_region {

namespace {

inline double sq(double x) { return x * x; }
inline double dist_sq(double x1, double y1, double x2, double y2) {
    return sq(x2 - x1) + sq(y2 - y1);
}
inline double dist(double x1, double y1, double x2, double y2) {
    return std::sqrt(dist_sq(x1, y1, x2, y2));
}

struct SegData {
    int nx = 0, ny = 0;
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    std::vector<int> labels;  // row-major labels[j*nx+i]

    int get_label(double x, double y) const {
        double xi = (nx > 1) ? (x - xmin) / (xmax - xmin) * (nx - 1) : 0;
        double yj = (ny > 1) ? (y - ymin) / (ymax - ymin) * (ny - 1) : 0;
        int i = static_cast<int>(std::round(std::max(0.0, std::min(static_cast<double>(nx - 1), xi))));
        int j = static_cast<int>(std::round(std::max(0.0, std::min(static_cast<double>(ny - 1), yj))));
        return labels[j * nx + i];
    }

    double edge_cost(double x1, double y1, double x2, double y2, double k, int n_samples) const {
        const double dx = (x2 - x1) / n_samples;
        const double dy = (y2 - y1) / n_samples;
        const double seg_len = dist(x1, y1, x2, y2) / n_samples;
        double cost = 0;
        for (int i = 0; i < n_samples; ++i) {
            double x = x1 + (i + 0.5) * dx;
            double y = y1 + (i + 0.5) * dy;
            int L = get_label(x, y);
            cost += seg_len * (L == 1 ? 1.0 : k);
        }
        return cost;
    }
};

bool load_segmentation(const std::string& path, SegData& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;
    std::istringstream ss(line);
    if (!(ss >> out.nx >> out.ny >> out.xmin >> out.xmax >> out.ymin >> out.ymax)) return false;
    if (!std::getline(f, line)) return false;  // skip terrain
    if (!std::getline(f, line)) return false;
    out.labels.clear();
    std::istringstream ss2(line);
    int v;
    while (ss2 >> v) out.labels.push_back(v);
    if (static_cast<int>(out.labels.size()) != out.nx * out.ny) return false;
    return true;
}

// Spatial grid: bucket index (bx, by) -> list of node indices
struct NodeGrid {
    int nbx, nby;
    double xmin, xmax, ymin, ymax;
    std::vector<std::vector<int>> buckets;  // buckets[bx + by*nbx]

    NodeGrid(int n, double xmin_, double xmax_, double ymin_, double ymax_)
        : nbx(n), nby(n), xmin(xmin_), xmax(xmax_), ymin(ymin_), ymax(ymax_) {
        buckets.resize(nbx * nby);
    }

    int bucket(double x, double y) const {
        int bx = static_cast<int>((x - xmin) / (xmax - xmin + 1e-10) * nbx);
        int by = static_cast<int>((y - ymin) / (ymax - ymin + 1e-10) * nby);
        bx = std::max(0, std::min(nbx - 1, bx));
        by = std::max(0, std::min(nby - 1, by));
        return bx + by * nbx;
    }

    void add(int node_idx, double x, double y) {
        int b = bucket(x, y);
        buckets[b].push_back(node_idx);
    }

    // Nearest node index by squared distance (avoid sqrt)
    int nearest(const std::vector<double>& xs, const std::vector<double>& ys,
                double x, double y, double max_radius_sq) const {
        int best = -1;
        double best_sq = max_radius_sq;
        const int bx0 = static_cast<int>((x - xmin) / (xmax - xmin + 1e-10) * nbx);
        const int by0 = static_cast<int>((y - ymin) / (ymax - ymin + 1e-10) * nby);
        const int rad = std::max(2, std::min(nbx, nby) / 4);
        for (int dy = -rad; dy <= rad; ++dy) {
            for (int dx = -rad; dx <= rad; ++dx) {
                int bx = bx0 + dx, by = by0 + dy;
                if (bx < 0 || bx >= nbx || by < 0 || by >= nby) continue;
                const auto& list = buckets[bx + by * nbx];
                for (int idx : list) {
                    double d2 = dist_sq(xs[idx], ys[idx], x, y);
                    if (d2 < best_sq) { best_sq = d2; best = idx; }
                }
            }
        }
        return best;
    }

    void nodes_in_radius(const std::vector<double>& xs, const std::vector<double>& ys,
                         double x, double y, double r, std::vector<int>& out) const {
        out.clear();
        const double r2 = r * r;
        const int bx0 = static_cast<int>((x - xmin) / (xmax - xmin + 1e-10) * nbx);
        const int by0 = static_cast<int>((y - ymin) / (ymax - ymin + 1e-10) * nby);
        const int cell_rad = std::max(1, static_cast<int>(std::ceil(r / ((xmax - xmin) / nbx))));
        for (int dy = -cell_rad; dy <= cell_rad; ++dy) {
            for (int dx = -cell_rad; dx <= cell_rad; ++dx) {
                int bx = bx0 + dx, by = by0 + dy;
                if (bx < 0 || bx >= nbx || by < 0 || by >= nby) continue;
                for (int idx : buckets[bx + by * nbx])
                    if (dist_sq(xs[idx], ys[idx], x, y) <= r2) out.push_back(idx);
            }
        }
    }
};

void steer(double from_x, double from_y, double to_x, double to_y, double step,
          double& out_x, double& out_y) {
    double d = dist(from_x, from_y, to_x, to_y);
    if (d <= step) { out_x = to_x; out_y = to_y; return; }
    double t = step / d;
    out_x = from_x + t * (to_x - from_x);
    out_y = from_y + t * (to_y - from_y);
}

}  // namespace

std::vector<std::pair<double, double>> run(
    const std::string& data_file,
    const Params& params,
    double* out_time_ms) {

    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<double, double>> path;

    SegData seg;
    if (!load_segmentation(data_file, seg)) return path;

    const double xmin = seg.xmin, xmax = seg.xmax, ymin = seg.ymin, ymax = seg.ymax;
    const double k = params.lift_cost_weight;
    const double step = params.step_length;
    const double goal_r = params.goal_radius;
    const double rewire_r = params.rewire_radius;
    const double rewire_r_sq = rewire_r * rewire_r;
    const int edge_samples = params.edge_samples;
    const int n_buckets = params.grid_buckets;

    std::mt19937 rng(0);
    std::uniform_real_distribution<double> u01(0, 1);
    std::uniform_real_distribution<double> ux(xmin, xmax);
    std::uniform_real_distribution<double> uy(ymin, ymax);

    std::vector<double> xs, ys;
    std::vector<int> parent;
    std::vector<double> cost_from_root;
    xs.push_back(params.start_x);
    ys.push_back(params.start_y);
    parent.push_back(-1);
    cost_from_root.push_back(0.0);

    NodeGrid grid(n_buckets, xmin, xmax, ymin, ymax);
    grid.add(0, params.start_x, params.start_y);

    const double max_nearest_sq = (step * 3) * (step * 3);
    int goal_idx = -1;

    for (int iter = 0; iter < params.max_iter; ++iter) {
        double rnd_x, rnd_y;
        if (u01(rng) < 0.1) {
            rnd_x = params.goal_x;
            rnd_y = params.goal_y;
        } else {
            rnd_x = ux(rng);
            rnd_y = uy(rng);
        }

        int nearest_idx = grid.nearest(xs, ys, rnd_x, rnd_y, max_nearest_sq);
        if (nearest_idx < 0) {
            double best_sq = 1e30;
            for (size_t i = 0; i < xs.size(); ++i) {
                double d2 = dist_sq(xs[i], ys[i], rnd_x, rnd_y);
                if (d2 < best_sq) { best_sq = d2; nearest_idx = static_cast<int>(i); }
            }
            if (nearest_idx < 0) nearest_idx = 0;
        }

        double new_x, new_y;
        steer(xs[nearest_idx], ys[nearest_idx], rnd_x, rnd_y, step, new_x, new_y);

        double c_edge = seg.edge_cost(xs[nearest_idx], ys[nearest_idx], new_x, new_y, k, edge_samples);
        double new_cost = cost_from_root[nearest_idx] + c_edge;

        if (new_x < xmin || new_x > xmax || new_y < ymin || new_y > ymax) continue;

        std::vector<int> in_radius;
        grid.nodes_in_radius(xs, ys, new_x, new_y, rewire_r, in_radius);
        for (int i : in_radius) {
            double d2 = dist_sq(xs[i], ys[i], new_x, new_y);
            if (d2 <= rewire_r_sq) {
                double c = seg.edge_cost(xs[i], ys[i], new_x, new_y, k, edge_samples);
                if (cost_from_root[i] + c < new_cost) {
                    new_cost = cost_from_root[i] + c;
                    nearest_idx = i;
                }
            }
        }

        int new_idx = static_cast<int>(xs.size());
        xs.push_back(new_x);
        ys.push_back(new_y);
        parent.push_back(nearest_idx);
        cost_from_root.push_back(new_cost);
        grid.add(new_idx, new_x, new_y);

        for (int i : in_radius) {
            if (i == new_idx) continue;
            double d2 = dist_sq(xs[i], ys[i], new_x, new_y);
            if (d2 <= rewire_r_sq) {
                double c = seg.edge_cost(new_x, new_y, xs[i], ys[i], k, edge_samples);
                if (new_cost + c < cost_from_root[i]) {
                    parent[i] = new_idx;
                    cost_from_root[i] = new_cost + c;
                }
            }
        }

        if (dist_sq(new_x, new_y, params.goal_x, params.goal_y) <= goal_r * goal_r) {
            double c_end = seg.edge_cost(new_x, new_y, params.goal_x, params.goal_y, k, edge_samples);
            double total = new_cost + c_end;
            if (goal_idx < 0 || total < cost_from_root[goal_idx]) {
                goal_idx = static_cast<int>(xs.size());
                xs.push_back(params.goal_x);
                ys.push_back(params.goal_y);
                parent.push_back(new_idx);
                cost_from_root.push_back(total);
                grid.add(goal_idx, params.goal_x, params.goal_y);

                path.clear();
                int idx = goal_idx;
                while (idx >= 0) {
                    path.emplace_back(xs[idx], ys[idx]);
                    idx = parent[idx];
                }
                std::reverse(path.begin(), path.end());
            }
        }
    }

    if (out_time_ms) {
        auto t1 = std::chrono::high_resolution_clock::now();
        *out_time_ms = 1e-6 * std::chrono::duration<double, std::nano>(t1 - t0).count();
    }
    return path;
}

// --- Public API: SegDataView, load_segmentation, print_switch_waypoints, write_visualization_png ---

int SegDataView::get_label(double x, double y) const {
    if (nx == 0 || ny == 0) return 0;
    double xi = (nx > 1) ? (x - xmin) / (xmax - xmin + 1e-10) * (nx - 1) : 0;
    double yj = (ny > 1) ? (y - ymin) / (ymax - ymin + 1e-10) * (ny - 1) : 0;
    int i = static_cast<int>(std::round(std::max(0.0, std::min(static_cast<double>(nx - 1), xi))));
    int j = static_cast<int>(std::round(std::max(0.0, std::min(static_cast<double>(ny - 1), yj))));
    return labels[j * nx + i];
}

bool load_segmentation(const std::string& path, SegDataView& out) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;
    std::istringstream ss(line);
    if (!(ss >> out.nx >> out.ny >> out.xmin >> out.xmax >> out.ymin >> out.ymax)) return false;
    if (!std::getline(f, line)) return false;
    if (!std::getline(f, line)) return false;
    out.labels.clear();
    std::istringstream ss2(line);
    int v;
    while (ss2 >> v) out.labels.push_back(v);
    return static_cast<int>(out.labels.size()) == out.nx * out.ny;
}

void print_switch_waypoints(const std::vector<std::pair<double, double>>& path,
                            const SegDataView& seg) {
    if (path.size() < 2) return;
    int prev = seg.get_label(path[0].first, path[0].second);
    for (size_t i = 1; i < path.size(); ++i) {
        int cur = seg.get_label(path[i].first, path[i].second);
        if (cur != prev) {
            std::cout << "  switch waypoint [" << i << "] (" << path[i].first << ", " << path[i].second
                      << ") " << (prev == 1 ? "rollable" : "lift") << " -> "
                      << (cur == 1 ? "rollable" : "lift") << std::endl;
            prev = cur;
        }
    }
}

namespace {
void draw_line(std::vector<unsigned char>& rgb, int W, int H,
               double xmin, double xmax, double ymin, double ymax,
               double x1, double y1, double x2, double y2,
               unsigned char r, unsigned char g, unsigned char b) {
    auto to_px = [&](double x, double y) -> std::pair<int, int> {
        int px = static_cast<int>(std::round((x - xmin) / (xmax - xmin + 1e-10) * (W - 1)));
        int py = static_cast<int>(std::round((y - ymin) / (ymax - ymin + 1e-10) * (H - 1)));
        py = H - 1 - py;
        return {std::max(0, std::min(W - 1, px)), std::max(0, std::min(H - 1, py))};
    };
    auto [p0x, p0y] = to_px(x1, y1);
    auto [p1x, p1y] = to_px(x2, y2);
    int dx = std::abs(p1x - p0x), sx = p0x < p1x ? 1 : -1;
    int dy = -std::abs(p1y - p0y), sy = p0y < p1y ? 1 : -1;
    int err = dx + dy;
    const int thick = 2;
    for (;;) {
        for (int dy2 = -thick; dy2 <= thick; ++dy2)
            for (int dx2 = -thick; dx2 <= thick; ++dx2) {
                int qx = p0x + dx2, qy = p0y + dy2;
                if (qx >= 0 && qx < W && qy >= 0 && qy < H) {
                    size_t idx = (static_cast<size_t>(qy) * W + qx) * 3;
                    rgb[idx] = r; rgb[idx + 1] = g; rgb[idx + 2] = b;
                }
            }
        if (p0x == p1x && p0y == p1y) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; p0x += sx; }
        if (e2 <= dx) { err += dx; p0y += sy; }
    }
}
void draw_circle(std::vector<unsigned char>& rgb, int W, int H,
                 double xmin, double xmax, double ymin, double ymax,
                 double cx, double cy, int radius,
                 unsigned char r, unsigned char g, unsigned char b) {
    int px = static_cast<int>(std::round((cx - xmin) / (xmax - xmin + 1e-10) * (W - 1)));
    int py = static_cast<int>(std::round((cy - ymin) / (ymax - ymin + 1e-10) * (H - 1)));
    py = H - 1 - py;
    for (int dy = -radius; dy <= radius; ++dy)
        for (int dx = -radius; dx <= radius; ++dx)
            if (dx * dx + dy * dy <= radius * radius) {
                int qx = px + dx, qy = py + dy;
                if (qx >= 0 && qx < W && qy >= 0 && qy < H) {
                    size_t idx = (static_cast<size_t>(qy) * W + qx) * 3;
                    rgb[idx] = r; rgb[idx + 1] = g; rgb[idx + 2] = b;
                }
            }
}
}  // namespace

bool write_visualization_png(const std::string& out_path,
                             const SegDataView& seg,
                             const std::vector<std::pair<double, double>>& path,
                             double start_x, double start_y, double goal_x, double goal_y,
                             bool open_after) {
    const int W = seg.nx * 2, H = seg.ny * 2;
    if (W <= 0 || H <= 0) return false;
    std::vector<unsigned char> rgb(static_cast<size_t>(W) * H * 3);
    const unsigned char red[] = {0xc0, 0x39, 0x2b}, green[] = {0x27, 0xae, 0x60};
    for (int py = 0; py < H; ++py) {
        double y = seg.ymax - (py + 0.5) / H * (seg.ymax - seg.ymin);
        for (int px = 0; px < W; ++px) {
            double x = seg.xmin + (px + 0.5) / W * (seg.xmax - seg.xmin);
            int L = seg.get_label(x, y);
            size_t idx = (static_cast<size_t>(py) * W + px) * 3;
            rgb[idx] = L == 1 ? green[0] : red[0];
            rgb[idx + 1] = L == 1 ? green[1] : red[1];
            rgb[idx + 2] = L == 1 ? green[2] : red[2];
        }
    }
    for (size_t i = 0; i + 1 < path.size(); ++i)
        draw_line(rgb, W, H, seg.xmin, seg.xmax, seg.ymin, seg.ymax,
                  path[i].first, path[i].second, path[i + 1].first, path[i + 1].second,
                  0, 0, 255);
    draw_circle(rgb, W, H, seg.xmin, seg.xmax, seg.ymin, seg.ymax, start_x, start_y, 8, 0, 255, 0);
    draw_circle(rgb, W, H, seg.xmin, seg.xmax, seg.ymin, seg.ymax, goal_x, goal_y, 10, 255, 165, 0);
    if (!stbi_write_png(out_path.c_str(), W, H, 3, rgb.data(), 0)) return false;
    if (open_after) {
        std::string cmd = "xdg-open \"" + out_path + "\" 2>/dev/null || open \"" + out_path + "\" 2>/dev/null || true";
        (void)std::system(cmd.c_str());
    }
    return true;
}

bool save_waypoints_segmented(const std::string& out_path,
                              const std::vector<std::pair<double, double>>& path,
                              const SegDataView& seg) {
    std::ofstream f(out_path);
    if (!f) return false;
    f << std::fixed << std::setprecision(6);
    for (const auto& p : path) {
        int label = seg.get_label(p.first, p.second);
        f << p.first << " " << p.second << " " << label << "\n";
    }
    return !!f;
}

}  // namespace rrt_star_region
