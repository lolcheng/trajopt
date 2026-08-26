#include <stddef.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "IpIpoptApplication.hpp"
#include "basis_function.hpp"
#include "chebyshev_basis.hpp"
#include "constraints/boundary_constraint_10var.hpp"
#include "constraints/coefficient_norm_constraint.hpp"
#include "constraints/penetration_based_force_constraint.hpp"
#include "constraints/piecewise_z_continuity_constraint.hpp"
#include "constraints/trajectory_bounds_constraint_8var.hpp"
#include "constraints/wheel_base_distance_terrain_constraint.hpp"
#include "constraints/wheel_terrain_penetration_constraint.hpp"
#include "constraints/wheels_distance_terrain_constraint.hpp"
#include "json_reader.hpp"
#include "objectives/base_height_stability_objective.hpp"
#include "objectives/initial_guess_proximity_objective.hpp"
#include "objectives/swing_wheel_lift_objective.hpp"
#include "objectives/yaw_path_alignment_objective.hpp"
#include "polynomial.hpp"
#include "rbf_terrain.hpp"
#include "seg_trajopt_nlp.hpp"
#include "waypoint_reader.hpp"

static double liftProfile(double s, double peak) {
    const double ss = std::sin(M_PI * s);
    return peak * ss * ss;
}

static bool projectToFrictionSafeXY(double& x, double& y, const RBFTerrain& terrain,
                                    double mu, double search_r = 0.25, int grid_n = 7) {
    double hx0, hy0;
    terrain.gradient(x, y, hx0, hy0);
    if (std::hypot(hx0, hy0) <= mu) return true;
    double best_x = x, best_y = y, best_s = 1e18;
    const double step = (2.0 * search_r) / std::max(1, grid_n - 1);
    for (int ix = 0; ix < grid_n; ++ix) {
        for (int iy = 0; iy < grid_n; ++iy) {
            const double tx = x - search_r + step * ix;
            const double ty = y - search_r + step * iy;
            double hx, hy;
            terrain.gradient(tx, ty, hx, hy);
            const double s = std::hypot(hx, hy);
            if (s < best_s) {
                best_s = s;
                best_x = tx;
                best_y = ty;
            }
        }
    }
    x = best_x;
    y = best_y;
    return best_s <= mu;
}

static void fitCoefficientsFromWaypoints(
    int n_coeff, int n_seg_z, const std::vector<double>& s,
    const std::vector<double>& xb, const std::vector<double>& yb,
    const std::vector<double>& zb, const std::vector<double>& psi,
    const std::vector<double>& xl, const std::vector<double>& yl, const std::vector<double>& zl,
    const std::vector<double>& xr, const std::vector<double>& yr, const std::vector<double>& zr,
    const std::shared_ptr<BasisFunction>& basis,
    std::vector<double>& cbx, std::vector<double>& cby, std::vector<double>& cbz, std::vector<double>& cpsi,
    std::vector<double>& clx, std::vector<double>& cly, std::vector<double>& clz,
    std::vector<double>& crx, std::vector<double>& cry, std::vector<double>& crz) {
    const int m = static_cast<int>(s.size());
    const int n = n_coeff;
    const int nz = n_coeff * std::max(1, n_seg_z);
    const double eps = 1e-12;
    auto lsq = [&](const std::vector<double>& s_fit,
                   const std::vector<double>& y,
                   std::vector<double>& coeff) {
        const int m_fit = static_cast<int>(s_fit.size());
        std::vector<double> Phi(m_fit * n);
        basis->computeBasisMatrix(s_fit.data(), m_fit, n, Phi.data());
        std::vector<double> Q(m_fit * n, 0.0), R(n * n, 0.0);
        for (int j = 0; j < n; ++j) {
            std::vector<double> v(m_fit, 0.0);
            for (int i = 0; i < m_fit; ++i) v[i] = Phi[i * n + j];
            for (int k = 0; k < j; ++k) {
                double rkj = 0.0;
                for (int i = 0; i < m_fit; ++i) rkj += Q[i * n + k] * v[i];
                R[k * n + j] = rkj;
                for (int i = 0; i < m_fit; ++i) v[i] -= rkj * Q[i * n + k];
            }
            double nv = 0.0;
            for (int i = 0; i < m_fit; ++i) nv += v[i] * v[i];
            nv = std::sqrt(nv);
            R[j * n + j] = nv;
            if (nv > eps) for (int i = 0; i < m_fit; ++i) Q[i * n + j] = v[i] / nv;
        }
        std::vector<double> qTy(n, 0.0);
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m_fit; ++i) qTy[j] += Q[i * n + j] * y[i];
        }
        coeff.assign(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            coeff[i] = qTy[i];
            for (int j = i + 1; j < n; ++j) coeff[i] -= R[i * n + j] * coeff[j];
            coeff[i] = (std::abs(R[i * n + i]) > eps) ? coeff[i] / R[i * n + i] : 0.0;
        }
    };

    lsq(s, xb, cbx); lsq(s, yb, cby); lsq(s, zb, cbz); lsq(s, psi, cpsi);
    lsq(s, xl, clx); lsq(s, yl, cly);
    lsq(s, xr, crx); lsq(s, yr, cry);

    clz.assign(nz, 0.0);
    crz.assign(nz, 0.0);
    const int seg_count = std::max(1, n_seg_z);
    for (int seg = 0; seg < seg_count; ++seg) {
        std::vector<double> s_local, zl_local, zr_local;
        const double lo = static_cast<double>(seg) / seg_count;
        const double hi = static_cast<double>(seg + 1) / seg_count;
        for (int i = 0; i < m; ++i) {
            const bool in_seg = (seg == seg_count - 1) ? (s[i] >= lo && s[i] <= hi) : (s[i] >= lo && s[i] < hi);
            if (!in_seg) continue;
            const double sl = std::clamp(s[i] * seg_count - seg, 0.0, 1.0);
            s_local.push_back(sl);
            zl_local.push_back(zl[i]);
            zr_local.push_back(zr[i]);
        }
        if (s_local.empty()) continue;
        std::vector<double> clz_seg, crz_seg;
        lsq(s_local, zl_local, clz_seg);
        lsq(s_local, zr_local, crz_seg);
        std::copy(clz_seg.begin(), clz_seg.end(), clz.begin() + seg * n_coeff);
        std::copy(crz_seg.begin(), crz_seg.end(), crz.begin() + seg * n_coeff);
    }
}

static bool writeTrajectoryFile(const std::string& path,
                                const std::vector<double>& Xb, const std::vector<double>& Yb,
                                const std::vector<double>& Zb, const std::vector<double>& Psib,
                                const std::vector<double>& Xl, const std::vector<double>& Yl, const std::vector<double>& Zl,
                                const std::vector<double>& Xr, const std::vector<double>& Yr, const std::vector<double>& Zr) {
    const int n = static_cast<int>(Xb.size());
    std::ofstream out(path);
    if (!out || n == 0) return false;
    out << std::fixed << std::setprecision(10);
    out << "s xb yb zb psi xl yl zl xr yr zr\n";
    for (int i = 0; i < n; ++i) {
        const double s = (n <= 1) ? 0.5 : static_cast<double>(i) / (n - 1);
        out << s << " " << Xb[i] << " " << Yb[i] << " " << Zb[i] << " " << Psib[i] << " "
            << Xl[i] << " " << Yl[i] << " " << Zl[i] << " "
            << Xr[i] << " " << Yr[i] << " " << Zr[i] << "\n";
    }
    return true;
}

int main(int argc, char** argv) {
    const int N_COEFF = 8;
    const int N_SAMPLE = 21;
    const double CLEAR_Z = 0.8;
    const double WHEELS_NOMINAL = 0.5;
    const double WHEEL_BASE_MIN = 0.2, WHEEL_BASE_MAX = 1.2;
    const double WHEELS_MIN = 0.31622776601683794, WHEELS_MAX = 1.0;
    const double TERRAIN_MARGIN = 0.2;
    const double COEFF_NORM_MAX = 10000.0;
    const double MAX_PENETRATION = 0.01;
    const double CONTACT_STIFFNESS = 5000.0;
    const double STANCE_FORCE_MIN = 5.0;
    const double STANCE_FORCE_MAX = 700.0;
    const double SWING_FORCE_MAX = 8.0;
    const double STANCE_CLEARANCE = 0.005;
    const double LIFT_PEAK = 0.24;
    const double SWING_NO_CONTACT_MARGIN = 0.03;
    const double FRICTION_MU = 0.55;
    const int INIT_DENSE_FACTOR = 12;
    const double STAIR_HEIGHT_EST = 0.20;
    const double BOUNDARY_BASE_XY_TOL = 1e-3;
    const double BOUNDARY_BASE_Z_TOL = 1e-3;
    const double BOUNDARY_BASE_YAW_TOL = 1e-3;
    const double BOUNDARY_WHEEL_XY_TOL = 5e-2;
    const double BOUNDARY_WHEEL_Z_TOL = 2e-2;
    const double BOUNDARY_BASE_Z_MIN = -2.0;
    const double BOUNDARY_BASE_Z_MAX = 3.0;
    const double BOUNDARY_WHEEL_Z_MIN = -2.0;
    const double BOUNDARY_WHEEL_Z_MAX = 3.0;
    const double PIECEWISE_Z_C0_TOL = 5e-3;
    const double INITIAL_GUESS_PROX_WEIGHT = 2e-2;
    const double BASE_HEIGHT_STABILITY_WEIGHT = 2.0 / static_cast<double>(N_SAMPLE);
    const double YAW_PATH_ALIGNMENT_WEIGHT = 8.0 / static_cast<double>(N_SAMPLE);
    const double SWING_LIFT_OBJECTIVE_WEIGHT = 6.0 / static_cast<double>(N_SAMPLE);

    auto t_total_start = std::chrono::high_resolution_clock::now();
    std::string data_dir = (argc > 1) ? argv[1] : std::string("terrain_res");
    if (data_dir.back() != '/') data_dir += '/';
    std::string out_dir = (argc > 2) ? argv[2] : std::string("lift_out");
    if (out_dir.back() != '/') out_dir += '/';
    std::string dir_for_mkdir = out_dir;
    if (!dir_for_mkdir.empty() && dir_for_mkdir.back() == '/') dir_for_mkdir.pop_back();
    if (!dir_for_mkdir.empty()) {
#ifdef _WIN32
        _mkdir(dir_for_mkdir.c_str());
#else
        mkdir(dir_for_mkdir.c_str(), 0755);
#endif
    }

    std::vector<double> centers_x, centers_y, weights;
    double sigma;
    const std::string rbf_file = data_dir + "rbf.json";
    if (!JSONReader::readRBFParams(rbf_file, centers_x, centers_y, weights, sigma)) {
        std::cerr << "Failed to read RBF from " << rbf_file << "\n";
        return -1;
    }
    auto terrain = std::make_shared<RBFTerrain>(centers_x, centers_y, weights, sigma);

    std::vector<double> wp_x, wp_y;
    const std::string waypoints_file = data_dir + "waypoints_segmented.txt";
    if (!WaypointReader::readFirstLiftLegSegment(waypoints_file, wp_x, wp_y)) {
        std::cerr << "Warning: no lift-leg segment found, fallback to first rollable segment.\n";
        if (!WaypointReader::readFirstRollableSegment(waypoints_file, wp_x, wp_y)) {
            std::cerr << "Failed to read segment from " << waypoints_file << "\n";
            return -1;
        }
    }
    const int n_wp_raw = static_cast<int>(wp_x.size());
    if (n_wp_raw < 2) {
        std::cerr << "Need at least 2 waypoints.\n";
        return -1;
    }

    const double x0 = wp_x.front(), y0 = wp_y.front();
    const double x1 = wp_x.back(), y1 = wp_y.back();
    const int n_wp = std::max(n_wp_raw * INIT_DENSE_FACTOR, 25);
    std::vector<double> wp_x_dense(n_wp, 0.0), wp_y_dense(n_wp, 0.0);
    for (int i = 0; i < n_wp; ++i) {
        const double u = (n_wp == 1) ? 0.0 : static_cast<double>(i) * (n_wp_raw - 1) / (n_wp - 1);
        const int k0 = static_cast<int>(std::floor(u));
        const int k1 = std::min(k0 + 1, n_wp_raw - 1);
        const double a = u - k0;
        wp_x_dense[i] = (1.0 - a) * wp_x[k0] + a * wp_x[k1];
        wp_y_dense[i] = (1.0 - a) * wp_y[k0] + a * wp_y[k1];
    }
    std::shared_ptr<BasisFunction> basis = std::make_shared<ChebyshevBasis>();
    Polynomial::set_basis(basis);

    std::vector<double> s_wp(n_wp), xb_wp(n_wp), yb_wp(n_wp), zb_wp(n_wp), psi_wp(n_wp);
    std::vector<double> xl_wp(n_wp), yl_wp(n_wp), zl_wp(n_wp);
    std::vector<double> xr_wp(n_wp), yr_wp(n_wp), zr_wp(n_wp);
    for (int i = 0; i < n_wp; ++i) {
        s_wp[i] = (n_wp <= 1) ? 0.5 : static_cast<double>(i) / (n_wp - 1);
        xb_wp[i] = wp_x_dense[i];
        yb_wp[i] = wp_y_dense[i];
        zb_wp[i] = terrain->height(wp_x_dense[i], wp_y_dense[i]) + CLEAR_Z;
        const double dx = (i < n_wp - 1) ? (wp_x_dense[i + 1] - wp_x_dense[i]) : (wp_x_dense[i] - wp_x_dense[i - 1]);
        const double dy = (i < n_wp - 1) ? (wp_y_dense[i + 1] - wp_y_dense[i]) : (wp_y_dense[i] - wp_y_dense[i - 1]);
        psi_wp[i] = std::atan2(dy, dx);
        const double ux = -std::sin(psi_wp[i]);
        const double uy = std::cos(psi_wp[i]);
        xl_wp[i] = xb_wp[i] + 0.5 * WHEELS_NOMINAL * ux;
        yl_wp[i] = yb_wp[i] + 0.5 * WHEELS_NOMINAL * uy;
        xr_wp[i] = xb_wp[i] - 0.5 * WHEELS_NOMINAL * ux;
        yr_wp[i] = yb_wp[i] - 0.5 * WHEELS_NOMINAL * uy;
    }

    const double h_start = terrain->height(x0, y0);
    const double h_end = terrain->height(x1, y1);
    const double dz_base = h_end - h_start;
    const int stair_cycles = std::max(0, static_cast<int>(std::ceil(std::max(0.0, dz_base) / STAIR_HEIGHT_EST)));
    const int N_SEG_Z = std::max(1, stair_cycles);
    std::vector<int> left_contact_wp(n_wp, 1), right_contact_wp(n_wp, 1);
    std::vector<double> left_swing_z(n_wp, -1e20), right_swing_z(n_wp, -1e20);
    std::vector<int> touchdown_left, touchdown_right;
    const int swing_half_width = std::max(1, n_wp / 14);
    const int swing_shift = swing_half_width + std::max(1, n_wp / 20);
    auto applySwingWindow = [&](int center, bool swing_left, double h0_cycle, double h1_cycle, double lift_peak_local) {
        const int i0 = std::max(0, center - swing_half_width);
        const int i1 = std::min(n_wp - 1, center + swing_half_width);
        for (int i = i0; i <= i1; ++i) {
            const double a = (i1 == i0) ? 1.0 : static_cast<double>(i - i0) / static_cast<double>(i1 - i0);
            const double z_ref = (1.0 - a) * h0_cycle + a * h1_cycle +
                                 SWING_NO_CONTACT_MARGIN + liftProfile(a, lift_peak_local);
            if (swing_left) {
                left_contact_wp[i] = 0;
                left_swing_z[i] = std::max(left_swing_z[i], z_ref);
            } else {
                right_contact_wp[i] = 0;
                right_swing_z[i] = std::max(right_swing_z[i], z_ref);
            }
        }
        const int td = std::min(n_wp - 1, i1 + 1);
        if (swing_left) touchdown_left.push_back(td);
        else touchdown_right.push_back(td);
    };
    bool next_swing_left = true;
    for (int c = 0; c < stair_cycles; ++c) {
        const int t = static_cast<int>(std::round(static_cast<double>(c + 1) * (n_wp - 1) / static_cast<double>(stair_cycles + 1)));
        const int c1 = std::max(0, std::min(n_wp - 1, t - swing_shift));
        const int c2 = std::max(0, std::min(n_wp - 1, t + swing_shift));
        const double h0_cycle = h_start + dz_base * static_cast<double>(c) / std::max(1, stair_cycles);
        const double h1_cycle = h_start + dz_base * static_cast<double>(c + 1) / std::max(1, stair_cycles);
        const double step_h = std::max(0.0, h1_cycle - h0_cycle);
        const double lift_peak_local = std::max(LIFT_PEAK, step_h + SWING_NO_CONTACT_MARGIN + 0.03);
        applySwingWindow(c1, next_swing_left, h0_cycle, h1_cycle, lift_peak_local);
        applySwingWindow(c2, !next_swing_left, h0_cycle, h1_cycle, lift_peak_local);
        next_swing_left = !next_swing_left;
    }

    int friction_fix_count = 0;
    for (int i = 0; i < n_wp; ++i) {
        if (left_contact_wp[i] > 0) {
            double x = xl_wp[i], y = yl_wp[i];
            if (projectToFrictionSafeXY(x, y, *terrain, FRICTION_MU)) {
                xl_wp[i] = x; yl_wp[i] = y; ++friction_fix_count;
            }
        }
        if (right_contact_wp[i] > 0) {
            double x = xr_wp[i], y = yr_wp[i];
            if (projectToFrictionSafeXY(x, y, *terrain, FRICTION_MU)) {
                xr_wp[i] = x; yr_wp[i] = y; ++friction_fix_count;
            }
        }
    }

    for (int i = 0; i < n_wp; ++i) {
        const double h_l = terrain->height(xl_wp[i], yl_wp[i]);
        const double h_r = terrain->height(xr_wp[i], yr_wp[i]);
        if (left_contact_wp[i] > 0) {
            zl_wp[i] = h_l + STANCE_CLEARANCE;
        } else {
            zl_wp[i] = std::max(h_l + SWING_NO_CONTACT_MARGIN, left_swing_z[i]);
        }
        if (right_contact_wp[i] > 0) {
            zr_wp[i] = h_r + STANCE_CLEARANCE;
        } else {
            zr_wp[i] = std::max(h_r + SWING_NO_CONTACT_MARGIN, right_swing_z[i]);
        }
    }
    for (int idx_td : touchdown_left) {
        idx_td = std::max(0, std::min(n_wp - 1, idx_td));
        left_contact_wp[idx_td] = 1;
        double x = xl_wp[idx_td], y = yl_wp[idx_td];
        (void)projectToFrictionSafeXY(x, y, *terrain, FRICTION_MU);
        xl_wp[idx_td] = x;
        yl_wp[idx_td] = y;
        zl_wp[idx_td] = terrain->height(x, y) + STANCE_CLEARANCE;
    }
    for (int idx_td : touchdown_right) {
        idx_td = std::max(0, std::min(n_wp - 1, idx_td));
        right_contact_wp[idx_td] = 1;
        double x = xr_wp[idx_td], y = yr_wp[idx_td];
        (void)projectToFrictionSafeXY(x, y, *terrain, FRICTION_MU);
        xr_wp[idx_td] = x;
        yr_wp[idx_td] = y;
        zr_wp[idx_td] = terrain->height(x, y) + STANCE_CLEARANCE;
    }

    std::vector<int> left_contact_sample(N_SAMPLE, 1), right_contact_sample(N_SAMPLE, 1);
    for (int k = 0; k < N_SAMPLE; ++k) {
        const int iw = static_cast<int>(std::round((n_wp - 1) * (static_cast<double>(k) / std::max(1, N_SAMPLE - 1))));
        left_contact_sample[k] = left_contact_wp[std::max(0, std::min(n_wp - 1, iw))];
        right_contact_sample[k] = right_contact_wp[std::max(0, std::min(n_wp - 1, iw))];
    }
    bool swing_left = true;
    int left_swing_cnt = 0, right_swing_cnt = 0;
    for (int i = 0; i < n_wp; ++i) {
        if (left_contact_wp[i] == 0) ++left_swing_cnt;
        if (right_contact_wp[i] == 0) ++right_swing_cnt;
    }
    if (left_swing_cnt == 0 && right_swing_cnt > 0) swing_left = false;
    if (right_swing_cnt == 0 && left_swing_cnt > 0) swing_left = true;

    const double psi0 = psi_wp.front();
    const double psi1 = psi_wp.back();
    std::vector<double> cbx(N_COEFF), cby(N_COEFF), cbz(N_COEFF), cpsi(N_COEFF);
    std::vector<double> clx(N_COEFF), cly(N_COEFF), clz(N_COEFF * N_SEG_Z);
    std::vector<double> crx(N_COEFF), cry(N_COEFF), crz(N_COEFF * N_SEG_Z);
    fitCoefficientsFromWaypoints(N_COEFF, N_SEG_Z, s_wp, xb_wp, yb_wp, zb_wp, psi_wp,
                                 xl_wp, yl_wp, zl_wp, xr_wp, yr_wp, zr_wp,
                                 basis, cbx, cby, cbz, cpsi, clx, cly, clz, crx, cry, crz);

    Ipopt::SmartPtr<SegTrajOptNLP> nlp = new SegTrajOptNLP(N_COEFF, N_SAMPLE, basis, N_SEG_Z);
    nlp->setTerrain(terrain);
    std::vector<Ipopt::Number> x0_vec(nlp->getNumVars(), 0.0);
    int off = 0;
    std::copy(cbx.begin(), cbx.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cby.begin(), cby.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cbz.begin(), cbz.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cpsi.begin(), cpsi.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(clx.begin(), clx.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cly.begin(), cly.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(clz.begin(), clz.end(), x0_vec.begin() + off); off += static_cast<int>(clz.size());
    std::copy(crx.begin(), crx.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cry.begin(), cry.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(crz.begin(), crz.end(), x0_vec.begin() + off);
    nlp->setInitialGuess(x0_vec.data());

    std::vector<double> Xb0, Yb0, Zb0, Psib0, Xl0, Yl0, Zl0, Xr0, Yr0, Zr0;
    nlp->extractTrajectory(x0_vec.data(), Xb0, Yb0, Zb0, Psib0, Xl0, Yl0, Zl0, Xr0, Yr0, Zr0);
    writeTrajectoryFile(out_dir + "initial_trajectory.txt", Xb0, Yb0, Zb0, Psib0, Xl0, Yl0, Zl0, Xr0, Yr0, Zr0);

    const double xmin = *std::min_element(centers_x.begin(), centers_x.end()) - TERRAIN_MARGIN;
    const double xmax = *std::max_element(centers_x.begin(), centers_x.end()) + TERRAIN_MARGIN;
    const double ymin = *std::min_element(centers_y.begin(), centers_y.end()) - TERRAIN_MARGIN;
    const double ymax = *std::max_element(centers_y.begin(), centers_y.end()) + TERRAIN_MARGIN;
    nlp->addConstraint(std::make_shared<BoundaryConstraint10Var>(
        N_COEFF, N_SEG_Z, x0, y0, psi0, x1, y1, psi1, CLEAR_Z, *terrain, basis, WHEELS_NOMINAL,
        BOUNDARY_BASE_XY_TOL, BOUNDARY_BASE_Z_TOL, BOUNDARY_BASE_YAW_TOL, BOUNDARY_WHEEL_XY_TOL, BOUNDARY_WHEEL_Z_TOL,
        BOUNDARY_BASE_Z_MIN, BOUNDARY_BASE_Z_MAX, BOUNDARY_WHEEL_Z_MIN, BOUNDARY_WHEEL_Z_MAX));
    nlp->addConstraint(std::make_shared<TrajectoryBoundsConstraint8Var>(N_COEFF, N_SAMPLE, xmin, xmax, ymin, ymax));
    nlp->addConstraint(std::make_shared<CoefficientNormConstraint>(nlp->getNumVars(), COEFF_NORM_MAX));
    nlp->addConstraint(std::make_shared<WheelBaseDistanceTerrainConstraint>(N_COEFF, N_SAMPLE, N_SEG_Z, WHEEL_BASE_MIN, WHEEL_BASE_MAX));
    nlp->addConstraint(std::make_shared<WheelsDistanceTerrainConstraint>(N_COEFF, N_SAMPLE, N_SEG_Z, WHEELS_MIN, WHEELS_MAX));
    nlp->addConstraint(std::make_shared<WheelTerrainPenetrationConstraint>(N_COEFF, N_SAMPLE, N_SEG_Z, terrain, MAX_PENETRATION));
    nlp->addConstraint(std::make_shared<PenetrationBasedForceConstraint>(
        N_COEFF, N_SAMPLE, N_SEG_Z, terrain, CONTACT_STIFFNESS, left_contact_sample, right_contact_sample,
        STANCE_FORCE_MIN, STANCE_FORCE_MAX, SWING_FORCE_MAX));
    if (N_SEG_Z > 1) {
        nlp->addConstraint(std::make_shared<PiecewiseZContinuityConstraint>(
            N_COEFF, N_SEG_Z, basis, PIECEWISE_Z_C0_TOL));
    }

    std::vector<std::pair<std::string, std::shared_ptr<ObjectiveBase>>> objective_terms;
    auto obj_x0 = std::make_shared<InitialGuessProximityObjective>(std::vector<double>(x0_vec.begin(), x0_vec.end()),
                                                                    INITIAL_GUESS_PROX_WEIGHT);
    auto obj_bh = std::make_shared<BaseHeightStabilityObjective>(N_COEFF, N_SAMPLE, terrain, CLEAR_Z,
                                                                 BASE_HEIGHT_STABILITY_WEIGHT);
    auto obj_yaw = std::make_shared<YawPathAlignmentObjective>(N_COEFF, N_SAMPLE, YAW_PATH_ALIGNMENT_WEIGHT);
    auto obj_lift = std::make_shared<SwingWheelLiftObjective>(N_COEFF, N_SAMPLE, N_SEG_Z, terrain, swing_left,
                                                              LIFT_PEAK, STANCE_CLEARANCE, SWING_LIFT_OBJECTIVE_WEIGHT);
    objective_terms.push_back({"InitialGuessProximity", obj_x0});
    objective_terms.push_back({"BaseHeightStability", obj_bh});
    objective_terms.push_back({"YawPathAlignment", obj_yaw});
    objective_terms.push_back({"SwingWheelLift", obj_lift});
    nlp->addObjective(obj_x0);
    nlp->addObjective(obj_bh);
    nlp->addObjective(obj_yaw);
    nlp->addObjective(obj_lift);

    std::cout << "--- Lift-Leg NLP ---\n";
    std::cout << "segment waypoints: " << n_wp
              << ", base dz: " << dz_base
              << ", stair cycles(20cm): " << stair_cycles << "\n";
    std::cout << "swing profile dominant wheel: " << (swing_left ? "left" : "right")
              << ", friction projection updates: " << friction_fix_count << "\n";
    std::cout << "friction-cone proxy: |grad(h)| <= mu, mu=" << FRICTION_MU << "\n";
    std::cout << "vars: " << (8 * N_COEFF + 2 * N_COEFF * N_SEG_Z)
              << " (z segments: " << N_SEG_Z << "), samples: " << N_SAMPLE << "\n";

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetStringValue("linear_solver", "mumps");
    app->Options()->SetIntegerValue("print_level", 3);
    app->Options()->SetIntegerValue("max_iter", 220);
    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetIntegerValue("limited_memory_max_history", 8);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    if (app->Initialize() != Ipopt::Solve_Succeeded) {
        std::cerr << "Ipopt initialization failed\n";
        return -1;
    }

    const auto t_solve_start = std::chrono::high_resolution_clock::now();
    const Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(nlp);
    const auto t_solve_done = std::chrono::high_resolution_clock::now();
    if (nlp->hasSolution()) {
        std::vector<double> Xb1, Yb1, Zb1, Psib1, Xl1, Yl1, Zl1, Xr1, Yr1, Zr1;
        nlp->extractTrajectory(nullptr, Xb1, Yb1, Zb1, Psib1, Xl1, Yl1, Zl1, Xr1, Yr1, Zr1);
        writeTrajectoryFile(out_dir + "final_trajectory.txt", Xb1, Yb1, Zb1, Psib1, Xl1, Yl1, Zl1, Xr1, Yr1, Zr1);
        const auto& sol = nlp->getSolution();
        double f_sum = 0.0;
        for (const auto& term : objective_terms) {
            const double val = term.second->eval_f(sol.data());
            f_sum += val;
            std::cout << term.first << ": " << val << "\n";
        }
        std::cout << "objective total: " << f_sum << "\n";
    }

    const auto t_total_done = std::chrono::high_resolution_clock::now();
    using Ms = std::chrono::duration<double, std::milli>;
    std::cout << std::fixed << std::setprecision(4)
              << "solve time: " << (std::chrono::duration_cast<Ms>(t_solve_done - t_solve_start).count() / 1000.0) << " s\n"
              << "total time: " << (std::chrono::duration_cast<Ms>(t_total_done - t_total_start).count() / 1000.0) << " s\n";
    std::cout << "Ipopt return: " << static_cast<int>(status) << "\n";
    bool ok = (status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level);
    if (!ok && nlp->hasSolution()) ok = true;
    return ok ? 0 : 1;
}
