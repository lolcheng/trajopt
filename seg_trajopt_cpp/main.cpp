/**
 * 分段轨迹优化：第一段 rollable 路径，8 组切比雪夫系数（基座 xyzψ + 左右轮 xy），
 * 轮心高度由 RBF 地形隐式给定。Ipopt 求解。
 */
#include <stddef.h>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "IpIpoptApplication.hpp"
#include "seg_trajopt_nlp.hpp"
#include "constraints/boundary_constraint_8var.hpp"
#include "constraints/wheel_base_distance_terrain_constraint.hpp"
#include "constraints/wheels_distance_terrain_constraint.hpp"
#include "constraints/trajectory_bounds_constraint_8var.hpp"
#include "constraints/rollable_region_constraint.hpp"
#include "constraints/coefficient_norm_constraint.hpp"
#include "objectives/initial_guess_proximity_objective.hpp"
#include "objectives/base_height_stability_objective.hpp"
#include "objectives/yaw_path_alignment_objective.hpp"
#include "rbf_terrain.hpp"
#include "json_reader.hpp"
#include "basis_function.hpp"
#include "chebyshev_basis.hpp"
#include "polynomial.hpp"
#include "waypoint_reader.hpp"
#include "segmentation_reader.hpp"

static void fitCoefficientsFromWaypoints(
    int n_coeff,
    const std::vector<double>& s,
    const std::vector<double>& xb, const std::vector<double>& yb,
    const std::vector<double>& zb, const std::vector<double>& psi,
    const std::vector<double>& xl, const std::vector<double>& yl,
    const std::vector<double>& xr, const std::vector<double>& yr,
    const std::shared_ptr<BasisFunction>& basis,
    std::vector<double>& cbx, std::vector<double>& cby, std::vector<double>& cbz, std::vector<double>& cpsi,
    std::vector<double>& clx, std::vector<double>& cly,
    std::vector<double>& crx, std::vector<double>& cry) {

    const int m = static_cast<int>(s.size());
    const int n = n_coeff;
    const double eps = 1e-12;

    std::vector<double> Phi(m * n);
    basis->computeBasisMatrix(s.data(), m, n, Phi.data());

    auto lsq = [&](const std::vector<double>& y, std::vector<double>& coeff) {
        std::vector<double> Q(m * n, 0.0), R(n * n, 0.0);
        for (int j = 0; j < n; ++j) {
            std::vector<double> v(m, 0.0);
            for (int i = 0; i < m; ++i) v[i] = Phi[i * n + j];
            for (int k = 0; k < j; ++k) {
                double rkj = 0.0;
                for (int i = 0; i < m; ++i) rkj += Q[i * n + k] * v[i];
                R[k * n + j] = rkj;
                for (int i = 0; i < m; ++i) v[i] -= rkj * Q[i * n + k];
            }
            double nv = 0.0;
            for (int i = 0; i < m; ++i) nv += v[i] * v[i];
            nv = std::sqrt(nv);
            R[j * n + j] = nv;
            if (nv > eps)
                for (int i = 0; i < m; ++i) Q[i * n + j] = v[i] / nv;
        }
        std::vector<double> Qt_y(n, 0.0);
        for (int j = 0; j < n; ++j) {
            double dot = 0.0;
            for (int i = 0; i < m; ++i) dot += Q[i * n + j] * y[i];
            Qt_y[j] = dot;
        }
        coeff.assign(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            coeff[i] = Qt_y[i];
            for (int j = i + 1; j < n; ++j) coeff[i] -= R[i * n + j] * coeff[j];
            if (std::abs(R[i * n + i]) > eps) coeff[i] /= R[i * n + i];
            else coeff[i] = 0.0;
        }
    };

    lsq(xb, cbx); lsq(yb, cby); lsq(zb, cbz); lsq(psi, cpsi);
    lsq(xl, clx); lsq(yl, cly);
    lsq(xr, crx); lsq(yr, cry);
}

static bool writeTrajectoryFile(const std::string& path,
                                const std::vector<double>& Xb, const std::vector<double>& Yb,
                                const std::vector<double>& Zb, const std::vector<double>& Psib,
                                const std::vector<double>& Xl, const std::vector<double>& Yl, const std::vector<double>& Zl,
                                const std::vector<double>& Xr, const std::vector<double>& Yr, const std::vector<double>& Zr) {
    const int n = static_cast<int>(Xb.size());
    if (n == 0 || Yb.size() != static_cast<size_t>(n)) return false;
    std::ofstream out(path);
    if (!out) return false;
    out << std::fixed << std::setprecision(10);
    out << "s xb yb zb psi xl yl zl xr yr zr\n";
    for (int i = 0; i < n; ++i) {
        double s = (n <= 1) ? 0.5 : static_cast<double>(i) / (n - 1);
        out << s << " " << Xb[i] << " " << Yb[i] << " " << Zb[i] << " " << Psib[i] << " "
            << Xl[i] << " " << Yl[i] << " " << Zl[i] << " "
            << Xr[i] << " " << Yr[i] << " " << Zr[i] << "\n";
    }
    return true;
}

static bool writeTerrainGrid(const std::string& path, const RBFTerrain& terrain,
                             double xmin, double xmax, double ymin, double ymax,
                             int nx, int ny) {
    if (nx < 2 || ny < 2) return false;
    std::ofstream out(path);
    if (!out) return false;
    out << std::fixed << std::setprecision(10);
    out << nx << " " << ny << " " << xmin << " " << xmax << " " << ymin << " " << ymax << "\n";
    for (int j = 0; j < ny; ++j) {
        double y = (ny == 1) ? ymin : ymin + (ymax - ymin) * static_cast<double>(j) / (ny - 1);
        for (int i = 0; i < nx; ++i) {
            double x = (nx == 1) ? xmin : xmin + (xmax - xmin) * static_cast<double>(i) / (nx - 1);
            out << terrain.height(x, y) << "\n";
        }
    }
    return true;
}

int main(int argc, char** argv) {
    const int N_COEFF = 5;
    const int N_SAMPLE = 15;
    const double CLEAR_Z = 0.8;
    const double WHEELS_NOMINAL = 0.5;
    const double WHEEL_BASE_MIN = 0.2, WHEEL_BASE_MAX = 1.2;
    const double WHEELS_MIN = 0.31622776601683794, WHEELS_MAX = 1.0;  // d_lr^2 in [0.1, 1]
    const double TERRAIN_MARGIN = 0.2;
    const double ROLLABLE_HARD_THRESHOLD = 0.1;
    const double ROLLABLE_HARD_BUFFER = 0.1;
    const bool ROLLABLE_INCLUDE_WHEELS = true;
    const double COEFF_NORM_MAX = 10000.0;
    const double BOUNDARY_BASE_XY_TOL = 1e-3;
    const double BOUNDARY_BASE_Z_TOL = 1e-3;
    const double BOUNDARY_BASE_YAW_TOL = 1e-3;
    const double BOUNDARY_WHEEL_XY_TOL = 5e-2;
    const double INITIAL_GUESS_PROX_WEIGHT = 1e-2;
    const double BASE_HEIGHT_STABILITY_WEIGHT = 2.0 / static_cast<double>(N_SAMPLE);
    const double YAW_PATH_ALIGNMENT_WEIGHT = 10.0 / static_cast<double>(N_SAMPLE);

    auto t_total_start = std::chrono::high_resolution_clock::now();

    std::string data_dir = (argc > 1) ? argv[1] : std::string("terrain_res");
    if (data_dir.back() != '/') data_dir += '/';
    std::string rbf_file = data_dir + std::string("rbf.json");
    std::string waypoints_file = data_dir + "waypoints_segmented.txt";
    std::string out_dir = (argc > 2) ? argv[2] : std::string("seg_out");
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
    if (!JSONReader::readRBFParams(rbf_file, centers_x, centers_y, weights, sigma)) {
        std::cerr << "Failed to read RBF from " << rbf_file << std::endl;
        return -1;
    }
    auto terrain = std::make_shared<RBFTerrain>(centers_x, centers_y, weights, sigma);

    std::vector<double> wp_x, wp_y;
    if (!WaypointReader::readFirstRollableSegment(waypoints_file, wp_x, wp_y)) {
        std::cerr << "Failed to read first rollable segment from " << waypoints_file << std::endl;
        return -1;
    }
    const int n_wp = static_cast<int>(wp_x.size());
    std::cout << "First rollable segment: " << n_wp << " waypoints" << std::endl;

    double x0 = wp_x.front(), y0 = wp_y.front();
    double x1 = wp_x.back(), y1 = wp_y.back();

    std::shared_ptr<BasisFunction> basis = std::make_shared<ChebyshevBasis>();
    Polynomial::set_basis(basis);

    std::vector<double> s_wp(n_wp), xb_wp(n_wp), yb_wp(n_wp), zb_wp(n_wp), psi_wp(n_wp);
    std::vector<double> xl_wp(n_wp), yl_wp(n_wp);
    std::vector<double> xr_wp(n_wp), yr_wp(n_wp);

    for (int i = 0; i < n_wp; ++i) {
        s_wp[i] = (n_wp <= 1) ? 0.5 : static_cast<double>(i) / (n_wp - 1);
        xb_wp[i] = wp_x[i];
        yb_wp[i] = wp_y[i];
        zb_wp[i] = terrain->height(wp_x[i], wp_y[i]) + CLEAR_Z;

        double dx = (i < n_wp - 1) ? (wp_x[i + 1] - wp_x[i]) : (wp_x[i] - wp_x[i - 1]);
        double dy = (i < n_wp - 1) ? (wp_y[i + 1] - wp_y[i]) : (wp_y[i] - wp_y[i - 1]);
        psi_wp[i] = std::atan2(dy, dx);

        double ux = -std::sin(psi_wp[i]);
        double uy = std::cos(psi_wp[i]);
        xl_wp[i] = xb_wp[i] + (WHEELS_NOMINAL / 2.0) * ux;
        yl_wp[i] = yb_wp[i] + (WHEELS_NOMINAL / 2.0) * uy;
        xr_wp[i] = xb_wp[i] - (WHEELS_NOMINAL / 2.0) * ux;
        yr_wp[i] = yb_wp[i] - (WHEELS_NOMINAL / 2.0) * uy;
    }

    double psi0 = psi_wp.front();
    double psi1 = psi_wp.back();

    std::vector<double> cbx(N_COEFF), cby(N_COEFF), cbz(N_COEFF), cpsi(N_COEFF);
    std::vector<double> clx(N_COEFF), cly(N_COEFF);
    std::vector<double> crx(N_COEFF), cry(N_COEFF);
    fitCoefficientsFromWaypoints(N_COEFF, s_wp, xb_wp, yb_wp, zb_wp, psi_wp,
                                 xl_wp, yl_wp, xr_wp, yr_wp,
                                 basis, cbx, cby, cbz, cpsi, clx, cly, crx, cry);

    Ipopt::SmartPtr<SegTrajOptNLP> nlp = new SegTrajOptNLP(N_COEFF, N_SAMPLE, basis);
    nlp->setTerrain(terrain);
    std::vector<std::pair<std::string, std::shared_ptr<ObjectiveBase>>> objective_terms;

    std::vector<Ipopt::Number> x0_vec(nlp->getNumVars(), 0.0);
    int off = 0;
    std::copy(cbx.begin(), cbx.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cby.begin(), cby.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cbz.begin(), cbz.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cpsi.begin(), cpsi.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(clx.begin(), clx.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cly.begin(), cly.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(crx.begin(), crx.end(), x0_vec.begin() + off); off += N_COEFF;
    std::copy(cry.begin(), cry.end(), x0_vec.begin() + off);
    nlp->setInitialGuess(x0_vec.data());

    auto t_init_done = std::chrono::high_resolution_clock::now();

    std::vector<double> Xb0, Yb0, Zb0, Psib0, Xl0, Yl0, Zl0, Xr0, Yr0, Zr0;
    nlp->extractTrajectory(x0_vec.data(), Xb0, Yb0, Zb0, Psib0, Xl0, Yl0, Zl0, Xr0, Yr0, Zr0);
    if (writeTrajectoryFile(out_dir + "initial_trajectory.txt", Xb0, Yb0, Zb0, Psib0, Xl0, Yl0, Zl0, Xr0, Yr0, Zr0))
        std::cout << "Initial trajectory written to " << out_dir << "initial_trajectory.txt" << std::endl;

    double xmin = *std::min_element(centers_x.begin(), centers_x.end()) - TERRAIN_MARGIN;
    double xmax = *std::max_element(centers_x.begin(), centers_x.end()) + TERRAIN_MARGIN;
    double ymin = *std::min_element(centers_y.begin(), centers_y.end()) - TERRAIN_MARGIN;
    double ymax = *std::max_element(centers_y.begin(), centers_y.end()) + TERRAIN_MARGIN;

    if (writeTerrainGrid(out_dir + "terrain_grid.txt", *terrain, xmin, xmax, ymin, ymax, 100, 100))
        std::cout << "Terrain grid written to " << out_dir << "terrain_grid.txt" << std::endl;

    auto boundary_constraint = std::make_shared<BoundaryConstraint8Var>(
        N_COEFF, x0, y0, psi0, x1, y1, psi1, CLEAR_Z, *terrain, basis, WHEELS_NOMINAL,
        BOUNDARY_BASE_XY_TOL, BOUNDARY_BASE_Z_TOL, BOUNDARY_BASE_YAW_TOL, BOUNDARY_WHEEL_XY_TOL);
    nlp->addConstraint(boundary_constraint);

    nlp->addConstraint(std::make_shared<TrajectoryBoundsConstraint8Var>(
        N_COEFF, N_SAMPLE, xmin, xmax, ymin, ymax));
    nlp->addConstraint(std::make_shared<CoefficientNormConstraint>(
        nlp->getNumVars(), COEFF_NORM_MAX));

    auto seg_reader = std::make_shared<SegmentationReader>();
    const std::string seg_file = data_dir + "terrain_segmentation_data.txt";
    const bool seg_loaded = seg_reader->load(seg_file);
    if (seg_loaded) {
        nlp->addConstraint(std::make_shared<RollableRegionConstraint>(
            N_COEFF, N_SAMPLE, seg_reader, ROLLABLE_HARD_THRESHOLD, ROLLABLE_INCLUDE_WHEELS, ROLLABLE_HARD_BUFFER));
        std::cout << "Rollable constraint: " << seg_file << std::endl;
    } else {
        std::cerr << "Warning: no segmentation file (" << seg_file << "), skipping rollable hard constraint.\n";
    }

    nlp->addConstraint(std::make_shared<WheelBaseDistanceTerrainConstraint>(
        N_COEFF, N_SAMPLE, terrain, WHEEL_BASE_MIN, WHEEL_BASE_MAX));
    nlp->addConstraint(std::make_shared<WheelsDistanceTerrainConstraint>(
        N_COEFF, N_SAMPLE, terrain, WHEELS_MIN, WHEELS_MAX));

    {
        auto obj_x0 = std::make_shared<InitialGuessProximityObjective>(
            std::vector<double>(x0_vec.begin(), x0_vec.end()), INITIAL_GUESS_PROX_WEIGHT);
        objective_terms.push_back({"InitialGuessProximity", obj_x0});
        nlp->addObjective(obj_x0);
    }
    {
        auto obj_bh = std::make_shared<BaseHeightStabilityObjective>(
            N_COEFF, N_SAMPLE, terrain, CLEAR_Z, BASE_HEIGHT_STABILITY_WEIGHT);
        objective_terms.push_back({"BaseHeightStability", obj_bh});
        nlp->addObjective(obj_bh);
    }
    {
        auto obj_yaw = std::make_shared<YawPathAlignmentObjective>(
            N_COEFF, N_SAMPLE, YAW_PATH_ALIGNMENT_WEIGHT);
        objective_terms.push_back({"YawPathAlignment", obj_yaw});
        nlp->addObjective(obj_yaw);
    }

    std::cout << "\n--- NLP summary ---\n"
              << "  variables: 8 x " << N_COEFF << " Chebyshev coeffs; samples: " << N_SAMPLE << "\n"
              << "  hard: boundary (box), traj bounds, coeff mean-square <= " << COEFF_NORM_MAX << "^2"
              << (seg_loaded ? ", rollable" : "") << ", wheel-base, wheel-wheel\n"
              << "  objectives: initial-guess proximity, base height vs terrain, yaw-path alignment\n\n";

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetStringValue("linear_solver", "mumps");
    app->Options()->SetIntegerValue("print_level", 3);
    app->Options()->SetIntegerValue("max_iter", 200);
    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("acceptable_tol", 0.2);
    app->Options()->SetIntegerValue("acceptable_iter", 5);
    app->Options()->SetNumericValue("acceptable_dual_inf_tol", 1e9);
    app->Options()->SetNumericValue("acceptable_constr_viol_tol", 1e-4);
    app->Options()->SetNumericValue("acceptable_compl_inf_tol", 1e-2);
    app->Options()->SetNumericValue("acceptable_obj_change_tol", 1e-4);
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetIntegerValue("limited_memory_max_history", 6);
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    app->Options()->SetNumericValue("nlp_scaling_max_gradient", 50.0);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetIntegerValue("watchdog_shortened_iter_trigger", 4);

    if (app->Initialize() != Ipopt::Solve_Succeeded) {
        std::cerr << "Ipopt initialization failed" << std::endl;
        return -1;
    }

    std::cout << "=== Solving ===" << std::endl;
    auto t_solve_start = std::chrono::high_resolution_clock::now();
    Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(nlp);
    auto t_solve_done = std::chrono::high_resolution_clock::now();

    auto t_write_start = std::chrono::high_resolution_clock::now();
    if (nlp->hasSolution()) {
        std::vector<double> Xb1, Yb1, Zb1, Psib1, Xl1, Yl1, Zl1, Xr1, Yr1, Zr1;
        nlp->extractTrajectory(nullptr, Xb1, Yb1, Zb1, Psib1, Xl1, Yl1, Zl1, Xr1, Yr1, Zr1);
        if (writeTrajectoryFile(out_dir + "final_trajectory.txt", Xb1, Yb1, Zb1, Psib1, Xl1, Yl1, Zl1, Xr1, Yr1, Zr1))
            std::cout << "Final trajectory written to " << out_dir << "final_trajectory.txt" << std::endl;

        const auto& sol = nlp->getSolution();
        if (!sol.empty()) {
            double f_total = 0.0;
            std::cout << "--- Objectives (at solution) ---\n";
            for (const auto& term : objective_terms) {
                double f = term.second->eval_f(sol.data());
                f_total += f;
                std::cout << "  " << term.first << ": " << f << "\n";
            }
            std::cout << "  total: " << f_total << "\n";

            std::vector<double> Xb_d, Yb_d, Zb_d, Psib_d, Xl_d, Yl_d, Zl_d, Xr_d, Yr_d, Zr_d;
            nlp->extractTrajectory(sol.data(), Xb_d, Yb_d, Zb_d, Psib_d, Xl_d, Yl_d, Zl_d, Xr_d, Yr_d, Zr_d);
            double max_v_boundary = 0.0;
            double max_v_bounds = 0.0;
            double max_v_roll = 0.0;
            double max_v_wb = 0.0;
            double max_v_ww = 0.0;

            auto upd_eq = [](double v, double& m) { m = std::max(m, std::abs(v)); };
            auto upd_lbub = [](double v, double lb, double ub, double& m) {
                if (v < lb) m = std::max(m, lb - v);
                if (v > ub) m = std::max(m, v - ub);
            };

            upd_eq(Xb_d.front() - x0, max_v_boundary);
            upd_eq(Yb_d.front() - y0, max_v_boundary);
            upd_eq(Zb_d.front() - (terrain->height(x0, y0) + CLEAR_Z), max_v_boundary);
            upd_eq(Psib_d.front() - psi0, max_v_boundary);
            upd_eq(Xb_d.back() - x1, max_v_boundary);
            upd_eq(Yb_d.back() - y1, max_v_boundary);
            upd_eq(Zb_d.back() - (terrain->height(x1, y1) + CLEAR_Z), max_v_boundary);
            upd_eq(Psib_d.back() - psi1, max_v_boundary);
            auto wheel_xy_res = [&](int i, bool left) {
                double ux = -std::sin(Psib_d[i]);
                double uy = std::cos(Psib_d[i]);
                double tx = Xb_d[i] + (left ? 1.0 : -1.0) * (WHEELS_NOMINAL / 2.0) * ux;
                double ty = Yb_d[i] + (left ? 1.0 : -1.0) * (WHEELS_NOMINAL / 2.0) * uy;
                double xw = left ? Xl_d[i] : Xr_d[i];
                double yw = left ? Yl_d[i] : Yr_d[i];
                upd_eq(xw - tx, max_v_boundary);
                upd_eq(yw - ty, max_v_boundary);
            };
            wheel_xy_res(0, true);
            wheel_xy_res(0, false);
            wheel_xy_res(nlp->get_n_sample() - 1, true);
            wheel_xy_res(nlp->get_n_sample() - 1, false);

            for (int i = 0; i < nlp->get_n_sample(); ++i) {
                upd_lbub(Xb_d[i], xmin, xmax, max_v_bounds);
                upd_lbub(Yb_d[i], ymin, ymax, max_v_bounds);
                upd_lbub(Xl_d[i], xmin, xmax, max_v_bounds);
                upd_lbub(Yl_d[i], ymin, ymax, max_v_bounds);
                upd_lbub(Xr_d[i], xmin, xmax, max_v_bounds);
                upd_lbub(Yr_d[i], ymin, ymax, max_v_bounds);
            }

            if (seg_loaded && seg_reader->getNx() > 1 && seg_reader->getNy() > 1) {
                for (int i = 0; i < nlp->get_n_sample(); ++i) {
                    double ell = seg_reader->interp(Xb_d[i], Yb_d[i]);
                    double g_roll = ROLLABLE_HARD_THRESHOLD - ell;
                    if (g_roll > ROLLABLE_HARD_BUFFER)
                        max_v_roll = std::max(max_v_roll, g_roll - ROLLABLE_HARD_BUFFER);
                    if (ROLLABLE_INCLUDE_WHEELS) {
                        double g_l = ROLLABLE_HARD_THRESHOLD - seg_reader->interp(Xl_d[i], Yl_d[i]);
                        double g_r = ROLLABLE_HARD_THRESHOLD - seg_reader->interp(Xr_d[i], Yr_d[i]);
                        if (g_l > ROLLABLE_HARD_BUFFER)
                            max_v_roll = std::max(max_v_roll, g_l - ROLLABLE_HARD_BUFFER);
                        if (g_r > ROLLABLE_HARD_BUFFER)
                            max_v_roll = std::max(max_v_roll, g_r - ROLLABLE_HARD_BUFFER);
                    }
                }
            }

            for (int i = 0; i < nlp->get_n_sample(); ++i) {
                double dxl = Xb_d[i] - Xl_d[i], dyl = Yb_d[i] - Yl_d[i], dzl = Zb_d[i] - Zl_d[i];
                double dxr = Xb_d[i] - Xr_d[i], dyr = Yb_d[i] - Yr_d[i], dzr = Zb_d[i] - Zr_d[i];
                double d2l = dxl * dxl + dyl * dyl + dzl * dzl;
                double d2r = dxr * dxr + dyr * dyr + dzr * dzr;
                upd_lbub(d2l, WHEEL_BASE_MIN * WHEEL_BASE_MIN, WHEEL_BASE_MAX * WHEEL_BASE_MAX, max_v_wb);
                upd_lbub(d2r, WHEEL_BASE_MIN * WHEEL_BASE_MIN, WHEEL_BASE_MAX * WHEEL_BASE_MAX, max_v_wb);
                double dx = Xr_d[i] - Xl_d[i], dy = Yr_d[i] - Yl_d[i], dz = Zr_d[i] - Zl_d[i];
                double d2w = dx * dx + dy * dy + dz * dz;
                upd_lbub(d2w, WHEELS_MIN * WHEELS_MIN, WHEELS_MAX * WHEELS_MAX, max_v_ww);
            }

            std::cout << "--- Constraint violation (diagnostic, endpoint-style boundary) ---\n"
                      << "  boundary (endpoint eq): " << max_v_boundary << "\n"
                      << "  trajectory bounds:     " << max_v_bounds << "\n"
                      << "  rollable hard:         ";
            if (seg_loaded)
                std::cout << max_v_roll;
            else
                std::cout << "n/a";
            std::cout << "\n  wheel-base interval:   " << max_v_wb << "\n"
                      << "  wheel-wheel interval:  " << max_v_ww << "\n";
        }
    }

    auto t_total_done = std::chrono::high_resolution_clock::now();
    using Ms = std::chrono::duration<double, std::milli>;
    std::cout << "--- Time ---\n"
              << std::fixed << std::setprecision(4)
              << "  init:  " << (std::chrono::duration_cast<Ms>(t_init_done - t_total_start).count() / 1000.0) << " s\n"
              << "  solve: " << (std::chrono::duration_cast<Ms>(t_solve_done - t_solve_start).count() / 1000.0) << " s\n"
              << "  post:  " << (std::chrono::duration_cast<Ms>(t_total_done - t_write_start).count() / 1000.0) << " s\n"
              << "  total: " << (std::chrono::duration_cast<Ms>(t_total_done - t_total_start).count() / 1000.0) << " s\n";
    std::cout << "Ipopt return: " << static_cast<int>(status) << "\n";
    std::cout << "Visualize: python3 scripts/plot_seg_trajectory.py " << out_dir << "\n";

    bool ok = (status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level);
    if (!ok && nlp->hasSolution()) {
        ok = true;
        std::cout << "(Using feasible iterate though Ipopt status is not success.)\n";
    }
    return ok ? 0 : 1;
}
