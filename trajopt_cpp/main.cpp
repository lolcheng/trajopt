#include <stddef.h>  // Must be before Ipopt headers (C-style header)
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "trajopt_nlp.hpp"
#include "constraints/boundary_constraint.hpp"
#include "constraints/wheel_base_distance_constraint.hpp"
#include "constraints/wheels_distance_constraint.hpp"
#include "constraints/anti_crossing_constraint.hpp"
#include "constraints/trajectory_bounds_constraint.hpp"
#include "constraints/normal_force_upper_bound_constraint.hpp"
#include "constraints/wheel_airborne_constraint.hpp"
#include "objectives/wheel_terrain_distance_objective.hpp"
#include "objectives/wheel_airborne_penalty_objective.hpp"
#include "objectives/path_length_objective.hpp"
#include "objectives/forward_progress_objective.hpp"
#include "objectives/regularization_objective.hpp"
#include "objectives/smoothing_objective.hpp"
#include "objectives/wheel_base_distance_nominal_objective.hpp"
#include "objectives/wheels_distance_nominal_objective.hpp"
#include "objectives/yaw_path_alignment_objective.hpp"
#include "objectives/contact_stiffness_objective.hpp"
#include "objectives/force_balance_objective.hpp"
#include "objectives/torque_balance_objective.hpp"
#include "objectives/friction_cone_objective.hpp"
#include "rbf_terrain.hpp"
#include <memory>
#include "json_reader.hpp"
#include "terrain_visualizer.hpp"
#include "trajectory_visualizer.hpp"
#include "force_data_saver.hpp"
#include "basis_function.hpp"
#include "polynomial_basis.hpp"
#include "chebyshev_basis.hpp"

int main(int argc, char** argv) {
    // 参数设置
    const int N_COEFF = 5;
    const int N_SAMPLE = 25;  // 采样点数量：25个点，间距约20cm（5m/25=0.2m）
    
    // 调试设置
    const bool DEBUG_MODE = false;  // 是否输出详细调试信息（设为true会输出大量信息）
    const bool VERBOSE_IPOPT = false;  // 是否输出Ipopt详细日志（设为true会输出大量迭代信息）
    
    // 可视化设置
    const bool VISUALIZATION_AUTO_CLOSE = false;  // 是否自动关闭可视化窗口
    const double VISUALIZATION_CLOSE_TIME = 5.0;  // 自动关闭时间（秒）
    
    // 起点和终点
    const double XB0 = 2.5, YB0 = 1.0, PSI0 = -M_PI / 2.0;
    const double XB1 = 2.5, YB1 = -3.5, PSI1 = -M_PI / 2.0;
    const double CLEAR_Z = 0.8;
    
    // 从JSON文件读取RBF参数
    std::string json_file = "/home/yizhe/trajopt/result/RBF/1753771011_476223707.json";
    if (argc > 1) {
        json_file = argv[1];  // 允许通过命令行参数指定JSON文件
    }
    
    std::vector<double> centers_x, centers_y, weights;
    double RBF_SIGMA;
    
    std::cout << "Reading RBF parameters from: " << json_file << std::endl;
    if (!JSONReader::readRBFParams(json_file, centers_x, centers_y, weights, RBF_SIGMA)) {
        std::cerr << "Failed to read RBF parameters from JSON file!" << std::endl;
        return -1;
    }
    
    // 创建RBF地形
    auto terrain = std::make_shared<RBFTerrain>(centers_x, centers_y, weights, RBF_SIGMA);
    
    // 可视化RBF地形
    if (DEBUG_MODE) {
        std::cout << "\n=== Visualizing RBF Terrain ===" << std::endl;
    }
    const double TERRAIN_MARGIN = 0.2;
    TerrainVisualizer::visualize(*terrain, 161, 161, TERRAIN_MARGIN,
                                 VISUALIZATION_AUTO_CLOSE, VISUALIZATION_CLOSE_TIME);
    if (DEBUG_MODE) {
        std::cout << "=== Visualization Complete ===\n" << std::endl;
    }
    
    // 参数定义（需要在计算初始猜测之前定义）
    const double WHEELS_NOMINAL = 0.5;
    const double MASS = 1.0;
    const double GRAVITY = 9.81;
    
    // 统一使用切比雪夫基
    std::shared_ptr<BasisFunction> basis = std::make_shared<ChebyshevBasis>();
    std::cout << "Using Chebyshev basis functions" << std::endl;
    // 关键：同步设置 Polynomial 包装器的全局基函数，
    // 这样所有仍通过 Polynomial:: 接口实现的约束/目标都会与当前基函数一致。
    Polynomial::set_basis(basis);
    
    // 创建NLP问题（使用Ipopt的SmartPtr）
    Ipopt::SmartPtr<TrajOptNLP> nlp = new TrajOptNLP(N_COEFF, N_SAMPLE, basis);
    nlp->setTerrain(terrain);
    
    // 计算初始猜测（基于边界条件和线性插值）
    if (DEBUG_MODE) {
        std::cout << "\n=== Computing Initial Guess ===" << std::endl;
    }
    const int n_init_points = 7;
    std::vector<double> s_init(n_init_points);
    // 使用切比雪夫-Lobatto节点（映射到[0,1]）初始化，
    // 相比均匀采样，端点附近分辨率更高，更匹配切比雪夫基。
    for (int i = 0; i < n_init_points; ++i) {
        double theta = M_PI * static_cast<double>(i) / static_cast<double>(n_init_points - 1);
        double x_cheb = std::cos(theta);         // [-1, 1]
        s_init[i] = 0.5 * (1.0 - x_cheb);        // [0, 1], 单调递增且包含端点
    }
    
    // Base轨迹初始猜测
    std::vector<double> xb_init(n_init_points), yb_init(n_init_points), zb_init(n_init_points), psi_init(n_init_points);
    for (int i = 0; i < n_init_points; ++i) {
        xb_init[i] = XB0 + (XB1 - XB0) * s_init[i];
        yb_init[i] = YB0 + (YB1 - YB0) * s_init[i];
        psi_init[i] = PSI0 + (PSI1 - PSI0) * s_init[i];
        // Z坐标：使用地形高度 + CLEAR_Z
        double xb = xb_init[i];
        double yb = yb_init[i];
        double h = terrain->height(xb, yb);
        zb_init[i] = h + CLEAR_Z;
    }
    
    // 轮子轨迹初始猜测
    std::vector<double> xl_init(n_init_points), yl_init(n_init_points), zl_init(n_init_points);
    std::vector<double> xr_init(n_init_points), yr_init(n_init_points), zr_init(n_init_points);
    for (int i = 0; i < n_init_points; ++i) {
        double ux = -std::sin(psi_init[i]);
        double uy = std::cos(psi_init[i]);
        xl_init[i] = xb_init[i] + (WHEELS_NOMINAL / 2.0) * ux;
        yl_init[i] = yb_init[i] + (WHEELS_NOMINAL / 2.0) * uy;
        zl_init[i] = terrain->height(xl_init[i], yl_init[i]);
        xr_init[i] = xb_init[i] - (WHEELS_NOMINAL / 2.0) * ux;
        yr_init[i] = yb_init[i] - (WHEELS_NOMINAL / 2.0) * uy;
        zr_init[i] = terrain->height(xr_init[i], yr_init[i]);
    }
    
    // 使用最小二乘拟合计算系数
    // 创建基矩阵（使用当前选定的基函数）
    std::vector<double> Phi_init(n_init_points * N_COEFF);
    basis->computeBasisMatrix(s_init.data(), n_init_points, N_COEFF, Phi_init.data());
    
    // 求解最小二乘问题：Phi * coeff ≈ y
    // 使用改进的Gram-Schmidt QR分解（比法方程更稳定）
    std::vector<double> cbx_init(N_COEFF), cby_init(N_COEFF), cbz_init(N_COEFF), cpsi_init(N_COEFF);
    std::vector<double> clx_init(N_COEFF), cly_init(N_COEFF), clz_init(N_COEFF);
    std::vector<double> crx_init(N_COEFF), cry_init(N_COEFF), crz_init(N_COEFF);
    
    auto lsq_fit = [&](const std::vector<double>& y, std::vector<double>& coeff) {
        const int m = n_init_points;
        const int n = N_COEFF;
        const double eps_rank = 1e-12;

        // Q: m x n（列正交），R: n x n（上三角）
        std::vector<double> Q(m * n, 0.0);
        std::vector<double> R(n * n, 0.0);

        // 改进Gram-Schmidt：对 Phi_init 的列做正交化
        for (int j = 0; j < n; ++j) {
            std::vector<double> v(m, 0.0);
            for (int i = 0; i < m; ++i) {
                v[i] = Phi_init[i * n + j];
            }

            for (int k = 0; k < j; ++k) {
                double r_kj = 0.0;
                for (int i = 0; i < m; ++i) {
                    r_kj += Q[i * n + k] * v[i];
                }
                R[k * n + j] = r_kj;
                for (int i = 0; i < m; ++i) {
                    v[i] -= r_kj * Q[i * n + k];
                }
            }

            double norm_v = 0.0;
            for (int i = 0; i < m; ++i) {
                norm_v += v[i] * v[i];
            }
            norm_v = std::sqrt(norm_v);
            R[j * n + j] = norm_v;

            if (norm_v > eps_rank) {
                for (int i = 0; i < m; ++i) {
                    Q[i * n + j] = v[i] / norm_v;
                }
            } else {
                // 退化列（近线性相关）时置零，回代阶段会稳定处理
                for (int i = 0; i < m; ++i) {
                    Q[i * n + j] = 0.0;
                }
            }
        }

        // 计算 Q^T y
        std::vector<double> Qt_y(n, 0.0);
        for (int j = 0; j < n; ++j) {
            double dot = 0.0;
            for (int i = 0; i < m; ++i) {
                dot += Q[i * n + j] * y[i];
            }
            Qt_y[j] = dot;
        }

        // 回代解 R c = Q^T y
        coeff.assign(N_COEFF, 0.0);
        for (int i = N_COEFF - 1; i >= 0; --i) {
            coeff[i] = Qt_y[i];
            for (int j = i + 1; j < N_COEFF; ++j) {
                coeff[i] -= R[i * N_COEFF + j] * coeff[j];
            }
            double pivot = R[i * N_COEFF + i];
            if (std::abs(pivot) > eps_rank) {
                coeff[i] /= pivot;
            } else {
                // 近奇异时稳定回退
                coeff[i] = 0.0;
            }

            if (!std::isfinite(coeff[i])) {
                std::cerr << "WARNING: NaN/Inf in coefficient fitting, setting to 0" << std::endl;
                coeff[i] = 0.0;
            }
        }
    };
    
    lsq_fit(xb_init, cbx_init);
    lsq_fit(yb_init, cby_init);
    lsq_fit(zb_init, cbz_init);
    lsq_fit(psi_init, cpsi_init);
    lsq_fit(xl_init, clx_init);
    lsq_fit(yl_init, cly_init);
    lsq_fit(zl_init, clz_init);
    lsq_fit(xr_init, crx_init);
    lsq_fit(yr_init, cry_init);
    lsq_fit(zr_init, crz_init);
    
    // 构建完整的初始猜测向量
    int n_vars = nlp->getNumVars();
    std::vector<Ipopt::Number> x0(n_vars, 0.0);
    
    // 设置多项式系数
    int offset = 0;
    std::copy(cbx_init.begin(), cbx_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(cby_init.begin(), cby_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(cbz_init.begin(), cbz_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(cpsi_init.begin(), cpsi_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(clx_init.begin(), clx_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(cly_init.begin(), cly_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(clz_init.begin(), clz_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(crx_init.begin(), crx_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(cry_init.begin(), cry_init.end(), x0.begin() + offset); offset += N_COEFF;
    std::copy(crz_init.begin(), crz_init.end(), x0.begin() + offset); offset += N_COEFF;
    
    // 设置接触力初始值（每个轮子承担一半重力，法向向上）
    double f_nominal = MASS * GRAVITY / 2.0;  // 每个轮子的法向力
    
    // 接触力初始值：fLz和fRz设为f_nominal，其他为0
    for (int i = 0; i < N_SAMPLE; ++i) {
        x0[offset + i] = 0.0;  // fLx
    }
    offset += N_SAMPLE;
    for (int i = 0; i < N_SAMPLE; ++i) {
        x0[offset + i] = 0.0;  // fLy
    }
    offset += N_SAMPLE;
    for (int i = 0; i < N_SAMPLE; ++i) {
        x0[offset + i] = f_nominal;  // fLz
    }
    offset += N_SAMPLE;
    for (int i = 0; i < N_SAMPLE; ++i) {
        x0[offset + i] = 0.0;  // fRx
    }
    offset += N_SAMPLE;
    for (int i = 0; i < N_SAMPLE; ++i) {
        x0[offset + i] = 0.0;  // fRy
    }
    offset += N_SAMPLE;
    for (int i = 0; i < N_SAMPLE; ++i) {
        x0[offset + i] = f_nominal;  // fRz
    }
    
    // 检查初始猜测是否有NaN/Inf
    bool has_invalid = false;
    for (size_t i = 0; i < x0.size(); ++i) {
        if (!std::isfinite(x0[i])) {
            std::cerr << "ERROR: Invalid number in initial guess x0[" << i << "] = " << x0[i] << std::endl;
            has_invalid = true;
        }
    }
    if (has_invalid) {
        std::cerr << "ERROR: Initial guess contains NaN or Inf values!" << std::endl;
        return -1;
    }
    
    // 设置初始猜测
    nlp->setInitialGuess(x0.data());
    if (DEBUG_MODE) {
        std::cout << "Initial guess computed and set." << std::endl;
        // 输出初始猜测的一些统计信息
        double max_coeff = 0.0;
        for (size_t i = 0; i < 10 * N_COEFF; ++i) {
            max_coeff = std::max(max_coeff, std::abs(x0[i]));
        }
        std::cout << "  Max absolute coefficient value: " << max_coeff << std::endl;
    }
    
    // 参数定义
    const double WHEEL_BASE_MIN = 0.4;
    const double WHEEL_BASE_MAX = 1.2;
    const double WHEEL_BASE_NOMINAL = 0.8;
    const double WHEELS_MIN = 0.1;
    const double WHEELS_MAX = 0.8;

    // 地形图显示范围（与TerrainVisualizer一致），用于硬边界约束
    double xmin = *std::min_element(centers_x.begin(), centers_x.end()) - TERRAIN_MARGIN;
    double xmax = *std::max_element(centers_x.begin(), centers_x.end()) + TERRAIN_MARGIN;
    double ymin = *std::min_element(centers_y.begin(), centers_y.end()) - TERRAIN_MARGIN;
    double ymax = *std::max_element(centers_y.begin(), centers_y.end()) + TERRAIN_MARGIN;
    
    // 添加边界约束（包括轮子初末点约束）
    auto boundary_constraint = std::make_shared<BoundaryConstraint>(
        N_COEFF, XB0, YB0, PSI0, XB1, YB1, PSI1, CLEAR_Z, *terrain, basis, WHEELS_NOMINAL);
    nlp->addConstraint(boundary_constraint);

    // 添加轨迹平面边界约束：严格不超过地形图范围
    auto traj_bounds_constraint = std::make_shared<TrajectoryBoundsConstraint>(
        N_COEFF, N_SAMPLE, xmin, xmax, ymin, ymax);
    nlp->addConstraint(traj_bounds_constraint);
    
    // 添加轮子-中心距离约束
    auto wheel_base_constraint = std::make_shared<WheelBaseDistanceConstraint>(
        N_COEFF, N_SAMPLE, WHEEL_BASE_MIN, WHEEL_BASE_MAX);
    nlp->addConstraint(wheel_base_constraint);
    
    // 添加左右轮距离约束
    auto wheels_distance_constraint = std::make_shared<WheelsDistanceConstraint>(
        N_COEFF, N_SAMPLE, WHEELS_MIN, WHEELS_MAX);
    nlp->addConstraint(wheels_distance_constraint);
    
    // 添加防交叉约束
    auto anti_crossing_constraint = std::make_shared<AntiCrossingConstraint>(
        N_COEFF, N_SAMPLE);
    nlp->addConstraint(anti_crossing_constraint);
    
    // 添加目标函数项
    // 平滑项
    auto smoothing_obj = std::make_shared<SmoothingObjective>(N_COEFF, N_SAMPLE);
    nlp->addObjective(smoothing_obj);
    
    // 正则化项
    auto regularization_obj = std::make_shared<RegularizationObjective>(N_COEFF);
    nlp->addObjective(regularization_obj);
    
    // 轮子-中心距离标称项
    auto wheel_base_nominal_obj = std::make_shared<WheelBaseDistanceNominalObjective>(
        N_COEFF, N_SAMPLE, WHEEL_BASE_NOMINAL);
    nlp->addObjective(wheel_base_nominal_obj);
    
    // 左右轮距离标称项
    auto wheels_nominal_obj = std::make_shared<WheelsDistanceNominalObjective>(
        N_COEFF, N_SAMPLE, WHEELS_NOMINAL);
    nlp->addObjective(wheels_nominal_obj);
    
    // 轮子-地形距离项（降低权重，避免与离地惩罚冲突）
    // 注意：这个目标函数会惩罚所有距离（包括下陷），可能导致优化器为避免下陷而让轮子离地
    // 因此降低权重，主要依靠离地惩罚项来控制离地
    const double W_WHEEL_TERRAIN = 20.0;  // 降低权重，避免与离地惩罚冲突
    auto wheel_terrain_obj = std::make_shared<WheelTerrainDistanceObjective>(
        N_COEFF, N_SAMPLE, *terrain, W_WHEEL_TERRAIN / N_SAMPLE);
    nlp->addObjective(wheel_terrain_obj);
    
    // 轮子离地惩罚项（惩罚轮子离地，让轮子离地时间尽可能短）
    // 当轮子离地超过5cm时（Zl - Hl > 0.05），惩罚离地距离的平方
    // 这是主要控制离地的目标函数，使用较大权重
    const double W_AIRBORNE_PENALTY = 1000.0;  // 大幅增大权重，强烈减少轮子离地时间
    const double AIRBORNE_THRESHOLD = 0.05;  // 离地判据阈值：5cm
    auto wheel_airborne_penalty_obj = std::make_shared<WheelAirbornePenaltyObjective>(
        N_COEFF, N_SAMPLE, terrain, W_AIRBORNE_PENALTY / N_SAMPLE, AIRBORNE_THRESHOLD);
    nlp->addObjective(wheel_airborne_penalty_obj);
    
    // 轨迹长度惩罚项（惩罚轨迹长度，让轨迹尽可能短）
    const double W_PATH_LENGTH = 10.0;  // 权重，让轨迹尽可能短
    auto path_length_obj = std::make_shared<PathLengthObjective>(
        N_COEFF, N_SAMPLE, W_PATH_LENGTH);
    nlp->addObjective(path_length_obj);

    // 前向进度惩罚项（减少回头和蛇形转弯）
    // 仅惩罚沿起点->终点方向的负进度
    const double W_FORWARD_PROGRESS = 300.0;
    auto forward_progress_obj = std::make_shared<ForwardProgressObjective>(
        N_COEFF, N_SAMPLE, XB0, YB0, XB1, YB1, W_FORWARD_PROGRESS / N_SAMPLE);
    nlp->addObjective(forward_progress_obj);
    
    // 轮子离地距离约束（允许下陷和离地，用于上楼梯等场景）
    // 下界：允许下陷（min_ground_distance < 0，用于弹簧模型）
    // 上界：允许离地（max_airborne_distance 较大，用于上楼梯）
    // 注意：这个约束主要用于防止极端情况，实际通过目标函数让轮子尽量贴地
    const double MIN_GROUND_DISTANCE = -0.05;   // -5cm（允许轻微下陷，用于弹簧模型）
    const double MAX_AIRBORNE_DISTANCE = 0.15;   // 15cm（进一步收紧上界，减少允许的离地距离）
    auto wheel_airborne_constraint = std::make_shared<WheelAirborneConstraint>(
        N_COEFF, N_SAMPLE, terrain, MIN_GROUND_DISTANCE, MAX_AIRBORNE_DISTANCE);
    nlp->addConstraint(wheel_airborne_constraint);
    
    // 偏航-路径对齐项
    auto yaw_alignment_obj = std::make_shared<YawPathAlignmentObjective>(N_COEFF, N_SAMPLE);
    nlp->addObjective(yaw_alignment_obj);
    
    // 法向力上限约束
    const double MAX_ACCEL = 1.0;
    auto normal_force_constraint = std::make_shared<NormalForceUpperBoundConstraint>(
        N_COEFF, N_SAMPLE, terrain, MASS, GRAVITY, MAX_ACCEL);
    nlp->addConstraint(normal_force_constraint);
    
    // 接触刚度项（确保接触力与穿透量符合弹簧模型）
    const double K_CONTACT = 500.0;
    const double W_CONTACT_STIFF = 100.0;  // 增大权重，确保接触力合理
    auto contact_stiffness_obj = std::make_shared<ContactStiffnessObjective>(
        N_COEFF, N_SAMPLE, terrain, K_CONTACT, W_CONTACT_STIFF / N_SAMPLE);
    nlp->addObjective(contact_stiffness_obj);
    
    // 力平衡项（软约束，使用较大权重）
    // 权重：W_FORCE_BAL / n_sample，确保合力接近零
    const double W_FORCE_BAL = 200.0;  // 进一步增大权重，确保力平衡
    auto force_balance_obj = std::make_shared<ForceBalanceObjective>(
        N_COEFF, N_SAMPLE, MASS, GRAVITY, W_FORCE_BAL / N_SAMPLE);
    nlp->addObjective(force_balance_obj);
    
    // 力矩平衡项（软约束，使用较大权重）
    // 权重：W_TORQUE_BAL / n_sample，确保合矩接近零
    const double W_TORQUE_BAL = 200.0;  // 进一步增大权重，确保力矩平衡
    auto torque_balance_obj = std::make_shared<TorqueBalanceObjective>(
        N_COEFF, N_SAMPLE, W_TORQUE_BAL / N_SAMPLE);
    nlp->addObjective(torque_balance_obj);
    
    // 摩擦锥违反项（软约束，使用较大权重）
    // 权重：W_FRIC_CONE / n_sample，确保摩擦锥满足
    const double MU_FRICTION = 0.6;
    const double W_FRIC_CONE = 500.0;  // 适当权重，确保摩擦锥满足但不至于过强
    auto friction_cone_obj = std::make_shared<FrictionConeObjective>(
        N_COEFF, N_SAMPLE, terrain, MU_FRICTION, W_FRIC_CONE / N_SAMPLE);
    nlp->addObjective(friction_cone_obj);
    
    // 测试：在第一次调用eval_f之前检查初始猜测（仅在DEBUG模式下）
    if (DEBUG_MODE) {
        std::cout << "\n=== Testing initial guess evaluation ===" << std::endl;
        std::vector<Ipopt::Number> test_x = x0;
        Ipopt::Number test_obj = 0.0;
        if (nlp->eval_f(static_cast<Ipopt::Index>(test_x.size()), test_x.data(), true, test_obj)) {
            std::cout << "  Initial objective value: " << test_obj << std::endl;
        } else {
            std::cerr << "  ERROR: Failed to evaluate objective at initial guess!" << std::endl;
            return -1;
        }
        
        // 测试约束
        Ipopt::Index n_vars_test = static_cast<Ipopt::Index>(test_x.size());
        Ipopt::Index n_cons_test = nlp->get_n_constraints();
        std::vector<Ipopt::Number> test_g(n_cons_test);
        if (nlp->eval_g(n_vars_test, test_x.data(), true, n_cons_test, test_g.data())) {
            std::cout << "  Initial constraint evaluation: OK" << std::endl;
            // 检查约束值范围
            double max_g = 0.0;
            for (Ipopt::Index i = 0; i < n_cons_test; ++i) {
                max_g = std::max(max_g, std::abs(test_g[i]));
            }
            std::cout << "  Max constraint value: " << max_g << std::endl;
        } else {
            std::cerr << "  ERROR: Failed to evaluate constraints at initial guess!" << std::endl;
            return -1;
        }
        std::cout << "=== Initial guess test complete ===\n" << std::endl;
    }
    
    // 创建Ipopt应用
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    
    // 设置选项
    app->Options()->SetStringValue("linear_solver", "mumps");
    app->Options()->SetIntegerValue("print_level", VERBOSE_IPOPT ? 12 : 5);  // 根据开关控制详细输出
    app->Options()->SetIntegerValue("max_iter", 1000);
    app->Options()->SetNumericValue("tol", 1e-4);
    app->Options()->SetNumericValue("acceptable_tol", 5e-4);
    app->Options()->SetIntegerValue("acceptable_iter", 10);
    app->Options()->SetStringValue("mu_strategy", "adaptive");  // 自适应障碍参数
    app->Options()->SetNumericValue("mu_init", 1e-1);  // 初始障碍参数
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");  // 使用L-BFGS近似
    
    // 启用自动缩放（归一化）以提高数值稳定性
    // 这对于不同量纲的变量（位置、力、角度等）混合的问题很重要
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");  // 基于梯度的自动缩放
    app->Options()->SetNumericValue("nlp_scaling_max_gradient", 100.0);  // 最大梯度缩放因子
    app->Options()->SetStringValue("nlp_scaling_obj_target_gradient", "yes");  // 目标函数梯度缩放
    app->Options()->SetNumericValue("obj_scaling_factor", 1.0);  // 目标函数缩放因子（1.0表示不缩放）
    
    // 变量和约束的自动缩放
    app->Options()->SetStringValue("nlp_scaling_direction", "gradient-based");  // 基于梯度的缩放方向
    
    // 初始化
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cerr << "Error during initialization!" << std::endl;
        return -1;
    }
    
    // 求解
    std::cout << "\n=== Solving NLP ===" << std::endl;
    
    // 记录开始时间
    auto start_time = std::chrono::high_resolution_clock::now();
    
    status = app->OptimizeTNLP(nlp);
    
    // 记录结束时间并计算耗时
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double elapsed_seconds = duration_ms.count() / 1000.0;
    
    // 输出时间统计
    std::cout << "\n=== Optimization Time Statistics ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    if (elapsed_seconds >= 60.0) {
        int minutes = static_cast<int>(elapsed_seconds / 60.0);
        double seconds = elapsed_seconds - minutes * 60.0;
        std::cout << "Total optimization time: " << minutes << "m " << seconds << "s" << std::endl;
    } else {
        std::cout << "Total optimization time: " << elapsed_seconds << " seconds" << std::endl;
    }
    std::cout << "Total optimization time: " << duration_ms.count() << " ms" << std::endl;
    std::cout << "Total optimization time: " << duration_us.count() << " μs" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    // 获取求解器状态码（SolverReturn，这是实际的状态码）
    // 注意：ApplicationReturnStatus和SolverReturn是不同的枚举类型
    Ipopt::SolverReturn solver_status = nlp->getSolverStatus();
    int solver_status_int = static_cast<int>(solver_status);
    int app_status_int = static_cast<int>(status);
    
    // 使用SolverReturn状态码进行判断（这是Ipopt内部实际使用的状态码）
    if (solver_status == Ipopt::Solve_Succeeded || solver_status == Ipopt::Solved_To_Acceptable_Level) {
        std::cout << "\n*** The problem solved!" << std::endl;
    } else {
        std::cout << "\n*** The problem FAILED!" << std::endl;
        std::cout << "Solver status code: " << solver_status_int << " (SolverReturn)" << std::endl;
        if (DEBUG_MODE) {
            std::cout << "Application status code: " << app_status_int << " (ApplicationReturnStatus)" << std::endl;
        }
        
        // 输出状态码含义（使用SolverReturn状态码）
        if (DEBUG_MODE) {
            switch (solver_status) {
                case Ipopt::Solve_Succeeded:
                    std::cout << "  Meaning: Solve succeeded" << std::endl;
                    break;
                case Ipopt::Solved_To_Acceptable_Level:
                    std::cout << "  Meaning: Solved to acceptable level" << std::endl;
                    break;
                case Ipopt::Invalid_Number_Detected:
                    std::cout << "  Meaning: Invalid number detected in problem (NaN or Inf)" << std::endl;
                    std::cout << "  This usually means:" << std::endl;
                    std::cout << "    - Division by zero or very small numbers" << std::endl;
                    std::cout << "    - Square root of negative numbers" << std::endl;
                    std::cout << "    - Numerical overflow in calculations" << std::endl;
                    std::cout << "    - Check terrain normal calculations and path tangent computations" << std::endl;
                    break;
                case Ipopt::Insufficient_Memory:
                    std::cout << "  Meaning: Insufficient memory" << std::endl;
                    break;
                case Ipopt::Restoration_Failed:
                    std::cout << "  Meaning: Restoration phase failed" << std::endl;
                    break;
                case Ipopt::Search_Direction_Becomes_Too_Small:
                    std::cout << "  Meaning: Search direction becomes too small" << std::endl;
                    std::cout << "  This usually means:" << std::endl;
                    std::cout << "    - Constraints may be too strict or conflicting" << std::endl;
                    std::cout << "    - Initial guess may be far from feasible region" << std::endl;
                    std::cout << "    - Problem may be ill-conditioned" << std::endl;
                    break;
                case Ipopt::Maximum_Iterations_Exceeded:
                    std::cout << "  Meaning: Maximum iterations exceeded" << std::endl;
                    break;
                case Ipopt::Diverging_Iterates:
                    std::cout << "  Meaning: Iterates are diverging" << std::endl;
                    break;
                case 5:  // Converged_To_Infeasible_Point (局部不可行)
                    std::cout << "  Meaning: Converged to a point of local infeasibility" << std::endl;
                    std::cout << "  This usually means:" << std::endl;
                    std::cout << "    - Problem may be infeasible" << std::endl;
                    std::cout << "    - Constraints may be too strict or conflicting" << std::endl;
                    std::cout << "    - Dual infeasibility is large (check Ipopt output above)" << std::endl;
                    std::cout << "    - Try relaxing constraints or adjusting weights" << std::endl;
                    break;
                default:
                    std::cout << "  Meaning: Unknown error (status " << status << ")" << std::endl;
                    std::cout << "  Common status codes:" << std::endl;
                    std::cout << "    0 = Solve_Succeeded" << std::endl;
                    std::cout << "    1 = Solved_To_Acceptable_Level" << std::endl;
                    std::cout << "    2 = Invalid_Number_Detected" << std::endl;
                    std::cout << "    3 = Search_Direction_Becomes_Too_Small" << std::endl;
                    std::cout << "    4 = Maximum_Iterations_Exceeded" << std::endl;
                    std::cout << "    5 = Converged_To_Infeasible_Point" << std::endl;
                    break;
            }
        } else {
            // 简化错误信息（按优先级判断，使用SolverReturn状态码）
            // 状态5：局部不可行（优先级最高，因为这是最常见的问题）
            if (solver_status_int == 5) {
                std::cout << "  Error: Converged to a point of local infeasibility" << std::endl;
                std::cout << "  Problem may be infeasible. Analysis:" << std::endl;
                std::cout << "    - Dual infeasibility: Check Ipopt output above (should be small)" << std::endl;
                std::cout << "    - Constraint violation: Check Ipopt output above (should be near zero)" << std::endl;
                std::cout << "  This means:" << std::endl;
                std::cout << "    - Primal constraints are satisfied, but dual problem is infeasible" << std::endl;
                std::cout << "    - Some Lagrange multipliers are very large" << std::endl;
                std::cout << "  Possible causes:" << std::endl;
                std::cout << "    1. Contact force objectives (W_CONTACT_STIFF, W_FRIC_CONE) may be too large" << std::endl;
                std::cout << "    2. Force/torque balance objectives may be too large" << std::endl;
                std::cout << "    3. Constraints may be conflicting (wheel-terrain distance vs force balance)" << std::endl;
                std::cout << "    4. Initial guess may be too far from feasible region" << std::endl;
                std::cout << "  Suggestions:" << std::endl;
                std::cout << "    - Current weights: W_CONTACT_STIFF=10, W_FRIC_CONE=10, W_FORCE_BAL=1, W_TORQUE_BAL=1" << std::endl;
                std::cout << "    - If still failing, try further reducing weights or removing some constraints" << std::endl;
                std::cout << "    - Check if wheel-terrain distance constraints are too strict" << std::endl;
                std::cout << "    - Verify initial contact forces satisfy force balance" << std::endl;
            } else if (solver_status_int == 2 || solver_status == Ipopt::Invalid_Number_Detected) {
                std::cout << "  Error: Invalid number detected (NaN or Inf)" << std::endl;
                std::cout << "  Check error messages above for which variable/constraint/objective failed" << std::endl;
            } else if (solver_status_int == 4 || solver_status == Ipopt::Search_Direction_Becomes_Too_Small) {
                std::cout << "  Error: Search direction too small" << std::endl;
            } else {
                std::cout << "  Error: Optimization failed (status " << solver_status_int << ")" << std::endl;
            }
        }
        
        // 即使失败，也尝试输出当前解的信息
        if (nlp->hasSolution()) {
            std::cout << "\nNote: A solution was found but may not be optimal." << std::endl;
        } else {
            std::cout << "\nWarning: No valid solution available." << std::endl;
        }
    }
    
    // 无论成功或失败，都尝试可视化轨迹（如果有解的话）
    if (nlp->hasSolution()) {
        // 提取轨迹系数并输出到终端（仅在DEBUG模式下详细输出）
        if (DEBUG_MODE) {
            std::cout << "\n=== Trajectory Coefficients ===" << std::endl;
        } else {
            std::cout << "\n=== Extracting Solution ===" << std::endl;
        }
        std::vector<double> cbx, cby, cbz, cpsi;
        std::vector<double> clx, cly, clz;
        std::vector<double> crx, cry, crz;
        
        nlp->extractCoefficients(nullptr, cbx, cby, cbz, cpsi, clx, cly, clz, crx, cry, crz);
        
        std::cout << "\nBase trajectory coefficients:" << std::endl;
        std::cout << "  cbx: ";
        for (size_t i = 0; i < cbx.size(); ++i) {
            std::cout << cbx[i];
            if (i < cbx.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "  cby: ";
        for (size_t i = 0; i < cby.size(); ++i) {
            std::cout << cby[i];
            if (i < cby.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "  cbz: ";
        for (size_t i = 0; i < cbz.size(); ++i) {
            std::cout << cbz[i];
            if (i < cbz.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "  cpsi: ";
        for (size_t i = 0; i < cpsi.size(); ++i) {
            std::cout << cpsi[i];
            if (i < cpsi.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        
        if (DEBUG_MODE) {
            std::cout << "\nLeft wheel trajectory coefficients:" << std::endl;
            std::cout << "  clx: ";
            for (size_t i = 0; i < clx.size(); ++i) {
                std::cout << clx[i];
                if (i < clx.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "  cly: ";
            for (size_t i = 0; i < cly.size(); ++i) {
                std::cout << cly[i];
                if (i < cly.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "  clz: ";
            for (size_t i = 0; i < clz.size(); ++i) {
                std::cout << clz[i];
                if (i < clz.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            
            std::cout << "\nRight wheel trajectory coefficients:" << std::endl;
            std::cout << "  crx: ";
            for (size_t i = 0; i < crx.size(); ++i) {
                std::cout << crx[i];
                if (i < crx.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "  cry: ";
            for (size_t i = 0; i < cry.size(); ++i) {
                std::cout << cry[i];
                if (i < cry.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
            std::cout << "  crz: ";
            for (size_t i = 0; i < crz.size(); ++i) {
                std::cout << crz[i];
                if (i < crz.size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        
        // 输出多项式表达式（仅在DEBUG模式下）
        if (DEBUG_MODE) {
            std::cout << "\n=== Polynomial Expressions (s ∈ [0, 1]) ===" << std::endl;
            std::cout << "\nBase X(s) = ";
            for (size_t i = 0; i < cbx.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << cbx[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Base Y(s) = ";
            for (size_t i = 0; i < cby.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << cby[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Base Z(s) = ";
            for (size_t i = 0; i < cbz.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << cbz[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Base Psi(s) = ";
            for (size_t i = 0; i < cpsi.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << cpsi[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "\nLeft Wheel X(s) = ";
            for (size_t i = 0; i < clx.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << clx[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Left Wheel Y(s) = ";
            for (size_t i = 0; i < cly.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << cly[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Left Wheel Z(s) = ";
            for (size_t i = 0; i < clz.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << clz[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "\nRight Wheel X(s) = ";
            for (size_t i = 0; i < crx.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << crx[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Right Wheel Y(s) = ";
            for (size_t i = 0; i < cry.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << cry[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
            
            std::cout << "Right Wheel Z(s) = ";
            for (size_t i = 0; i < crz.size(); ++i) {
                if (i > 0) std::cout << " + ";
                std::cout << crz[i];
                if (i > 0) std::cout << "*s^" << i;
            }
            std::cout << std::endl;
        }
        
        // 提取轨迹和接触力
        if (DEBUG_MODE) {
            std::cout << "\n=== Visualizing Trajectory ===" << std::endl;
        }
        
        std::vector<double> Xb, Yb, Zb, Psib;
        std::vector<double> Xl, Yl, Zl;
        std::vector<double> Xr, Yr, Zr;
        
        nlp->extractTrajectory(nullptr, Xb, Yb, Zb, Psib, Xl, Yl, Zl, Xr, Yr, Zr);
        
        TrajectoryVisualizer::visualize(*terrain, Xb, Yb, Zb, Psib,
                                       Xl, Yl, Zl, Xr, Yr, Zr,
                                       "/tmp/rbf_terrain_data.txt",
                                       VISUALIZATION_AUTO_CLOSE, VISUALIZATION_CLOSE_TIME);
        if (DEBUG_MODE) {
            std::cout << "=== Trajectory Visualization Complete ===\n" << std::endl;
        }
        
        // 提取接触力并可视化
        if (DEBUG_MODE) {
            std::cout << "\n=== Visualizing Forces ===" << std::endl;
        }
        std::vector<double> fLx, fLy, fLz;
        std::vector<double> fRx, fRy, fRz;
        
        nlp->extractContactForces(nullptr, fLx, fLy, fLz, fRx, fRy, fRz);
        
        const double MU_FRICTION = 0.6;
        
        ForceDataSaver::visualize(*terrain, fLx, fLy, fLz, fRx, fRy, fRz,
                                 Xb, Yb, Zb, Xl, Yl, Zl, Xr, Yr, Zr,
                                 MASS, GRAVITY, MU_FRICTION,
                                 VISUALIZATION_AUTO_CLOSE, VISUALIZATION_CLOSE_TIME);
        if (DEBUG_MODE) {
            std::cout << "=== Force Visualization Complete ===\n" << std::endl;
        }
    } else {
        std::cout << "\n*** No solution available for visualization" << std::endl;
    }
    
    return 0;
}

