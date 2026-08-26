#include "trajectory_visualizer.hpp"
#include "rbf_terrain.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>

void TrajectoryVisualizer::visualize(const RBFTerrain& terrain,
                                     const std::vector<double>& Xb,
                                     const std::vector<double>& Yb,
                                     const std::vector<double>& Zb,
                                     const std::vector<double>& Psib,
                                     const std::vector<double>& Xl,
                                     const std::vector<double>& Yl,
                                     const std::vector<double>& Zl,
                                     const std::vector<double>& Xr,
                                     const std::vector<double>& Yr,
                                     const std::vector<double>& Zr,
                                     const std::string& terrain_data_file,
                                     bool auto_close, double close_time) {
    int n_points = static_cast<int>(Xb.size());
    if (n_points == 0) {
        std::cerr << "Error: Empty trajectory data!" << std::endl;
        return;
    }
    
    std::cout << "Saving trajectory data..." << std::endl;
    
    // 保存轨迹数据到文件
    std::string traj_file = "/tmp/trajectory_data.txt";
    std::ofstream out(traj_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file " << traj_file << " for writing!" << std::endl;
        return;
    }
    
    out << std::fixed << std::setprecision(8);
    out << n_points << std::endl;
    
    // 计算轮子离地距离并输出轨迹数据
    std::vector<double> Hl(n_points), Hr(n_points);  // 地形高度
    std::vector<double> airborne_l(n_points), airborne_r(n_points);  // 离地距离
    
    // 计算地形高度和离地距离
    for (int i = 0; i < n_points; ++i) {
        Hl[i] = terrain.height(Xl[i], Yl[i]);
        Hr[i] = terrain.height(Xr[i], Yr[i]);
        airborne_l[i] = Zl[i] - Hl[i];  // 左轮离地距离（正数表示离地，负数表示下陷）
        airborne_r[i] = Zr[i] - Hr[i];  // 右轮离地距离
    }
    
    // 输出轨迹数据（包含离地距离）
    for (int i = 0; i < n_points; ++i) {
        out << Xb[i] << " " << Yb[i] << " " << Zb[i] << " " << Psib[i] << " "
            << Xl[i] << " " << Yl[i] << " " << Zl[i] << " "
            << Xr[i] << " " << Yr[i] << " " << Zr[i] << " "
            << airborne_l[i] << " " << airborne_r[i] << std::endl;
    }
    out.close();
    
    std::cout << "Trajectory data saved to: " << traj_file << std::endl;
    
    // 确定可视化范围（基于轨迹和地形）
    double xmin = *std::min_element(Xb.begin(), Xb.end());
    double xmax = *std::max_element(Xb.begin(), Xb.end());
    double ymin = *std::min_element(Yb.begin(), Yb.end());
    double ymax = *std::max_element(Yb.begin(), Yb.end());
    
    // 扩展范围以包含轮子
    xmin = std::min({xmin, *std::min_element(Xl.begin(), Xl.end()),
                     *std::min_element(Xr.begin(), Xr.end())});
    xmax = std::max({xmax, *std::max_element(Xl.begin(), Xl.end()),
                     *std::max_element(Xr.begin(), Xr.end())});
    ymin = std::min({ymin, *std::min_element(Yl.begin(), Yl.end()),
                     *std::min_element(Yr.begin(), Yr.end())});
    ymax = std::max({ymax, *std::max_element(Yl.begin(), Yl.end()),
                     *std::max_element(Yr.begin(), Yr.end())});
    
    // 添加边界
    double margin = 0.5;
    xmin -= margin;
    xmax += margin;
    ymin -= margin;
    ymax += margin;
    
    // 获取地形范围
    double terrain_xmin, terrain_xmax, terrain_ymin, terrain_ymax;
    terrain.getBounds(terrain_xmin, terrain_xmax, terrain_ymin, terrain_ymax);
    
    // 合并范围
    xmin = std::min(xmin, terrain_xmin);
    xmax = std::max(xmax, terrain_xmax);
    ymin = std::min(ymin, terrain_ymin);
    ymax = std::max(ymax, terrain_ymax);
    
    // 生成地形数据（如果还没有）
    std::string terrain_file = terrain_data_file.empty() ? "/tmp/rbf_terrain_data.txt" : terrain_data_file;
    bool need_terrain = true;
    if (!terrain_data_file.empty()) {
        std::ifstream test(terrain_file);
        if (test.good()) {
            need_terrain = false;
            test.close();
        }
    }
    
    if (need_terrain) {
        std::cout << "Generating terrain data for visualization..." << std::endl;
        int nx = 161, ny = 161;
        std::vector<double> xs(nx), ys(ny);
        for (int i = 0; i < nx; ++i) {
            xs[i] = xmin + (xmax - xmin) * i / (nx - 1);
        }
        for (int i = 0; i < ny; ++i) {
            ys[i] = ymin + (ymax - ymin) * i / (ny - 1);
        }
        
        std::ofstream terrain_out(terrain_file);
        terrain_out << std::fixed << std::setprecision(8);
        terrain_out << nx << " " << ny << std::endl;
        terrain_out << xmin << " " << xmax << " " << ymin << " " << ymax << std::endl;
        
        for (int i = 0; i < nx; ++i) {
            terrain_out << xs[i];
            if (i < nx - 1) terrain_out << " ";
        }
        terrain_out << std::endl;
        
        for (int i = 0; i < ny; ++i) {
            terrain_out << ys[i];
            if (i < ny - 1) terrain_out << " ";
        }
        terrain_out << std::endl;
        
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                terrain_out << terrain.height(xs[i], ys[j]);
                if (i < nx - 1) terrain_out << " ";
            }
            terrain_out << std::endl;
        }
        terrain_out.close();
    }
    
    std::cout << "\nTo visualize the trajectory, run:" << std::endl;
    std::cout << "  python3 trajopt_cpp/visualize_trajectory.py " << traj_file << " " << terrain_file << std::endl;
    
    // 尝试自动运行可视化
    std::cout << "\nAttempting to run visualization automatically..." << std::endl;
    
    std::vector<std::string> possible_paths = {
        "/home/yizhe/trajopt/trajopt_cpp/visualize_trajectory.py",
        "../visualize_trajectory.py",
        "../../trajopt_cpp/visualize_trajectory.py",
        "visualize_trajectory.py"
    };
    
    bool success = false;
    for (const auto& script_path : possible_paths) {
        std::string cmd = "python3 " + script_path + " " + traj_file + " " + terrain_file;
        if (auto_close) {
            cmd += " --auto-close " + std::to_string(close_time);
        }
        cmd += " 2>&1";
        int ret = system(cmd.c_str());
        if (ret == 0) {
            success = true;
            break;
        }
    }
    
    if (!success) {
        std::cout << "\nAutomatic visualization failed. Please run manually:" << std::endl;
        std::cout << "  python3 trajopt_cpp/visualize_trajectory.py " << traj_file << " " << terrain_file << std::endl;
    }
}

