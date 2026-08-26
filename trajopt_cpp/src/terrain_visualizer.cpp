#include "terrain_visualizer.hpp"
#include "rbf_terrain.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

// 获取RBF中心的范围
static void getTerrainBounds(const RBFTerrain& terrain, 
                             double& xmin, double& xmax,
                             double& ymin, double& ymax) {
    terrain.getBounds(xmin, xmax, ymin, ymax);
}

void TerrainVisualizer::visualize(const RBFTerrain& terrain, 
                                  int nx, int ny, double margin,
                                  bool auto_close, double close_time) {
    double xmin, xmax, ymin, ymax;
    getTerrainBounds(terrain, xmin, xmax, ymin, ymax);
    
    xmin -= margin;
    xmax += margin;
    ymin -= margin;
    ymax += margin;
    
    visualize(terrain, xmin, xmax, ymin, ymax, nx, ny, auto_close, close_time);
}

void TerrainVisualizer::visualize(const RBFTerrain& terrain,
                                  double xmin, double xmax,
                                  double ymin, double ymax,
                                  int nx, int ny,
                                  bool auto_close, double close_time) {
    std::cout << "Generating terrain visualization data..." << std::endl;
    
    // 创建网格
    std::vector<double> xs(nx), ys(ny);
    for (int i = 0; i < nx; ++i) {
        xs[i] = xmin + (xmax - xmin) * i / (nx - 1);
    }
    for (int i = 0; i < ny; ++i) {
        ys[i] = ymin + (ymax - ymin) * i / (ny - 1);
    }
    
    // 计算地形高度
    std::vector<std::vector<double>> Z(ny, std::vector<double>(nx));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Z[j][i] = terrain.height(xs[i], ys[j]);
        }
        if ((j + 1) % 20 == 0) {
            std::cout << "  Progress: " << (j + 1) << "/" << ny << " rows" << std::endl;
        }
    }
    
    // 输出数据到文件
    std::string data_file = "/tmp/rbf_terrain_data.txt";
    std::ofstream out(data_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file " << data_file << " for writing!" << std::endl;
        return;
    }
    
    out << std::fixed << std::setprecision(8);
    out << nx << " " << ny << std::endl;
    out << xmin << " " << xmax << " " << ymin << " " << ymax << std::endl;
    
    // 输出x坐标
    for (int i = 0; i < nx; ++i) {
        out << xs[i];
        if (i < nx - 1) out << " ";
    }
    out << std::endl;
    
    // 输出y坐标
    for (int i = 0; i < ny; ++i) {
        out << ys[i];
        if (i < ny - 1) out << " ";
    }
    out << std::endl;
    
    // 输出高度数据
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            out << Z[j][i];
            if (i < nx - 1) out << " ";
        }
        out << std::endl;
    }
    out.close();
    
    std::cout << "Terrain data saved to: " << data_file << std::endl;
    std::cout << "\nTo visualize the terrain, run:" << std::endl;
    std::cout << "  python3 visualize_terrain.py " << data_file << std::endl;
    std::cout << "\nOr from the project root:" << std::endl;
    std::cout << "  python3 trajopt_cpp/visualize_terrain.py " << data_file << std::endl;
    
    // 尝试自动运行可视化
    std::cout << "\nAttempting to run visualization automatically..." << std::endl;
    
    // 尝试多个可能的脚本路径
    std::vector<std::string> possible_paths = {
        "/home/yizhe/trajopt/trajopt_cpp/visualize_terrain.py",
        "../visualize_terrain.py",
        "../../trajopt_cpp/visualize_terrain.py",
        "visualize_terrain.py"
    };
    
    bool success = false;
    for (const auto& script_path : possible_paths) {
        std::string cmd = "python3 " + script_path + " " + data_file;
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
        std::cout << "  python3 trajopt_cpp/visualize_terrain.py " << data_file << std::endl;
    }
}

