#include "force_data_saver.hpp"
#include "rbf_terrain.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>
#include <unistd.h>  // for getpid()

void ForceDataSaver::visualize(const RBFTerrain& terrain,
                               const std::vector<double>& fLx, const std::vector<double>& fLy,
                               const std::vector<double>& fLz,
                               const std::vector<double>& fRx, const std::vector<double>& fRy,
                               const std::vector<double>& fRz,
                               const std::vector<double>& Xb, const std::vector<double>& Yb,
                               const std::vector<double>& Zb,
                               const std::vector<double>& Xl, const std::vector<double>& Yl,
                               const std::vector<double>& Zl,
                               const std::vector<double>& Xr, const std::vector<double>& Yr,
                               const std::vector<double>& Zr,
                               double MASS, double GRAVITY, double MU_FRICTION,
                               bool auto_close, double close_time) {
    int n_sample = static_cast<int>(fLx.size());
    if (n_sample == 0) {
        std::cerr << "Error: Empty force data!" << std::endl;
        return;
    }
    
    // 使用临时文件（程序结束后会自动删除）
    std::string tmp_file = "/tmp/force_data_" + std::to_string(getpid()) + ".txt";
    
    std::ofstream out(tmp_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create temporary file for force data!" << std::endl;
        return;
    }
    
    out << std::fixed << std::setprecision(8);
    
    // 写入采样点数
    out << n_sample << std::endl;
    
    // 写入接触力
    for (int i = 0; i < n_sample; ++i) {
        out << fLx[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << fLy[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << fLz[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << fRx[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << fRy[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << fRz[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    
    // 写入轨迹数据
    for (int i = 0; i < n_sample; ++i) {
        out << Xb[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Yb[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Zb[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Xl[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Yl[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Zl[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Xr[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Yr[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << Zr[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    
    // 计算并写入地形法向量
    std::vector<double> nLx(n_sample), nLy(n_sample), nLz(n_sample);
    std::vector<double> nRx(n_sample), nRy(n_sample), nRz(n_sample);
    
    for (int i = 0; i < n_sample; ++i) {
        double dhx_l, dhy_l;
        terrain.gradient(Xl[i], Yl[i], dhx_l, dhy_l);
        double nLx_raw = -dhx_l;
        double nLy_raw = -dhy_l;
        double nLz_raw = 1.0;
        double nL_norm = std::sqrt(nLx_raw * nLx_raw + nLy_raw * nLy_raw + nLz_raw * nLz_raw + 1e-6);
        nLx[i] = nLx_raw / nL_norm;
        nLy[i] = nLy_raw / nL_norm;
        nLz[i] = nLz_raw / nL_norm;
        
        double dhx_r, dhy_r;
        terrain.gradient(Xr[i], Yr[i], dhx_r, dhy_r);
        double nRx_raw = -dhx_r;
        double nRy_raw = -dhy_r;
        double nRz_raw = 1.0;
        double nR_norm = std::sqrt(nRx_raw * nRx_raw + nRy_raw * nRy_raw + nRz_raw * nRz_raw + 1e-6);
        nRx[i] = nRx_raw / nR_norm;
        nRy[i] = nRy_raw / nR_norm;
        nRz[i] = nRz_raw / nR_norm;
    }
    
    for (int i = 0; i < n_sample; ++i) {
        out << nLx[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << nLy[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << nLz[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << nRx[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << nRy[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    for (int i = 0; i < n_sample; ++i) {
        out << nRz[i];
        if (i < n_sample - 1) out << " ";
    }
    out << std::endl;
    
    // 写入参数
    out << MASS << std::endl;
    out << GRAVITY << std::endl;
    out << MU_FRICTION << std::endl;
    
    out.close();
    
    // 调用Python可视化脚本
    std::vector<std::string> possible_paths = {
        "/home/yizhe/trajopt/trajopt_cpp/visualize_forces.py",
        "../visualize_forces.py",
        "../../trajopt_cpp/visualize_forces.py",
        "visualize_forces.py"
    };
    
    bool success = false;
    for (const auto& script_path : possible_paths) {
        std::string cmd = "python3 " + script_path + " " + tmp_file;
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
    
    // 删除临时文件
    std::remove(tmp_file.c_str());
    
    if (!success) {
        std::cout << "\nAutomatic force visualization failed. Please check the Python script path." << std::endl;
    }
}

