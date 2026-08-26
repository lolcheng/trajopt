/**
 * 地形摩擦锥分割：使用 trajopt_cpp 生成的 RBF 文件，GPU 判断每点摩擦锥，调用 Python 可视化。
 * 用法: ./terrain_friction_segment [rbf.json] [output_base]
 * 默认 RBF 路径与 trajopt_cpp/main.cpp:64 一致。
 */

#include "json_reader.hpp"
#include "rbf_friction_cuda.hpp"
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <filesystem>

const bool VISUALIZATION_AUTO_CLOSE = false;
const double VISUALIZATION_CLOSE_TIME = 5.0;

namespace {

/** 在相同网格上计算 RBF 地形高度，行优先，与 GPU 分割网格一致 */
std::vector<double> computeTerrainHeightGrid(
    const std::vector<double>& cx, const std::vector<double>& cy,
    const std::vector<double>& w, double sigma,
    double xmin, double xmax, double ymin, double ymax, int nx, int ny) {
    const double sigma2 = sigma * sigma;
    const int n_centers = static_cast<int>(cx.size());
    std::vector<double> out(static_cast<size_t>(nx) * ny);
    for (int j = 0; j < ny; ++j) {
        double y = (ny > 1) ? (ymin + (ymax - ymin) * static_cast<double>(j) / (ny - 1)) : ymin;
        for (int i = 0; i < nx; ++i) {
            double x = (nx > 1) ? (xmin + (xmax - xmin) * static_cast<double>(i) / (nx - 1)) : xmin;
            double h = 0.0;
            for (int k = 0; k < n_centers; ++k) {
                double dx = x - cx[static_cast<size_t>(k)];
                double dy = y - cy[static_cast<size_t>(k)];
                double r2 = dx * dx + dy * dy;
                h += w[static_cast<size_t>(k)] * std::exp(-0.5 * r2 / sigma2);
            }
            out[static_cast<size_t>(j) * nx + i] = h;
        }
    }
    return out;
}

/** 写入可视化数据文件供 Python 脚本读取：第1行 nx ny xmin xmax ymin ymax，第2行地形高度，第3行标签 */
bool writeVisualizationData(const std::string& path,
                            int nx, int ny,
                            double xmin, double xmax, double ymin, double ymax,
                            const std::vector<double>& terrain_height,
                            const std::vector<int>& labels) {
    const size_t n = static_cast<size_t>(nx) * ny;
    if (terrain_height.size() != n || labels.size() != n) return false;
    std::ofstream out(path);
    if (!out) return false;
    out.precision(16);
    out << nx << " " << ny << " " << xmin << " " << xmax << " " << ymin << " " << ymax << "\n";
    for (size_t i = 0; i < n; ++i) {
        if (i) out << " ";
        out << terrain_height[i];
    }
    out << "\n";
    for (size_t i = 0; i < n; ++i) {
        if (i) out << " ";
        out << labels[i];
    }
    out << "\n";
    return !!out;
}

/** 调用 Python 可视化脚本（与 trajopt_cpp TrajectoryVisualizer 相同方式） */
void runVisualizationScript(const std::string& data_file, bool auto_close, double close_time) {
    std::vector<std::string> possible_paths = {
        "/home/yizhe/trajopt/terrain_friction_segment/scripts/visualize_terrain_segmentation.py",
        "scripts/visualize_terrain_segmentation.py",
        "../scripts/visualize_terrain_segmentation.py",
        "../../terrain_friction_segment/scripts/visualize_terrain_segmentation.py",
    };
    std::cout << "\nAttempting to run visualization..." << std::endl;
    bool success = false;
    for (const auto& script_path : possible_paths) {
        std::string cmd = "python3 " + script_path + " " + data_file;
        if (auto_close) cmd += " --auto-close " + std::to_string(close_time);
        cmd += " 2>&1";
        int ret = std::system(cmd.c_str());
        if (ret == 0) { success = true; break; }
    }
    if (!success) {
        std::cout << "\nAutomatic visualization failed. Run manually:" << std::endl;
        std::cout << "  python3 scripts/visualize_terrain_segmentation.py " << data_file << std::endl;
    }
}

} // namespace

int main(int argc, char** argv) {
    using Clock = std::chrono::high_resolution_clock;
    using Ms = std::chrono::duration<double, std::milli>;
    auto t_total_start = Clock::now();

    // 默认 RBF 与 trajopt_cpp/main.cpp:64 一致；可选第4参数为摩擦系数 mu（更小=更严，台阶/陡坡更易判为需抬腿）
    const std::string default_rbf = "/home/yizhe/trajopt/result/RBF/1753771011_476223707.json";
    std::string json_path = (argc > 1) ? argv[1] : default_rbf;
    const std::string default_out = "/home/yizhe/trajopt/test/terrain_res";
    const std::string out_dir = (argc > 2) ? argv[2] : default_out;
    std::filesystem::create_directories(out_dir);
    std::string out_base = out_dir + "/segmentation";
    double mu = 0.3;  // 默认 0.3（更严，便于识别台阶等陡坡）；可传 0.4/0.6
    if (argc > 3) {
        try { mu = std::stod(argv[3]); } catch (...) {}
        mu = std::max(0.1, std::min(1.0, mu));
    }

    std::vector<double> centers_x, centers_y, weights;
    double sigma;
    auto t_load_start = Clock::now();
    std::cout << "Reading RBF from: " << json_path << std::endl;
    if (!JSONReader::readRBFParams(json_path, centers_x, centers_y, weights, sigma)) {
        std::cerr << "Failed to load RBF. Try: " << argv[0] << " /path/to/trajopt_cpp/rbf.json" << std::endl;
        return 1;
    }
    auto t_load_end = Clock::now();

    double xmin = *std::min_element(centers_x.begin(), centers_x.end());
    double xmax = *std::max_element(centers_x.begin(), centers_x.end());
    double ymin = *std::min_element(centers_y.begin(), centers_y.end());
    double ymax = *std::max_element(centers_y.begin(), centers_y.end());
    const double margin = 0.3;
    xmin -= margin;
    xmax += margin;
    ymin -= margin;
    ymax += margin;

    rbf_friction_cuda::Params params;
    params.xmin = xmin;
    params.xmax = xmax;
    params.ymin = ymin;
    params.ymax = ymax;
    params.nx = 512;
    params.ny = 512;
    params.mu = mu;

    const double nz_min = 1.0 / std::sqrt(1.0 + mu * mu);
    std::cout << "Friction cone: mu=" << mu << " => n_z_min=" << std::fixed << std::setprecision(3) << nz_min
              << " (stricter for steps/ramps)" << std::endl;
    std::vector<int> labels;
    auto t_gpu_start = Clock::now();
    std::cout << "Running GPU friction cone segmentation (" << params.nx << "x" << params.ny << ")..." << std::endl;
    if (!rbf_friction_cuda::run(centers_x, centers_y, weights, sigma, params, labels)) {
        std::cerr << "GPU segmentation failed." << std::endl;
        return 1;
    }
    auto t_gpu_end = Clock::now();

    int rollable = 0;
    for (int v : labels) if (v) rollable++;
    std::cout << "Rollable: " << rollable << " / " << labels.size() << std::endl;

    auto t_io_start = Clock::now();
    // Copy RBF JSON to output dir for later use in fine planning
    std::string rbf_dest = out_dir + "/rbf.json";
    try {
        std::filesystem::copy_file(json_path, rbf_dest, std::filesystem::copy_options::overwrite_existing);
        std::cout << "Saved RBF copy: " << rbf_dest << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Failed to copy RBF to " << rbf_dest << ": " << e.what() << std::endl;
    }

    std::string txt_path = out_base + ".txt";
    if (!rbf_friction_cuda::saveSegmentation(txt_path,
                                             params.nx, params.ny,
                                             params.xmin, params.xmax, params.ymin, params.ymax,
                                             labels)) {
        std::cerr << "Failed to write " << txt_path << std::endl;
        return 1;
    }
    std::cout << "Saved " << txt_path << std::endl;

    std::string png_path = out_base + ".png";
    if (rbf_friction_cuda::saveSegmentationImage(png_path, params.nx, params.ny, labels, false))
        std::cout << "Saved " << png_path << " (green=rollable, red=lift leg)" << std::endl;
    else
        std::cerr << "Failed to write " << png_path << std::endl;

    auto t_terrain_start = Clock::now();
    std::vector<double> terrain_height = computeTerrainHeightGrid(
        centers_x, centers_y, weights, sigma,
        params.xmin, params.xmax, params.ymin, params.ymax,
        params.nx, params.ny);
    auto t_terrain_end = Clock::now();

    std::string comparison_path = out_base + "_comparison.png";
    if (rbf_friction_cuda::saveTerrainSegmentationComparison(
            comparison_path, params.nx, params.ny, terrain_height, labels, false))
        std::cout << "Saved " << comparison_path << std::endl;

    const std::string viz_data_file = out_dir + "/terrain_segmentation_data.txt";
    if (writeVisualizationData(viz_data_file,
                               params.nx, params.ny,
                               params.xmin, params.xmax, params.ymin, params.ymax,
                               terrain_height, labels)) {
        std::cout << "Visualization data saved to: " << viz_data_file << std::endl;
        std::cout << "\nTo visualize (terrain vs segmentation), run:" << std::endl;
        std::cout << "  python3 scripts/visualize_terrain_segmentation.py " << viz_data_file << std::endl;
        runVisualizationScript(viz_data_file, VISUALIZATION_AUTO_CLOSE, VISUALIZATION_CLOSE_TIME);
    }
    auto t_io_end = Clock::now();

    auto t_total_end = Clock::now();
    std::cout << "\n--- Run time (ms) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Load RBF:        " << Ms(t_load_end - t_load_start).count() << std::endl;
    std::cout << "  GPU segmentation:" << Ms(t_gpu_end - t_gpu_start).count() << std::endl;
    std::cout << "  Terrain grid:    " << Ms(t_terrain_end - t_terrain_start).count() << std::endl;
    std::cout << "  I/O + viz call:  " << Ms(t_io_end - t_io_start).count() << std::endl;
    std::cout << "  Total:           " << Ms(t_total_end - t_total_start).count() << std::endl;

    return 0;
}
