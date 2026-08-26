/**
 * RRT* with region cost (C++). Reads segmentation, plans path, writes path and waypoints to terrain_res.
 * Usage: rrt_star_region [data.txt] [start_x start_y goal_x goal_y] [path_out.txt]
 * Default: data/path/png under test/terrain_res/
 */

#include "rrt_star_region.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

int main(int argc, char** argv) {
    const std::string out_dir = "/home/yizhe/trajopt/test/terrain_res";
    std::filesystem::create_directories(out_dir);
    std::string data_file = out_dir + "/terrain_segmentation_data.txt";
    std::string path_file = out_dir + "/rrt_star_path.txt";
    std::string png_file = out_dir + "/rrt_star_path.png";
    std::string waypoints_file = out_dir + "/waypoints_segmented.txt";
    rrt_star_region::Params params;
    if (argc >= 2) data_file = argv[1];
    if (argc >= 6) {
        params.start_x = std::stod(argv[2]);
        params.start_y = std::stod(argv[3]);
        params.goal_x = std::stod(argv[4]);
        params.goal_y = std::stod(argv[5]);
    }
    if (argc >= 7) path_file = argv[6];

    double time_ms = 0;
    std::vector<std::pair<double, double>> path = rrt_star_region::run(data_file, params, &time_ms);

    if (path.empty()) {
        std::cerr << "No path found." << std::endl;
        return 1;
    }

    std::ofstream out(path_file);
    if (!out) {
        std::cerr << "Cannot write " << path_file << std::endl;
        return 1;
    }
    out << std::fixed << std::setprecision(6);
    for (const auto& p : path) out << p.first << " " << p.second << "\n";
    out.close();

    std::cout << "Path: " << path.size() << " waypoints, " << std::fixed << std::setprecision(1) << time_ms << " ms" << std::endl;
    std::cout << "Saved: " << path_file << std::endl;

    rrt_star_region::SegDataView seg;
    if (rrt_star_region::load_segmentation(data_file, seg)) {
        std::cout << "Switch waypoints (rollable <-> lift):" << std::endl;
        rrt_star_region::print_switch_waypoints(path, seg);
        if (rrt_star_region::save_waypoints_segmented(waypoints_file, path, seg)) {
            std::cout << "Waypoints (x y region): " << waypoints_file << std::endl;
        }
        if (rrt_star_region::write_visualization_png(png_file, seg, path,
                params.start_x, params.start_y, params.goal_x, params.goal_y, true)) {
            std::cout << "Visualization: " << png_file << std::endl;
        }
    }

    return 0;
}
