#ifndef RRT_STAR_REGION_HPP
#define RRT_STAR_REGION_HPP

#include <vector>
#include <string>

namespace rrt_star_region {

struct SegDataView {
    int nx = 0, ny = 0;
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    std::vector<int> labels;
    int get_label(double x, double y) const;
};

bool load_segmentation(const std::string& path, SegDataView& out);

struct Params {
    double start_x = 0.6, start_y = 1.2;
    double goal_x = 2.4, goal_y = -3.0;
    double lift_cost_weight = 5.0;
    int max_iter = 2000;
    double step_length = 0.15;
    double goal_radius = 0.2;
    double rewire_radius = 0.5;
    int edge_samples = 20;
    int grid_buckets = 64;  // spatial index: nbx = nby = grid_buckets
};

// Returns path as (x,y) pairs; empty if not found.
std::vector<std::pair<double, double>> run(
    const std::string& data_file,
    const Params& params,
    double* out_time_ms = nullptr);

// Print waypoints where path switches rollable <-> lift (uses seg.get_label along path).
void print_switch_waypoints(const std::vector<std::pair<double, double>>& path,
                            const SegDataView& seg);

// Write segmentation + path + start/goal to PNG and optionally open (no Python).
bool write_visualization_png(const std::string& out_path,
                             const SegDataView& seg,
                             const std::vector<std::pair<double, double>>& path,
                             double start_x, double start_y, double goal_x, double goal_y,
                             bool open_after = true);

// Save waypoints with region label (x y label per line; label 0=lift, 1=rollable) for fine planning.
bool save_waypoints_segmented(const std::string& out_path,
                              const std::vector<std::pair<double, double>>& path,
                              const SegDataView& seg);

}  // namespace rrt_star_region

#endif
