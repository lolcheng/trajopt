#include "segmentation_reader.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

bool SegmentationReader::load(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    if (!std::getline(f, line)) return false;
    std::istringstream ss(line);
    if (!(ss >> nx_ >> ny_ >> xmin_ >> xmax_ >> ymin_ >> ymax_)) return false;
    if (!std::getline(f, line)) return false;  // skip terrain line
    if (!std::getline(f, line)) return false;
    labels_.clear();
    std::istringstream ss2(line);
    int v;
    while (ss2 >> v) labels_.push_back(static_cast<double>(v));
    return static_cast<int>(labels_.size()) == nx_ * ny_;
}

double SegmentationReader::interp(double x, double y) const {
    if (nx_ < 2 || ny_ < 2 || labels_.empty()) return 0.0;
    double xi = (x - xmin_) / (xmax_ - xmin_ + 1e-12) * (nx_ - 1);
    double yj = (y - ymin_) / (ymax_ - ymin_ + 1e-12) * (ny_ - 1);
    xi = std::max(0.0, std::min(static_cast<double>(nx_ - 1), xi));
    yj = std::max(0.0, std::min(static_cast<double>(ny_ - 1), yj));
    int i0 = static_cast<int>(std::floor(xi));
    int j0 = static_cast<int>(std::floor(yj));
    i0 = std::max(0, std::min(nx_ - 2, i0));
    j0 = std::max(0, std::min(ny_ - 2, j0));
    int i1 = i0 + 1, j1 = j0 + 1;
    double tx = xi - i0, ty = yj - j0;
    double L00 = labels_[j0 * nx_ + i0];
    double L10 = labels_[j0 * nx_ + i1];
    double L01 = labels_[j1 * nx_ + i0];
    double L11 = labels_[j1 * nx_ + i1];
    return (1 - tx) * (1 - ty) * L00 + tx * (1 - ty) * L10 + (1 - tx) * ty * L01 + tx * ty * L11;
}

void SegmentationReader::interpGrad(double x, double y, double& ddx, double& ddy) const {
    ddx = 0.0;
    ddy = 0.0;
    if (nx_ < 2 || ny_ < 2 || labels_.empty()) return;
    double dx = (xmax_ - xmin_) / (nx_ - 1);
    double dy = (ymax_ - ymin_) / (ny_ - 1);
    if (dx < 1e-12 || dy < 1e-12) return;
    double xi = (x - xmin_) / (xmax_ - xmin_ + 1e-12) * (nx_ - 1);
    double yj = (y - ymin_) / (ymax_ - ymin_ + 1e-12) * (ny_ - 1);
    xi = std::max(0.0, std::min(static_cast<double>(nx_ - 1), xi));
    yj = std::max(0.0, std::min(static_cast<double>(ny_ - 1), yj));
    int i0 = static_cast<int>(std::floor(xi));
    int j0 = static_cast<int>(std::floor(yj));
    i0 = std::max(0, std::min(nx_ - 2, i0));
    j0 = std::max(0, std::min(ny_ - 2, j0));
    int i1 = i0 + 1, j1 = j0 + 1;
    double tx = xi - i0, ty = yj - j0;
    double L00 = labels_[j0 * nx_ + i0];
    double L10 = labels_[j0 * nx_ + i1];
    double L01 = labels_[j1 * nx_ + i0];
    double L11 = labels_[j1 * nx_ + i1];
    ddx = ((1 - ty) * (L10 - L00) + ty * (L11 - L01)) / dx;
    ddy = ((1 - tx) * (L01 - L00) + tx * (L11 - L10)) / dy;
}
