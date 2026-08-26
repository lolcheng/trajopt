#include "waypoint_reader.hpp"
#include <fstream>
#include <sstream>

namespace {
bool readFirstSegmentByType(const std::string& filepath, int desired_type,
                            std::vector<double>& x, std::vector<double>& y) {
    x.clear();
    y.clear();
    std::ifstream f(filepath);
    if (!f.is_open()) return false;

    std::string line;
    bool collecting = false;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        double xi, yi;
        int seg_type;
        if (!(ss >> xi >> yi >> seg_type)) continue;
        if (seg_type == desired_type) {
            collecting = true;
            x.push_back(xi);
            y.push_back(yi);
        } else if (collecting) {
            break;
        }
    }
    return !x.empty();
}
}  // namespace

bool WaypointReader::readFirstRollableSegment(const std::string& filepath,
                                              std::vector<double>& x,
                                              std::vector<double>& y) {
    return readFirstSegmentByType(filepath, 1, x, y);
}

bool WaypointReader::readFirstLiftLegSegment(const std::string& filepath,
                                             std::vector<double>& x,
                                             std::vector<double>& y) {
    return readFirstSegmentByType(filepath, 0, x, y);
}
