#include "waypoint_reader.hpp"
#include <fstream>
#include <sstream>

bool WaypointReader::readFirstRollableSegment(const std::string& filepath,
                                              std::vector<double>& x,
                                              std::vector<double>& y) {
    x.clear();
    y.clear();
    std::ifstream f(filepath);
    if (!f.is_open()) return false;

    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        double xi, yi;
        int seg_type;
        if (!(ss >> xi >> yi >> seg_type)) continue;
        if (seg_type != 1) {
            if (x.empty()) continue;
            break;
        }
        x.push_back(xi);
        y.push_back(yi);
    }
    return !x.empty();
}
