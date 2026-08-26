#ifndef WAYPOINT_READER_HPP
#define WAYPOINT_READER_HPP

#include <string>
#include <vector>

/**
 * 从 waypoints_segmented.txt 读取第一段 rollable 路径的 waypoints
 * 格式: 每行 "x y segment_type" (segment_type 1=rollable, 0=non-rollable)
 * 返回从起点开始连续 segment_type==1 的 (x, y) 序列
 */
class WaypointReader {
public:
    /**
     * 读取文件并提取第一段 rollable 的 waypoints
     * @param filepath 如 "terrain_res/waypoints_segmented.txt"
     * @param x 输出 x 坐标
     * @param y 输出 y 坐标
     * @return 是否成功且至少有一个点
     */
    static bool readFirstRollableSegment(const std::string& filepath,
                                         std::vector<double>& x,
                                         std::vector<double>& y);
};

#endif
