#ifndef SEGMENTATION_READER_HPP
#define SEGMENTATION_READER_HPP

#include <string>
#include <vector>

/**
 * 读取 terrain_segmentation_data.txt（nx ny xmin xmax ymin ymax + 地形行 + 标签行），
 * 提供标签场的双线性插值及梯度，用于可微的 rollable 区域约束。
 * 标签 0=lift, 1=rollable；插值后 [0,1]，约束时要求 >= 0.5 表示在 rollable 内。
 */
class SegmentationReader {
public:
    bool load(const std::string& path);

    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    double getXmin() const { return xmin_; }
    double getXmax() const { return xmax_; }
    double getYmin() const { return ymin_; }
    double getYmax() const { return ymax_; }

    /** 双线性插值标签 [0,1]，越接近 1 越在 rollable 内 */
    double interp(double x, double y) const;
    /** 插值对 x,y 的偏导 */
    void interpGrad(double x, double y, double& ddx, double& ddy) const;

private:
    int nx_ = 0, ny_ = 0;
    double xmin_ = 0, xmax_ = 0, ymin_ = 0, ymax_ = 0;
    std::vector<double> labels_;  // row-major labels_[j*nx_+i], 0 or 1
};

#endif
