#ifndef RBF_TERRAIN_HPP
#define RBF_TERRAIN_HPP

#include <vector>

/**
 * RBF (Radial Basis Function) 地形类
 * 用于计算地形高度和梯度
 */
class RBFTerrain {
public:
    /**
     * 构造函数
     * @param centers_x RBF中心x坐标
     * @param centers_y RBF中心y坐标
     * @param weights RBF权重
     * @param sigma RBF长度尺度参数
     */
    RBFTerrain(const std::vector<double>& centers_x,
               const std::vector<double>& centers_y,
               const std::vector<double>& weights,
               double sigma);
    
    /**
     * 计算地形高度
     * @param x x坐标
     * @param y y坐标
     * @return 地形高度
     */
    double height(double x, double y) const;
    
    /**
     * 批量计算地形高度
     * @param x x坐标数组
     * @param y y坐标数组
     * @param n 数组长度
     * @param h 输出的高度数组
     */
    void height(const double* x, const double* y, int n, double* h) const;
    
    /**
     * 计算地形梯度
     * @param x x坐标
     * @param y y坐标
     * @param hx 输出的∂H/∂x
     * @param hy 输出的∂H/∂y
     */
    void gradient(double x, double y, double& hx, double& hy) const;
    
    /**
     * 批量计算地形梯度
     * @param x x坐标数组
     * @param y y坐标数组
     * @param n 数组长度
     * @param hx 输出的∂H/∂x数组
     * @param hy 输出的∂H/∂y数组
     */
    void gradient(const double* x, const double* y, int n,
                  double* hx, double* hy) const;
    
    /**
     * 计算表面法向量（未归一化）
     * @param x x坐标
     * @param y y坐标
     * @param nx 输出的法向量x分量
     * @param ny 输出的法向量y分量
     * @param nz 输出的法向量z分量（通常为1）
     */
    void normal(double x, double y, double& nx, double& ny, double& nz) const;
    
    /**
     * 批量计算表面法向量（未归一化）
     */
    void normal(const double* x, const double* y, int n,
               double* nx, double* ny, double* nz) const;
    
    /**
     * 获取RBF中心数量
     */
    int numCenters() const { return static_cast<int>(centers_x_.size()); }
    
    /**
     * 获取RBF中心x坐标范围
     */
    void getXRange(double& xmin, double& xmax) const;
    
    /**
     * 获取RBF中心y坐标范围
     */
    void getYRange(double& ymin, double& ymax) const;
    
    /**
     * 获取RBF中心范围
     */
    void getBounds(double& xmin, double& xmax, double& ymin, double& ymax) const;

private:
    std::vector<double> centers_x_;
    std::vector<double> centers_y_;
    std::vector<double> weights_;
    double sigma_;
    double sigma2_;  // sigma^2，预计算以提高效率
};

#endif // RBF_TERRAIN_HPP

