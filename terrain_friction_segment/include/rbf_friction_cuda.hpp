#ifndef RBF_FRICTION_CUDA_HPP
#define RBF_FRICTION_CUDA_HPP

#include <vector>
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
void launch_friction_cone(
    const double* d_cx, const double* d_cy, const double* d_w,
    int n_centers, double sigma2,
    double xmin, double xmax, double ymin, double ymax,
    int nx, int ny, double nz_min,
    int* d_labels);
#ifdef __cplusplus
}
#endif

/**
 * 在 GPU 上对地形网格采样，根据摩擦锥约束判断每点是否可滚动。
 * 满足摩擦锥：n_z >= 1/sqrt(1+mu^2) 的区域标为可滚动(1)，否则为需抬腿(0)。
 */
namespace rbf_friction_cuda {

/**
 * 参数
 */
struct Params {
    double xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0;
    int nx = 256, ny = 256;
    double mu = 0.6;
};

/**
 * 运行 GPU 分割
 * @param centers_x, centers_y, weights RBF 参数
 * @param sigma RBF sigma
 * @param params 网格与摩擦系数
 * @param out_labels 输出 nx*ny，1=可滚动，0=需抬腿（由本函数 resize 并写入）
 * @return 是否成功
 */
bool run(const std::vector<double>& centers_x,
         const std::vector<double>& centers_y,
         const std::vector<double>& weights,
         double sigma,
         const Params& params,
         std::vector<int>& out_labels);

/**
 * 将分割结果保存为文本
 * 格式：第一行 nx ny xmin xmax ymin ymax，随后 nx*ny 个 0/1（行优先）
 */
bool saveSegmentation(const std::string& path,
                      int nx, int ny,
                      double xmin, double xmax, double ymin, double ymax,
                      const std::vector<int>& labels);

/**
 * 将分割结果保存为 PNG 图像并可选打开
 * 绿=可滚动，红=需抬腿
 */
bool saveSegmentationImage(const std::string& path,
                          int nx, int ny,
                          const std::vector<int>& labels,
                          bool open_after_save = true);

/**
 * 生成并保存「原始地形图 | 分割图」对比图
 * 左：RBF 地形高度（蓝低→绿→红高）；右：分割（绿=可滚动，红=需抬腿）
 * @param terrain_height 与 labels 同网格，行优先，nx*ny
 */
bool saveTerrainSegmentationComparison(
    const std::string& path,
    int nx, int ny,
    const std::vector<double>& terrain_height,
    const std::vector<int>& labels,
    bool open_after_save = true);

} // namespace rbf_friction_cuda

#endif
