#ifndef ROLLABLE_REGION_CONSTRAINT_HPP
#define ROLLABLE_REGION_CONSTRAINT_HPP

#include "constraint_base.hpp"
#include "segmentation_reader.hpp"
#include <memory>

/**
 * 约束路径落在 rollable 区域内。
 * 使用分割图的双线性插值，约束 threshold - interp(x,y) <= 0，即 interp >= threshold（在 rollable 内）。
 * threshold 默认 0.5；可设为 0.35 等略放宽边界。
 * 可选缓冲 buffer：约束 threshold - interp(x,y) <= buffer，
 * 等价于 interp >= threshold - buffer（允许少量越界，减少不可行）。
 * 默认仅约束基座（N_SAMPLE 条），也可选择同时约束左右轮（3*N_SAMPLE 条）。
 */
class RollableRegionConstraint : public ConstraintBase {
public:
    RollableRegionConstraint(int n_coeff, int n_sample,
                            std::shared_ptr<SegmentationReader> seg,
                            double threshold = 0.5,
                            bool include_wheels = false,
                            double buffer = 0.0);

    int getNumConstraints() const override { return include_wheels_ ? 3 * n_sample_ : n_sample_; }
    void eval_g(const double* x, double* g) override;
    void eval_jac_g(const double* x, double* values, int offset = 0) override;
    void getBounds(double* g_lb, double* g_ub, int offset = 0) override;
    int getJacobianStructure(int* iRow, int* jCol, int offset = 0) override;

private:
    int n_coeff_;
    int n_sample_;
    double threshold_;
    double buffer_;
    bool include_wheels_;
    std::shared_ptr<SegmentationReader> seg_;
    int idx_cbx_, idx_cby_, idx_clx_, idx_cly_, idx_crx_, idx_cry_;
    std::vector<double> s_samples_;
    std::vector<double> Phi_;
    void evalTrajectory(const double* coeff, double* traj) const;
};

#endif
