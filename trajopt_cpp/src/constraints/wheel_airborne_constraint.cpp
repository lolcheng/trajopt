#include "constraints/wheel_airborne_constraint.hpp"
#include "polynomial.hpp"
#include <cmath>

WheelAirborneConstraint::WheelAirborneConstraint(
    int n_coeff, int n_sample,
    std::shared_ptr<RBFTerrain> terrain,
    double min_ground_distance,
    double max_airborne_distance)
    : n_coeff_(n_coeff), n_sample_(n_sample),
      terrain_(terrain),
      min_ground_distance_(min_ground_distance),
      max_airborne_distance_(max_airborne_distance) {
    
    // 变量索引（根据TrajOptNLP中的布局）
    idx_clx_ = 4 * n_coeff_;
    idx_cly_ = 5 * n_coeff_;
    idx_clz_ = 6 * n_coeff_;
    idx_crx_ = 7 * n_coeff_;
    idx_cry_ = 8 * n_coeff_;
    idx_crz_ = 9 * n_coeff_;
    
    // 创建采样点 s ∈ [0, 1]
    s_samples_.resize(n_sample_);
    for (int i = 0; i < n_sample_; ++i) {
        s_samples_[i] = static_cast<double>(i) / (n_sample_ - 1);
    }
    
    // 预计算Vandermonde矩阵
    Phi_.resize(n_sample_ * n_coeff_);
    Polynomial::vandermonde(s_samples_.data(), n_sample_, n_coeff_, Phi_.data());
}

WheelAirborneConstraint::~WheelAirborneConstraint() = default;

int WheelAirborneConstraint::getNumConstraints() const {
    // 左轮和右轮各2个约束（下界和上界）
    // 下界约束：Zl - Hl >= min_ground_distance
    // 上界约束：Zl - Hl <= max_airborne_distance
    return 4 * n_sample_;
}

void WheelAirborneConstraint::evalTrajectory(const double* coeff, double* traj) const {
    Polynomial::eval(Phi_.data(), n_sample_, coeff, n_coeff_, traj);
}

void WheelAirborneConstraint::eval_g(const double* x, double* g) {
    // 提取系数
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轨迹
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(clz, Zl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    evalTrajectory(crz, Zr.data());
    
    // 计算地形高度
    std::vector<double> Hl(n_sample_), Hr(n_sample_);
    terrain_->height(Xl.data(), Yl.data(), n_sample_, Hl.data());
    terrain_->height(Xr.data(), Yr.data(), n_sample_, Hr.data());
    
    // 约束: min_ground_distance <= Zl - Hl <= max_airborne_distance
    //       min_ground_distance <= Zr - Hr <= max_airborne_distance
    // 转换为两个不等式约束：
    //   下界：g1 = min_ground_distance - (Zl - Hl) <= 0  (即 Zl - Hl >= min_ground_distance)
    //   上界：g2 = (Zl - Hl) - max_airborne_distance <= 0  (即 Zl - Hl <= max_airborne_distance)
    for (int i = 0; i < n_sample_; ++i) {
        double dist_l = Zl[i] - Hl[i];
        double dist_r = Zr[i] - Hr[i];
        
        // 左轮下界约束：Zl - Hl >= min_ground_distance
        g[i] = min_ground_distance_ - dist_l;
        // 左轮上界约束：Zl - Hl <= max_airborne_distance
        g[n_sample_ + i] = dist_l - max_airborne_distance_;
        
        // 右轮下界约束：Zr - Hr >= min_ground_distance
        g[2 * n_sample_ + i] = min_ground_distance_ - dist_r;
        // 右轮上界约束：Zr - Hr <= max_airborne_distance
        g[3 * n_sample_ + i] = dist_r - max_airborne_distance_;
    }
}

void WheelAirborneConstraint::eval_jac_g(const double* x, double* values, int offset) {
    // 提取系数
    const double* clx = x + idx_clx_;
    const double* cly = x + idx_cly_;
    const double* clz = x + idx_clz_;
    const double* crx = x + idx_crx_;
    const double* cry = x + idx_cry_;
    const double* crz = x + idx_crz_;
    
    // 评估轨迹
    std::vector<double> Xl(n_sample_), Yl(n_sample_), Zl(n_sample_);
    std::vector<double> Xr(n_sample_), Yr(n_sample_), Zr(n_sample_);
    
    evalTrajectory(clx, Xl.data());
    evalTrajectory(cly, Yl.data());
    evalTrajectory(clz, Zl.data());
    evalTrajectory(crx, Xr.data());
    evalTrajectory(cry, Yr.data());
    evalTrajectory(crz, Zr.data());
    
    // 计算地形高度和梯度
    std::vector<double> Hl(n_sample_), Hr(n_sample_);
    std::vector<double> Hlx(n_sample_), Hly(n_sample_);
    std::vector<double> Hrx(n_sample_), Hry(n_sample_);
    
    terrain_->height(Xl.data(), Yl.data(), n_sample_, Hl.data());
    terrain_->height(Xr.data(), Yr.data(), n_sample_, Hr.data());
    terrain_->gradient(Xl.data(), Yl.data(), n_sample_, Hlx.data(), Hly.data());
    terrain_->gradient(Xr.data(), Yr.data(), n_sample_, Hrx.data(), Hry.data());
    
    // 雅可比结构：
    // 下界约束 g1 = min_ground_distance - (Zl - Hl):
    //   ∂g1/∂clz = -Phi_row
    //   ∂g1/∂clx = Hlx[i] * Phi_row
    //   ∂g1/∂cly = Hly[i] * Phi_row
    // 上界约束 g2 = (Zl - Hl) - max_airborne_distance:
    //   ∂g2/∂clz = Phi_row
    //   ∂g2/∂clx = -Hlx[i] * Phi_row
    //   ∂g2/∂cly = -Hly[i] * Phi_row
    // 注意：offset是values数组的偏移，不是约束索引
    
    int idx = 0;
    (void)offset;  // offset在TrajOptNLP中已经处理，这里不需要使用
    
    // 左轮下界约束的雅可比（前n_sample个约束）
    for (int i = 0; i < n_sample_; ++i) {
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对clz的梯度: ∂g1/∂clz = -Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -Phi_row[j];
        }
        // 对clx的梯度: ∂g1/∂clx = Hlx[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = Hlx[i] * Phi_row[j];
        }
        // 对cly的梯度: ∂g1/∂cly = Hly[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = Hly[i] * Phi_row[j];
        }
    }
    
    // 左轮上界约束的雅可比（第2个n_sample个约束）
    for (int i = 0; i < n_sample_; ++i) {
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对clz的梯度: ∂g2/∂clz = Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = Phi_row[j];
        }
        // 对clx的梯度: ∂g2/∂clx = -Hlx[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -Hlx[i] * Phi_row[j];
        }
        // 对cly的梯度: ∂g2/∂cly = -Hly[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -Hly[i] * Phi_row[j];
        }
    }
    
    // 右轮下界约束的雅可比（第3个n_sample个约束）
    for (int i = 0; i < n_sample_; ++i) {
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对crz的梯度: ∂g1/∂crz = -Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -Phi_row[j];
        }
        // 对crx的梯度: ∂g1/∂crx = Hrx[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = Hrx[i] * Phi_row[j];
        }
        // 对cry的梯度: ∂g1/∂cry = Hry[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = Hry[i] * Phi_row[j];
        }
    }
    
    // 右轮上界约束的雅可比（第4个n_sample个约束）
    for (int i = 0; i < n_sample_; ++i) {
        const double* Phi_row = Phi_.data() + i * n_coeff_;
        
        // 对crz的梯度: ∂g2/∂crz = Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = Phi_row[j];
        }
        // 对crx的梯度: ∂g2/∂crx = -Hrx[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -Hrx[i] * Phi_row[j];
        }
        // 对cry的梯度: ∂g2/∂cry = -Hry[i] * Phi_row
        for (int j = 0; j < n_coeff_; ++j) {
            values[idx++] = -Hry[i] * Phi_row[j];
        }
    }
}

void WheelAirborneConstraint::getBounds(double* g_l, double* g_u, int offset) {
    // 约束1（下界）: g1 = min_ground_distance - (Zl - Hl) <= 0
    //   即: Zl - Hl >= min_ground_distance
    // 约束2（上界）: g2 = (Zl - Hl) - max_airborne_distance <= 0
    //   即: Zl - Hl <= max_airborne_distance
    // 注意：g_l和g_u指针已经在调用时偏移了offset，所以这里不需要再使用offset
    (void)offset;
    const double NEG_INF = -1e20;
    
    // 左轮下界约束（前n_sample个）
    for (int i = 0; i < n_sample_; ++i) {
        g_l[i] = NEG_INF;  // 下界：负无穷
        g_u[i] = 0.0;      // 上界：0
    }
    // 左轮上界约束（第2个n_sample个）
    for (int i = 0; i < n_sample_; ++i) {
        g_l[n_sample_ + i] = NEG_INF;  // 下界：负无穷
        g_u[n_sample_ + i] = 0.0;      // 上界：0
    }
    // 右轮下界约束（第3个n_sample个）
    for (int i = 0; i < n_sample_; ++i) {
        g_l[2 * n_sample_ + i] = NEG_INF;  // 下界：负无穷
        g_u[2 * n_sample_ + i] = 0.0;      // 上界：0
    }
    // 右轮上界约束（第4个n_sample个）
    for (int i = 0; i < n_sample_; ++i) {
        g_l[3 * n_sample_ + i] = NEG_INF;  // 下界：负无穷
        g_u[3 * n_sample_ + i] = 0.0;      // 上界：0
    }
}

int WheelAirborneConstraint::getJacobianStructure(int* iRow, int* jCol, int offset) {
    if (iRow == nullptr || jCol == nullptr) {
        // 只返回非零元素数量
        // 每个约束依赖3个变量组（clx, cly, clz 或 crx, cry, crz），每个变量组n_coeff个系数
        // 共4*n_sample个约束（左轮下界、左轮上界、右轮下界、右轮上界）
        return 4 * n_sample_ * 3 * n_coeff_;
    }
    
    // 每个约束依赖于对应的轮子轨迹系数
    // 左轮下界约束 g[i] 依赖于 clx, cly, clz (各n_coeff个系数)
    // 左轮上界约束 g[n_sample_ + i] 依赖于 clx, cly, clz (各n_coeff个系数)
    // 右轮下界约束 g[2*n_sample_ + i] 依赖于 crx, cry, crz (各n_coeff个系数)
    // 右轮上界约束 g[3*n_sample_ + i] 依赖于 crx, cry, crz (各n_coeff个系数)
    
    int nnz = 0;
    // 左轮下界约束
    for (int i = 0; i < n_sample_; ++i) {
        int row = offset + i;
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_clz_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_clx_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_cly_ + j;
            nnz++;
        }
    }
    // 左轮上界约束
    for (int i = 0; i < n_sample_; ++i) {
        int row = offset + n_sample_ + i;
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_clz_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_clx_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_cly_ + j;
            nnz++;
        }
    }
    // 右轮下界约束
    for (int i = 0; i < n_sample_; ++i) {
        int row = offset + 2 * n_sample_ + i;
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_crz_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_crx_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_cry_ + j;
            nnz++;
        }
    }
    // 右轮上界约束
    for (int i = 0; i < n_sample_; ++i) {
        int row = offset + 3 * n_sample_ + i;
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_crz_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_crx_ + j;
            nnz++;
        }
        for (int j = 0; j < n_coeff_; ++j) {
            iRow[nnz] = row;
            jCol[nnz] = idx_cry_ + j;
            nnz++;
        }
    }
    
    return nnz;
}

