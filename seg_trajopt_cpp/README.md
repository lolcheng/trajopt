# seg_trajopt_cpp：第一段 Rollable 路径轨迹优化

## 1. 路径参数与基函数

- 归一化路径参数 `s ∈ [0,1]`，在 `N_SAMPLE` 个均匀采样点上施加约束与目标。  
- 每一标量轨迹表示为  
`q(s) = sum_{k=0}^{N_COEFF-1} c_k * T_k(2s-1)`  
其中 `T_k` 为第一类切比雪夫多项式，自变量 `x = 2s-1 ∈ [-1,1]`。  
- 实现上通过 `BasisFunction::computeBasisMatrix` 与 `SegTrajOptNLP` 中存储的布局一致。

当前默认：`N_COEFF = 5`，`N_SAMPLE = 15`。

---

## 2. 决策变量（8 × N_COEFF 维）

向量 `x` 按块排列：

1. `cbx`, `cby`, `cbz`：基座位置 `(x_b, y_b, z_b)` 的系数
2. `cpsi`：偏航 `psi` 的系数
3. `clx`, `cly`：左轮接触点水平位置 `(x_l, y_l)` 的系数
4. `crx`, `cry`：右轮接触点水平位置 `(x_r, y_r)` 的系数

**轮心高度**不显式优化：`z_l = h(x_l, y_l)`，`z_r = h(x_r, y_r)`。

---

## 3. 初始猜测

1. 从 `waypoints_segmented.txt` 读取第一段 rollable 的 `(x,y)`。
2. `z_b = h(x,y) + CLEAR_Z`（默认 `CLEAR_Z = 0.8`）。
3. `psi` 由相邻 waypoint 差分方向得到：`atan2(dy, dx)`。
4. 左右轮：在垂直于前进方向的两侧偏移 `WHEELS_NOMINAL/2`（默认轮距名义 0.5 m）。
5. 对每个分量在 waypoints 上做**最小二乘**拟合切比雪夫系数（与优化用同一基矩阵）。

---

## 4. 目标函数

1. **InitialGuessProximity**
  `J_1 = w_x0 * ||x - x0||^2`  
   `w_x0 = 1e-2`（`x0` 为上一节系数向量）。
2. **BaseHeightStability**
  在采样点 `s_i`：  
   `e_i = z_b(s_i) - h(x_b(s_i), y_b(s_i)) - CLEAR_Z`  
   `J_2 = w_h * sum_i e_i^2`，`w_h = 2 / N_SAMPLE`。
3. **YawPathAlignment**（来自 `trajopt_cpp`）
  用 `(x_b,y_b)` 的离散切向与 `psi` 对齐，惩罚 `(1 - cos_alignment)^2` 的加权和，权重 `w_yaw = 10 / N_SAMPLE`。  
   注：该目标内部用 `Polynomial::vandermonde` 与主问题切比雪夫基混用时，仅作**近似一致**；若需严格同一基，需在源码中改为 Chebyshev 基矩阵。

**总目标**：`J = J_1 + J_2 + J_3`（各模块梯度相加，由 Ipopt 求和）。

---

## 5. 硬约束

### 1 BoundaryConstraint8Var（16 维，区间形式）

在 `s=0` 与 `s=1` 上：

- 基座：`x_b,y_b,z_b,psi` 与由端点 waypoint 与 `CLEAR_Z`、地形决定的期望值之差，限制在容差内（基座位置/偏航 `1e-3`，轮端 `(x,y)` 几何 `5e-2`）。  
- 轮心 `(x_l,y_l)`、`(x_r,y_r)` 与由基座、`psi` 及 `WHEELS_NOMINAL` 推出的几何目标之差，同上轮端容差。

`z` 端点目标为 `h(x_endpoint, y_endpoint) + CLEAR_Z`。

### 2 TrajectoryBoundsConstraint8Var

对所有采样点：`x_b,y_b,x_l,y_l,x_r,y_r` 落在 RBF 中心包围盒加 `TERRAIN_MARGIN` 的矩形内。

### 3 CoefficientNormConstraint（1 维）

`g = (1 / n_vars) * sum_k x_k^2 <= COEFF_NORM_MAX^2`  

默认 `COEFF_NORM_MAX = 10000`，`n_vars = 8 * N_COEFF`。

### 4 RollableRegionConstraint（若成功读取分割文件）

在采样点上对分割标签 `ell(x,y)`（双线性插值）施加：  
`ell >= threshold - buffer`（默认 `threshold = buffer = 0.1`）。  
默认 **基座 + 左右轮** 均约束（`ROLLABLE_INCLUDE_WHEELS = true`）。

### 5 WheelBaseDistanceTerrainConstraint

对每个采样点、左右两侧：

`d^2 = ||p_b - p_w||^2`（含 `z`，轮用 `h(x_w,y_w)`）  

`WHEEL_BASE_MIN^2 <= d^2 <= WHEEL_BASE_MAX^2`  

默认 `0.2 <= d <= 1.2`（米）。

### 6 WheelsDistanceTerrainConstraint

`d_{lr}^2 = ||p_r - p_l||^2`（含 `z_l, z_r` 由地形）  

`WHEELS_MIN^2 <= d_{lr}^2 <= WHEELS_MAX^2`  

默认 `sqrt(0.1) <= d_{lr} <= 1.0`（即平方在 `[0.1, 1]`）。

## 6. 运行
./run_and_visualize.sh ../../test/terrain_res seg_out
