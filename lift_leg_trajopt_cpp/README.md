# lift_leg_trajopt_cpp：台阶地形抬腿轨迹优化

本项目用于 `lift leg` 路段轨迹优化，目标是让轮式底盘在台阶/分层地形上以“交替抬腿”的方式通过。  
实现上复用了 `trajopt_cpp` 的基函数 + Ipopt NLP 框架，并在初始估计和约束建模上针对抬腿场景做了强化。

---

## 1. 轨迹参数化与优化变量

- 路径参数 `s ∈ [0,1]`，在 `N_SAMPLE` 个均匀点离散。
- 所有状态用切比雪夫基表示：  
  `q(s) = Σ c_k * T_k(2s-1)`。
- 当前决策变量为 **10 组系数**（共 `10 * N_COEFF`）：
  - 基座：`cbx, cby, cbz, cpsi`
  - 左轮：`clx, cly, clz`
  - 右轮：`crx, cry, crz`

默认参数（`main.cpp`）：

- `N_COEFF = 6`
- `N_SAMPLE = 21`

---

## 2. 台阶感知初始估计（核心）

初值生成不再是单次抬腿，而是自动根据地形分层做多次“上一级台阶”动作：

1. 读取 lift-leg 路段 waypoint（优先 `segment_type == 0`）。
2. 计算每个 waypoint 的基座地形高程 `h_b`。
3. 用 `LEVEL_MERGE_TOL` 对 `h_b` 聚类，得到离散层级 `level_id`。
4. 检测层级上升位置（`level_id[i] > level_id[i-1]`）作为“上台阶事件”。
5. 每个事件生成一个抬腿窗口，并左右腿交替摆动：
   - 摆动腿在窗口内采用抬腿曲线 + 层级过渡高度
   - 接触腿保持小离地间隙
6. 摆动腿初值强制离地：
   - `z_swing >= h(x,y) + SWING_NO_CONTACT_MARGIN`

这套逻辑对应你提的“每次抬腿都上一级台阶”。

---

## 3. 接触腿摩擦锥代理（初值阶段）

在初始估计阶段加入了一个轻量“摩擦锥可行性代理”：

- 使用局部坡度近似地面切向需求，要求 `|∇h| <= μ`（`μ = FRICTION_MU`）。
- 对接触腿足点 `(x,y)` 做局部搜索投影 `projectToFrictionSafeXY()`：
  - 若当前点坡度过大，自动挪到邻域内坡度更小的位置。

> 说明：这是初值可行性代理，不是完整刚体动力学下的严格摩擦锥 `||F_t|| <= μ F_n`。  
> 严格版需要显式切向力变量与动力学方程耦合。

---

## 4. 约束建模（无显式接触力变量）

### 4.1 边界与几何

- `BoundaryConstraint10Var`：
  - 基座起终点 `(x,y,z,psi)` 约束
  - 轮子起终点几何约束（相对基座）
  - 轮子起终点贴地约束（`z - h(x,y)`）
- `TrajectoryBoundsConstraint8Var`：`x/y` 落在地形包围盒附近
- `WheelBaseDistanceTerrainConstraint`：基座到左右轮距离区间约束
- `WheelsDistanceTerrainConstraint`：左右轮距区间约束
- `CoefficientNormConstraint`：系数范数上界

### 4.2 穿透力近似（替代显式接触力决策变量）

你提出“接触力加入变量会拖慢求解”，当前实现采用：

1. `WheelTerrainPenetrationConstraint`  
   约束轮地穿透上限：
   - `h(x_w,y_w) - z_w <= MAX_PENETRATION`

2. `PenetrationBasedForceConstraint`  
   用穿透估计法向力：
   - `F_n = k * smooth_positive(h - z)`
   - `k = CONTACT_STIFFNESS`

并基于接触时序掩码（`left_contact_mask/right_contact_mask`）施加分段约束：

- 接触相：`STANCE_FORCE_MIN <= F_n <= STANCE_FORCE_MAX`
- 摆动相：`F_n <= SWING_FORCE_MAX`

这样保留“力可解释性”，同时避免把 6 维接触力显式作为优化变量。

---

## 5. 目标函数

- `InitialGuessProximityObjective`：解不偏离结构化初值
- `BaseHeightStabilityObjective`：基座离地高度稳定在 `CLEAR_Z`
- `YawPathAlignmentObjective`：偏航与路径切向对齐
- `SwingWheelLiftObjective`：强化摆动腿抬升形状与接触腿稳定贴地

---

## 6. 运行

```bash
cd lift_leg_trajopt_cpp
cmake -S . -B build_lift_new
cmake --build build_lift_new -j4
./build_lift_new/lift_leg_trajopt_cpp ../test/terrain_res lift_leg_out_stair_init
```

---

## 7. 可视化

```bash
python3 scripts/plot_seg_trajectory.py lift_leg_out_stair_init --rbf ../test/terrain_res/rbf.json --margin 0.8
```

无界面模式：

```bash
python3 scripts/plot_seg_trajectory.py lift_leg_out_stair_init --rbf ../test/terrain_res/rbf.json --margin 0.8 --no-show
```

输出：

- `seg_trajectory_plot.png`：2D/3D 轨迹 + RBF 地形
- `seg_trajectory_kinematics.png`：运动学曲线

脚本在缺少 `terrain_grid.txt` 时会自动从 `rbf.json` 重建“轨迹附近”高程图。

---

## 8. 当前局限与后续方向

- 当前摩擦约束是“坡度代理”，非严格力锥。
- 若需要严格物理一致性，建议下一步扩展为：
  - 显式接触力变量 `(F_t, F_n)`
  - 严格摩擦锥 `||F_t|| <= μF_n`
  - 与机体动力学/力矩平衡联立优化
