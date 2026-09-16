# Experiment Protocol, Evaluation Methodology, Reproducibility Guide, and Ablation Plan

本文依据当前源码、CMake、README、脚本以及仓库内现存输入/输出文件编写，不把历史文件的存在解释为已经完成或验证过的实验结论。

全文使用三种状态：

- **Implemented and reproducible now**：当前源码已有对应程序和数据接口，可以按协议重新生成；是否数值收敛仍需记录实际运行状态。
- **Existing artifact but provenance uncertain**：文件确实存在，但缺少生成 commit、完整命令、参数、机器环境或 solver 日志，不能确认由当前源码生成。
- **Proposed experiment**：建议未来执行的实验或评估；本文没有运行，也没有结果。

# 1. Experimental Objectives

实验应分别回答以下问题，避免用单个成功截图代替定量验证：

1. RBF terrain 是否稳定、正确地提供规划所需的连续高度、梯度和法向？
2. 基于 \(\|\nabla f\|\le\mu\) 的 friction segmentation 是否产生合理且边界稳定的 rollable/lift 区域？
3. region-aware RRT* 相比无区域惩罚规划，是否降低 lift 区内的路径长度或占比，同时保持可接受的路径长度和规划时间？
4. rollable trajectory optimizer 是否生成平滑、满足端点和机构几何范围的轨迹？
5. lift-leg optimizer 是否能在 lift 区维持足够摆动间隙、限制穿透并保持分段轮高连续？
6. 输出轨迹是否满足所有实际注册的硬约束，而不是仅有较低目标值？
7. Ipopt 的状态、迭代时间、初值敏感性和参数敏感性如何？

这些问题目前没有一套仓库内自动汇总的实验报告。现有程序输出可覆盖其中一部分，其余属于 **Proposed experiment**。

# 2. Experiment Pipeline

```text
RBF JSON
  -> terrain_friction_segment
  -> terrain_segmentation_data.txt
  -> rrt_star_region
  -> rrt_star_path.txt + waypoints_segmented.txt
  -> first rollable/lift segment extraction
  -> seg_trajopt_cpp or lift_leg_trajopt_cpp
  -> metric evaluation
  -> non-interactive visualization
```

阶段边界如下：

| Stage | Input | Current executable/output | Status |
|---|---|---|---|
| Terrain reconstruction and segmentation | RBF JSON, \(\mu\) | `terrain_friction_segment`; height grid and labels | Implemented and reproducible now |
| Coarse planning | segmentation data, start/goal | `rrt_star_region`; path and labeled waypoints | Implemented, but isolated output routing is incomplete |
| Roll segment optimization | RBF, segmentation, labeled waypoints | `seg_trajopt_cpp`; initial/final trajectory | Implemented and reproducible now for the first label-1 segment |
| Lift segment optimization | RBF, labeled waypoints | `lift_leg_trajopt_cpp`; initial/final trajectory | Implemented and reproducible now for the first label-0 segment |
| Whole-route orchestration/stitching | all labeled segments | no implementation found | Proposed experiment/infrastructure |
| Metrics aggregation | stage outputs and logs | no unified evaluator found | Proposed experiment/infrastructure |

# 3. Current Test Assets

## 3.1 Inputs and coarse-planning artifacts

| Location | Content | Interpretation | Provenance status |
|---|---|---|---|
| `test/terrain_res/rbf.json` | RBF centers and weights | Concrete baseline input accepted by current readers. | Input asset exists; upstream sensor/run provenance uncertain |
| `test/terrain_res/segmentation.txt` | `512 x 512` binary label grid over `[-0.7,5.7] x [-3.7,2.7]` | Output format matches current segmentation writer. | Existing artifact but provenance uncertain |
| `test/terrain_res/terrain_segmentation_data.txt` | same grid header, RBF heights, labels | Direct input for current C++ RRT* and seg optimizer. | Existing artifact but provenance uncertain |
| `test/terrain_res/segmentation.png`, `segmentation_comparison.png` | segmentation visualizations | Qualitative figures only; no configuration embedded. | Existing artifact but provenance uncertain |
| `test/terrain_res/rrt_star_path.txt` | 16 planar path points from `(0.6,1.2)` to `(2.4,-3.0)` | Matches current path text contract. | Existing artifact but provenance uncertain |
| `test/terrain_res/waypoints_segmented.txt` | same 16 points with labels: 10 roll, 3 lift, 3 roll | Current readers can extract the first contiguous requested segment. | Existing artifact but provenance uncertain |
| `test/terrain_res/rrt_star_path.png` | coarse-path visualization | No node tree, seed, timing or command metadata embedded. | Existing artifact but provenance uncertain |

`terrain_friction_segment/segmentation.txt` 和 `segmentation.png` 也存在，且当前 SHA-256 分别与 `test/terrain_res` 中同名文件一致；这只能确认文件内容相同，不能确认生成 commit 或参数。根目录级 `terrain_friction_segment/segmentation_comparison.png` 同样是历史图像。

## 3.2 Fine-optimization artifacts

- `seg_trajopt_cpp/` 当前没有持久的 `seg_out` 或其他结果目录；只有源码、README、脚本和历史 build 目录。
- `lift_leg_trajopt_cpp/lift_leg_out` 目录存在但为空。
- 下列 12 个 lift 目录包含 `initial_trajectory.txt` 和 `final_trajectory.txt`：`lift_leg_out_rerun`、`lift_leg_out_rerun2` 至 `lift_leg_out_rerun9`、`lift_leg_out_smoke`、`lift_leg_out_stair_init`、`lift_leg_out_test`。
- 部分 lift 目录还包含 `seg_trajectory_plot.png`、`seg_trajectory_kinematics.png`、`initial_dense_samples_debug.csv/.png` 或 `chebyshev_fit_wheel_z.png`。

这些 lift 输出时间早于当前若干 lift 源码文件的修改时间，并且目录中没有 commit、完整命令、Ipopt 状态或参数快照。因此全部标记为 **Existing artifact but provenance uncertain**，不能作为当前版本通过测试的证据。

`trajopt_cpp/` 没有稳定的持久轨迹结果目录；当前程序主要写 `/tmp` 可视化桥接文件。仓库根目录的 `ipopt_log.txt` 和 `presentation.log` 缺少与某次配置的完整绑定，也应视为来源不确定的历史 artifact。

## 3.3 Test and visualization scripts

| File | What it actually does | Status |
|---|---|---|
| `test/casadi_test.py` | 求解标量 `(x-3)^2`，仅验证 CasADi 能调用 Ipopt。 | Implemented and reproducible now; smoke test, not project validation |
| `test/cupy_test.py` | 随机生成 `N=2000`,`M=20000` 数据，比较 NumPy、CuPy 和 fused CuPy RBF 的单次计时。无 seed、warm-up、重复统计或正确性误差输出。 | Implemented microbenchmark; not a reproducible performance claim |
| `terrain_friction_segment/scripts/*.py` | 分区/RRT 绘图以及独立 Python RRT* 备用实现。部分脚本调用 `plt.show()`。 | Implemented; headless behavior varies by script |
| `seg_trajopt_cpp/scripts/plot_seg_trajectory.py` | 从 initial/final trajectory 生成轨迹和运动学 PNG；支持 `--no-show`。 | Implemented and reproducible now |
| `lift_leg_trajopt_cpp/scripts/plot_seg_trajectory.py` | 同上，可额外读取 RBF；支持 `--rbf`、`--initial-only`、`--no-show`。 | Implemented and reproducible now |
| `trajopt_cpp/visualize_*.py` | 消费 `/tmp` 临时格式，显示地形、轨迹、接触力/平衡曲线。 | Implemented, but tied to transient files |

历史 build 目录和其中二进制来自何种工具链不可确认，不应纳入基线复现证据，也不应作为新实验的构建目录复用。

# 4. Environment

## 4.1 Confirmed dependencies

| Dependency | Required by | Version status |
|---|---|---|
| CMake | four C++ subprojects | minimum `3.10`; CUDA subproject minimum `3.18` |
| C++17 compiler | all C++ programs | exact compiler/version UNKNOWN |
| CUDA Toolkit/runtime | `terrain_friction_segment` segmentation target | exact version and GPU architecture UNKNOWN |
| Ipopt | three trajectory optimizers | exact version UNKNOWN |
| MUMPS | selected by all three C++ optimizer mains | availability depends on Ipopt installation; exact version UNKNOWN |
| Python 3 | prototypes, tests and visualization | exact version UNKNOWN |
| NumPy | Python prototypes/tests/plots | exact version UNKNOWN |
| Matplotlib | visualization | exact version UNKNOWN |
| CasADi | root Python NLP prototypes and smoke test | exact version UNKNOWN |
| CuPy | `test/cupy_test.py` only | exact version/CUDA compatibility UNKNOWN |

## 4.2 OS assumptions

当前代码和脚本明显面向 Linux/POSIX：使用 `python3`、Bash、`/tmp`、`/home/yizhe/...`、`unistd.h/getpid()`、`xdg-open`、链接库 `m` 以及 GCC 风格 `-Wall -Wextra -O3`。Windows/macOS 是否无需修改即可完整构建为 **UNKNOWN**。复现实验必须记录 OS、发行版、内核/驱动、编译器、CUDA、Ipopt 和 Python 包版本；仓库当前没有锁文件或容器定义。

# 5. Build Protocol

以下为 **Implemented and reproducible now** 的独立构建协议，使用新的 `_build_experiment/`，不复用仓库中的历史 `build*`。命令按 Linux shell 表示；本文没有执行这些命令。

```bash
cmake -S trajopt_cpp -B _build_experiment/trajopt_cpp -DCMAKE_BUILD_TYPE=Release
cmake --build _build_experiment/trajopt_cpp --parallel

cmake -S terrain_friction_segment -B _build_experiment/terrain_friction_segment -DCMAKE_BUILD_TYPE=Release
cmake --build _build_experiment/terrain_friction_segment --parallel

cmake -S seg_trajopt_cpp -B _build_experiment/seg_trajopt_cpp -DCMAKE_BUILD_TYPE=Release
cmake --build _build_experiment/seg_trajopt_cpp --parallel

cmake -S lift_leg_trajopt_cpp -B _build_experiment/lift_leg_trajopt_cpp -DCMAKE_BUILD_TYPE=Release
cmake --build _build_experiment/lift_leg_trajopt_cpp --parallel
```

构建日志应保存到新实验目录。多配置生成器的实际可执行路径可能含 `Release/`；这必须记录，不能假定。成功构建只验证编译与链接，不验证算法收敛。

# 6. Baseline End-to-End Experiment

## 6.1 Safe minimum baseline using existing coarse artifacts

以下协议不会覆盖 `test/terrain_res`；它从现有 baseline input/coarse artifacts 开始，把新轨迹写入唯一目录。`RUN_ID` 必须替换为不会与已有目录冲突的标识。

```bash
RUN=experiments/baseline_RUN_ID
mkdir -p "$RUN/trajectory/roll" "$RUN/trajectory/lift" "$RUN/figures" "$RUN/logs"

./_build_experiment/seg_trajopt_cpp/seg_trajopt_cpp \
  test/terrain_res "$RUN/trajectory/roll" \
  >"$RUN/logs/seg_trajopt.log" 2>&1

./_build_experiment/lift_leg_trajopt_cpp/lift_leg_trajopt_cpp \
  test/terrain_res "$RUN/trajectory/lift" \
  >"$RUN/logs/lift_trajopt.log" 2>&1

python3 seg_trajopt_cpp/scripts/plot_seg_trajectory.py \
  "$RUN/trajectory/roll" --no-show

python3 lift_leg_trajopt_cpp/scripts/plot_seg_trajectory.py \
  "$RUN/trajectory/lift" --rbf test/terrain_res/rbf.json --no-show
```

该协议可复现“第一段 roll + 第一段 lift”的求解调用，但使用的 segmentation/RRT artifacts 仍是 **Existing artifact but provenance uncertain**，所以不是从原始 RBF 全量重生成的严格端到端实验。

## 6.2 Terrain regeneration into a fresh directory

```bash
RUN=experiments/baseline_RUN_ID
mkdir -p "$RUN/input" "$RUN/segmentation" "$RUN/logs"
cp test/terrain_res/rbf.json "$RUN/input/rbf.json"

./_build_experiment/terrain_friction_segment/terrain_friction_segment \
  "$RUN/input/rbf.json" "$RUN/segmentation" 0.3 \
  >"$RUN/logs/terrain_segmentation.log" 2>&1

python3 terrain_friction_segment/scripts/visualize_terrain_segmentation.py \
  "$RUN/segmentation/terrain_segmentation_data.txt" --auto-close 1
```

注意：分区可执行程序本身会自动尝试启动可视化，且没有 `--no-viz` 开关；在无图形环境中可能阻塞或失败，必须纳入运行记录。

## 6.3 RRT* regeneration limitation

算法调用形式为：

```bash
./_build_experiment/terrain_friction_segment/rrt_star_region \
  "$RUN/segmentation/terrain_segmentation_data.txt" \
  0.6 1.2 2.4 -3.0 "$RUN/trajectory/rrt_star_path.txt"
```

但当前 CLI 只能重定向 path text。`waypoints_segmented.txt` 和 `rrt_star_path.png` 仍写到源码中硬编码的 `/home/yizhe/trajopt/test/terrain_res`，并会尝试打开 PNG。为了“不覆盖现有实验输出”，**不要在主工作副本中执行上述 RRT 命令**。在不修改代码的前提下，严格隔离的全量 RRT 复现目前不可用；需要在一次明确授权的后续任务中增加输出目录参数，或在隔离环境中准备与硬编码路径完全分离的仓库副本。

因此当前没有“一条命令从 RBF 生成完整轨迹”的实现；即使手工串接，seg/lift 也只优化第一段匹配标签，不遍历和拼接整条路径。

# 7. Evaluation Metrics

“Available”分为程序直接打印和可由当前持久输出计算两类；没有 evaluator 的项目不写成已经测量。

## 7.1 Terrain reconstruction / geometry

| Metric | Definition | Availability |
|---|---|---|
| Terrain height consistency | 在相同 query 点比较 CPU `RBFTerrain`、分区 CPU writer 和独立参考实现的高度，报告 max/RMS error。 | Recommended metric；当前没有一致性测试 |
| Gradient consistency | 解析梯度与中心差分梯度比较，报告 max/RMS relative error，并避开近零分母。 | Recommended metric；当前没有 checker |
| Normal unit error | \(\max|\|n\|_2-1|\)。 | Recommended metric；可从模型计算但未被现有程序汇总 |

## 7.2 Segmentation

| Metric | Definition | Availability |
|---|---|---|
| Rollable/lift area ratio | label `1`/`0` 数除以总网格数；均匀网格时等同面积占比。 | Available：分区程序直接打印 rollable count，文件可复算 |
| Boundary behavior | 标签边界长度、连通分量数、对微小 \(\mu\)/坐标扰动的 label flip rate。 | Recommended metric |
| Threshold consistency | 检查每格 `label == (sqrt(fx^2+fy^2) <= mu)`。 | Recommended metric |

## 7.3 RRT*

| Metric | Definition | Availability |
|---|---|---|
| Path length | \(\sum_i\|p_{i+1}-p_i\|_2\)。 | Available from `rrt_star_path.txt`; not directly printed |
| Planning time | `run()` wall time in ms。 | Available：直接打印 |
| Number of returned waypoints | path vector size。 | Available：直接打印 |
| Number of tree nodes | inserted tree nodes. | Recommended metric；未输出 |
| Roll/lift path ratio | 应按每条边的标签采样长度计算；仅按 waypoint 数是近似值。 | Recommended evaluator；point-label ratio 可由现有文件近似计算 |
| Region transition count | 相邻 waypoint label 改变次数。 | Available from labeled waypoint file |
| Region-weighted cost | 使用与 planner 相同的 20 点边积分和 `k_lift`。 | Recommended evaluator；最终 cost 未输出 |
| Success rate | 多 seed 达到 goal 的比例。 | Proposed；当前 seed 固定为 0 且 CLI 不暴露 seed |

## 7.4 Trajectory optimization

| Metric | Definition | Availability |
|---|---|---|
| Objective value/components | 总目标与分项目标。 | Available：seg/lift 直接打印；full TNLP 打印总目标，分项未统一持久化 |
| Ipopt status | ApplicationReturnStatus 和/或 SolverReturn。 | Available：直接打印；必须与轨迹文件一起保存 |
| Iteration count | Ipopt 迭代次数。 | Available only from captured Ipopt console/log;没有结构化字段 |
| Solve time | OptimizeTNLP wall time。 | Available：三个 C++ main 均打印；seg/lift 还打印部分阶段耗时 |
| Maximum constraint violation | \(\max_j\max(g_L-g_j,0,g_j-g_U)\)。 | Partially available：seg 自行报告若干分组最大值；full 依赖 Ipopt summary；lift 未汇总 |
| Mean constraint violation | 对正违反量求均值，并同时报告 violated count。 | Recommended metric |
| Trajectory path length | base/left/right 采样折线长度，分别报告。 | Available from trajectory text，尚无统一汇总脚本 |
| Smoothness | 二阶采样差分平方和；需明确按 `s` 还是物理时间。 | Available from text as an offline calculation；recommended to standardize |
| Minimum base terrain clearance | \(\min_i[z_b-f(x_b,y_b)]\)。 | Available from trajectory + RBF；未自动汇总 |
| Wheel clearance/penetration | \(z_w-f(x_w,y_w)\) 的 min/max/RMS。 | Available from trajectory + RBF；未自动汇总 |
| Wheel spacing | \(\|p_l-p_r\|\) 的 min/max/mean。 | Available from trajectory text；未自动汇总 |
| Base height stability | `z_b-f(x_b,y_b)-CLEAR_Z` 的 RMS/max。 | Available from trajectory + RBF；未自动汇总 |

## 7.5 Lift-leg and force metrics

| Metric | Definition | Availability |
|---|---|---|
| Maximum swing height | 最大 wheel-terrain clearance；应按实际 contact mask 分轮报告。 | Recommended metric；输出文件不保存 contact mask |
| Minimum swing clearance | swing mask 内最小 clearance。 | Recommended metric；需同步保存 mask |
| Maximum penetration | \(\max(f-z_w,0)\)。 | Available from lift trajectory + RBF；未自动汇总 |
| Piecewise continuity | 每个内部边界左右极限的 C0 jump。 | Recommended metric；采样输出未必恰好包含双侧极限，最好从系数计算 |
| Force balance residual | \(\|f_l+f_r+[0,0,-mg]\|\)。 | Available transiently in full force visualization；未持久保存为 metrics |
| Torque balance residual | \(\|(p_l-p_b)\times f_l+(p_r-p_b)\times f_r\|\)。 | Available transiently in visualization；未持久汇总 |
| Friction-cone violation | \(\max(\|f_t\|^2-\mu^2f_n^2,0)\) 和 negative-normal violation。 | Recommended metric；objective 可算但当前脚本未汇总数值 |

# 8. Ablation Experiments

本节全部是 **Proposed experiment**，没有运行结果。

| ID | Independent variable | Controlled variables | Primary metrics | Expected question |
|---|---|---|---|---|
| A. RBF sigma | reader/model \(\sigma\) | same centers, weights, query grid, start/goal | height/gradient change, label flip rate, path/constraint metrics | 地形平滑尺度如何传播到分区与轨迹？ |
| B. Friction coefficient \(\mu\) | segmentation \(\mu\) | same RBF/grid/planner seed | roll area ratio, transitions, lift path length | 摩擦假设如何改变模式分配？ |
| C. RRT* region penalty | `lift_cost_weight` including no penalty | same map, seed, iterations | path length, lift-length ratio, weighted cost, time | 绕行长度与 lift 暴露如何权衡？ |
| D. Chebyshev coefficients | K | same input, Ns, weights, tolerances | objective, violation, smoothness, time | 表达能力与病态性/计算量如何权衡？ |
| E. Collocation/sample points | Ns | same K and input | sampled/dense-grid violation, time | 稀疏采样是否漏掉中间违反？ |
| F. Smoothing objective | on/off or weight | full optimizer otherwise fixed | curvature proxy, length, violation | 平滑项对几何可行性和路径形状的贡献？ |
| G. Yaw alignment | on/off or weight | same initial guess and other terms | yaw-tangent error, constraint metrics | 航向项是否减少侧向姿态偏差？ |
| H. Initial guess proximity | seg/lift weight | same fitted xref | solve success, iterations, objective, deviation from coarse path | 邻近项提高收敛还是限制必要调整？ |
| I. Base height stability | weight | same clear height and other terms | clearance RMS/max, feasibility | 基座高度项对稳定性和几何约束的作用？ |
| J. Swing-wheel lift | weight/profile peak | same contact mask | swing clearance, penetration, path smoothness | 抬升目标是否足够跨越 lift 地形且不过度抬高？ |
| K. Force/contact objectives | contact/force/torque/friction terms individually removed | same full problem and initialization | balance/torque/cone residual, penetration, solve time | 各软力学项分别贡献什么？ |
| L. Segmentation vs none | region-aware pipeline vs `k_lift=1`/不使用标签 | same RBF/start/goal/budget | lift ratio, route length, transitions, fine-planner success | 模式分区是否优于无分区规划？ |

修改源码常量才能执行的消融必须在独立 branch/commit 中记录补丁；不能只写在实验笔记里。

# 9. Solver Sensitivity Experiments

以下均为 **Proposed experiment**：

| Factor | Protocol | Metrics |
|---|---|---|
| Initial guess | waypoint fit、扰动系数、零/线性参考；每种多次且保存 x0 | status, iterations, time, final objective/violation |
| Ipopt max iterations | 从小到大设置预算 | best violation/objective vs time, termination reason |
| Tolerance | 扫描 `tol` 和 acceptable tolerances，不能只改其中一个 | status, iterations, final primal/dual residual |
| Objective weights | 单项对数尺度扫描，其余固定 | per-term objective and physical metrics |
| Polynomial order | 与 Ablation D 相同，并用 dense validation grid | fit error, between-collocation violation, conditioning proxy |
| Collocation resolution | 与 Ablation E 相同 | solve cost and dense-grid feasibility |

因为当前 random RRT seed 固定，solver sensitivity 与 planner stochastic sensitivity 应分开设计。

# 10. Numerical Validation

## 10.1 Gradient finite-difference check

对每个 objective 在可行域内多个随机方向 \(v\) 比较

\[
g^Tv\quad\text{与}\quad\frac{J(x+hv)-J(x-hv)}{2h},
\]

并随 \(h\) 做尺度扫描。至少覆盖：

- RBF height gradient；
- Chebyshev 一、二阶导数；
- full 的 terrain distance、airborne、yaw、contact stiffness、friction cone 和 torque balance；
- seg/lift 的 base-height、yaw、initial proximity、swing-lift。

应特别针对当前已发现的实现事实建立回归检查：airborne `eval_f` 使用 threshold，而 `eval_grad_f` 没有减 threshold；contact/friction/normal-force 导数省略了部分地形法向位置依赖。检查结果在修复前应记录为 known discrepancy，不能调大容差掩盖。

## 10.2 Constraint Jacobian finite-difference check

逐列或用随机方向检查

\[
J_g(x)v\quad\text{与}\quad\frac{g(x+hv)-g(x-hv)}{2h}.
\]

至少覆盖 boundary、trajectory bounds、wheel-base distance、wheel distance、anti-crossing、wheel airborne、normal-force bound、rollable interpolation、penetration force、terrain penetration、piecewise-z continuity 和 coefficient norm。还需验证 `getNumConstraints()`、bounds、Jacobian sparse structure 与实际写入顺序完全一致。

当前仓库没有自动 finite-difference test harness；本节为 **Proposed experiment**。

# 11. Failure Cases

每次失败不得只记录“failed”。建议使用以下分类和必需证据：

| Failure class | Record at minimum |
|---|---|
| Infeasible | Ipopt status、最大/平均违反、最高违反的 constraint 名称与采样索引 |
| Ipopt divergence/stagnation | termination reason、迭代数、objective/violation 历史、初值 ID |
| NaN/Inf | 首个非有限变量/目标/约束、调用阶段、参数、栈或日志 |
| Terrain boundary | 越界主体、采样点、距 RBF center 包围盒距离 |
| Extreme slope | \(\|\nabla f\|\)、\(n_z\)、\(\mu\)、label 和位置 |
| RRT* unable to reach goal | seed、iterations、step/radii、tree node count、最近 goal 距离 |
| Segment discontinuity | 相邻段 base/wheel/yaw 及一阶导 jump |
| Wheel penetration | 最大 penetration、轮侧、采样位置、contact mode |
| Excessive lift height | 最大 clearance、参考 peak、对应轮与路径位置 |

失败运行也应保留配置和日志，但不要把未收敛的 `final_trajectory.txt` 标为成功结果。seg/lift 当前在 `finalize_solution` 后可能把非成功 status 的最终 iterate 当作可用输出，因此实验判定必须以保存的 status 和独立 violation 检查为准。

# 12. Result Directory Convention

以下是 **Proposed convention**；不要移动或重命名当前历史文件。

```text
experiments/
  <experiment_name>/
    manifest.yaml
    config/
      environment.txt
      build-options.txt
      planner.yaml
      optimizer.yaml
    input/
      rbf.json
      input-checksums.txt
    segmentation/
      terrain_segmentation_data.txt
      segmentation.txt
    planning/
      rrt_star_path.txt
      waypoints_segmented.txt
    trajectory/
      roll_<segment_id>/
      lift_<segment_id>/
      stitched/
    metrics/
      metrics.json
      constraint_residuals.csv
    figures/
    logs/
      build.log
      segmentation.log
      rrt.log
      ipopt_<segment_id>.log
```

目录必须是新建且唯一的；禁止用 `latest` 覆盖历史结果。输入应复制或只读引用并保存 checksum。

# 13. Experiment Record Template

```yaml
status: proposed | running | success | failed | historical-unknown
commit: <git commit plus dirty-worktree flag>
date_utc: <ISO-8601>
operator: <name-or-id>
machine:
  hostname: <value>
  os: <distribution/version>
  cpu: <model>
  ram: <capacity>
  gpu: <model>
  gpu_driver: <version>
  cuda: <version>
software:
  compiler: <name/version>
  cmake: <version>
  ipopt: <version>
  linear_solver: mumps
  python: <version>
  packages: <frozen package list>
input_rbf:
  path: <path>
  sha256: <hash>
start_goal: [sx, sy, gx, gy]
segmentation:
  sigma: <value>
  mu: <value>
  grid: [nx, ny]
planner:
  seed: <value>
  max_iter: <value>
  step_length: <value>
  goal_radius: <value>
  rewire_radius: <value>
  edge_samples: <value>
  lift_cost_weight: <value>
optimizer:
  kind: full | roll | lift
  n_coeff: <value>
  n_sample: <value>
  objective_weights: <mapping>
  constraint_parameters: <mapping>
  ipopt_options: <mapping>
result:
  ipopt_status: <numeric-and-symbolic>
  iterations: <value-or-UNKNOWN>
  metrics_file: <path>
notes: <free text>
```

# 14. Reproducibility Checklist

- [ ] 使用唯一的新实验目录，没有覆盖 `test/terrain_res`、`*_out*`、日志或图片。
- [ ] 记录 commit、dirty status 和所有本地补丁。
- [ ] 保存输入 RBF、分区数据和各输入文件 SHA-256。
- [ ] 记录 OS、CPU、RAM、GPU、驱动、CUDA、编译器、CMake、Ipopt/MUMPS 和 Python 包版本。
- [ ] 从空的新 build directory 配置并保存 CMake/build log。
- [ ] 保存 `sigma`、\(\mu\)、网格、start/goal 和全部 planner 参数。
- [ ] 保存 RRT seed；若代码仍固定 seed，明确记录该事实。
- [ ] 保存 K、Ns、变量布局、目标权重、约束参数和全部 Ipopt options。
- [ ] 捕获 stdout/stderr，包括 symbolic/numeric Ipopt status、迭代数和时间。
- [ ] 以独立 evaluator 在 dense grid 上重算约束违反，不以输出文件存在作为成功标准。
- [ ] 报告失败样本和总尝试数，不只保留成功运行。
- [ ] 区分 point-label ratio 与按路径长度计算的 region ratio。
- [ ] 所有消融一次只改变声明的 independent variable。
- [ ] 绘图使用无窗口模式；记录脚本和参数。
- [ ] 将 proposed、running、failed、success 和 historical-unknown 明确区分。

# 15. Recommended Paper Experiments

以下全部为 **Proposed experiment**，不是现有结果：

1. **Baseline**：固定公开 RBF 输入、start/goal 和参数，报告完整阶段时间、路径/区域指标、每段 solver status 和 dense-grid feasibility。
2. **Ablation**：执行第 8 节 A-L，给出均值、离散度、失败率和每项唯一变量控制。
3. **Terrain difficulty**：按最大坡度、坡度分布、局部曲率、lift 连通长度和高度变化分层测试。
4. **Friction changes**：统一扫描 segmentation 与 contact model 的 \(\mu\)，并额外研究两者不一致时的行为。
5. **Computational performance**：分开测 RBF evaluation、CUDA segmentation、RRT、NLP initialization/solve/postprocess；进行 warm-up 和重复试验。
6. **Robustness**：扰动 RBF weights/centers、start/goal、初值和 planner seed，报告成功率及置信区间。
7. **Trajectory feasibility**：用比优化 collocation 更密的采样验证几何、间隙、接触和摩擦残差。
8. **Real robot experiment**：在建立时间参数化、控制接口、安全状态机和同步传感器记录后，比较计划与实测轨迹/接触。当前仓库没有这些模块。

## Current Experimental Coverage

| Experiment | Code exists | Data exists | Result exists | Reproducible | Status |
|---|---|---|---|---|---|
| CasADi-to-Ipopt smoke test | Yes | Synthetic inline | Console only | Yes, given dependencies | Implemented and reproducible now |
| CuPy RBF timing microbenchmark | Yes | Random inline | No preserved result | Partially; no seed/repeats/environment capture | Implemented, methodology insufficient for a claim |
| RBF height reconstruction | Yes | `test/terrain_res/rbf.json` | Height grid embedded in segmentation data | Component yes; reference truth absent | Implemented and reproducible now; accuracy unvalidated |
| CUDA friction segmentation | Yes | RBF input exists | txt/PNG artifacts exist | Component yes on compatible CUDA system | Existing artifact provenance uncertain; implementation reproducible |
| Region-aware C++ RRT* | Yes | segmentation input exists | path/waypoints/PNG exist | Algorithm yes; safe isolated output routing no | Existing artifact provenance uncertain |
| First rollable segment optimization | Yes | required inputs exist | No current seg result directory | Yes to a new output directory | Implemented and reproducible now, not run in this review |
| First lift segment optimization | Yes | required inputs exist | 12 nonempty historical result directories | Yes to a new output directory | Existing artifacts provenance uncertain; implementation reproducible |
| General full trajectory optimization | Yes | RBF input exists | no stable persistent trajectory result | Runnable with dependencies; output capture incomplete | Implemented, no verified current result |
| Constraint finite-difference validation | No unified harness | Inputs exist | No report found | No | Proposed experiment |
| RRT penalty ablation | Parameters exist in code | Baseline data exists | No controlled result set found | Requires controlled runner/config changes | Proposed experiment |
| Objective/solver sensitivity suite | Objective/options exist | Baseline data exists | No controlled result set found | No automated suite | Proposed experiment |
| Whole-route roll/lift optimization and stitching | No | Coarse labeled route exists | No | No | Proposed infrastructure and experiment |
| Real-robot validation | No interface found | No synchronized robot dataset found | No | No | Proposed experiment |

