# 1. Deployment Scope

本仓库当前可部署的是一套 **offline / planning-side trajectory optimization** 工具链：输入离线 RBF 地形快照，经摩擦分区、二维粗规划和 segment-level 非线性优化，输出以归一化参数 `s in [0,1]` 表示的基座与双轮几何轨迹。

它不是完整机器人 runtime controller，也不是 RBF-LIO 本体。当前仓库没有 LiDAR/IMU 数据接入、在线 RBF 拟合、机器人状态估计、ROS 通信、时间参数化、低层控制或执行安全状态机。

未来系统边界应是：

```text
RBF-LIO / terrain estimator
  -> 带 frame、时间戳、单位和 sigma 的 RBF terrain snapshot
  -> 本仓库的 offline planner / optimizer
  -> 带时间、速度、加速度、接触模式的 trajectory
  -> robot-specific trajectory adapter
  -> low-level controller
  <- state feedback / execution status / stop or replan request
```

当前已实现的是中间的 terrain snapshot 到几何 segment trajectory。上游 RBF-LIO 接口和下游 robot controller 接口均为 **Not implemented in this repository**。因此，新 Ubuntu 机器可以复现离线候选轨迹生成，但仅靠当前仓库不能得到可安全直接下发实机的控制命令。

# 2. Supported Components

| Component | Language | Main dependency | GPU required? | Executable / entry | Input | Output |
|---|---|---|---|---|---|---|
| `terrain_friction_segment` | C++17 + CUDA | CUDA runtime | Yes | `terrain_friction_segment` | RBF JSON、输出目录、可选 `mu` | RBF 副本、label grid、高度+label 数据、PNG |
| `rrt_star_region` | C++17 | standard C++ | No | `rrt_star_region` | segmentation data、可选 start/goal/path output | 2D path、带标签 waypoints、PNG |
| `trajopt_cpp` | C++17 | Ipopt/MUMPS | No | `trajopt_cpp` | RBF JSON；起终点等在源码中 | stdout、`/tmp` 轨迹/地形/力可视化数据 |
| `seg_trajopt_cpp` | C++17 | Ipopt/MUMPS | No | `seg_trajopt_cpp` | 含 RBF、waypoints、segmentation 的 data dir；output dir | 第一段 rollable 的 initial/final trajectory、terrain grid |
| `lift_leg_trajopt_cpp` | C++17 | Ipopt/MUMPS | No | `lift_leg_trajopt_cpp` | 含 RBF、waypoints 的 data dir；output dir | 第一段 lift 的 initial/final trajectory |
| Python prototypes | Python | CasADi、NumPy | No；CasADi/Ipopt 需可用 | `solve_casadi.py`, `trajopt_sdf.py` | 当前主要为源码内路径与参数 | 原型求解与可视化；不串接 C++ pipeline |
| Visualization scripts | Python | NumPy、Matplotlib | No | 各模块 `scripts/*.py`、`visualize_*.py` | 文本轨迹、地形、分区或临时文件 | PNG 或交互窗口 |
| `test/cupy_test.py` | Python | CuPy | Yes | Python script | 随机合成数据 | CPU/GPU 单次计时；不参与部署链 |

只有 `terrain_friction_segment` 分区目标必须使用 CUDA。`rrt_star_region` 与三个 Ipopt 优化器是 CPU 程序。

# 3. System Dependencies

| Dependency | Classification | Used by | Notes confirmed from repository |
|---|---|---|---|
| C++17 compiler | Required | all C++ targets | CMake sets C++17；具体编译器/version 未锁定 |
| CMake | Required | all C++ targets | minimum 3.10；CUDA subproject minimum 3.18 |
| Ipopt development library | Required for trajectory optimizers | full/seg/lift | CMake uses `find_package`, `pkg-config`, then manual lookup |
| MUMPS | Required by current runtime configuration | full/seg/lift | every optimizer sets `linear_solver=mumps` |
| CUDA Toolkit | Required only for segmentation executable | terrain segmentation | CMake requires `CUDAToolkit`; default architectures `70;80;86` |
| Python 3 | Required for current visualization workflow; prototype-only otherwise | plots, prototypes, tests | exact version not specified |
| NumPy | Visualization/prototype | Python scripts | no locked version |
| Matplotlib | Visualization-only | plot scripts and C++-launched helpers | no locked version |
| CasADi | Prototype-only | root Python optimizers, smoke test | not required by C++ pipeline |
| CuPy | Optional test-only | `test/cupy_test.py` | must match installed CUDA; not required for C++ CUDA kernel |
| `pkg-config` | Build helper | Ipopt fallback discovery | required if `find_package(Ipopt)` fails |
| `stb_image_write.h` | Vendored | segmentation/RRT PNG output | already in repository |
| `matplotlibcpp.h` | Vendored | full optimizer support code | already in repository |

仓库没有依赖锁文件、Dockerfile、Conda environment 或统一安装脚本。版本兼容性必须在部署记录中显式保存。

# 4. Ubuntu Installation

以下命令是基于当前 CMake 与 `trajopt_cpp/INSTALL.md` 的推荐准备过程；**本文没有执行任何安装命令**。Ubuntu release 与包版本不应在未验证时固定。

## 4.1 Core build and optimization packages

```bash
sudo apt update
sudo apt install build-essential cmake pkg-config \
  coinor-libipopt-dev libmumps-dev \
  python3 python3-pip python3-venv
```

`coinor-libipopt-dev` 是三个 CMake 文件直接提示的 apt package。`libmumps-dev` 是当前 `linear_solver=mumps` 配置所需的求解器依赖；某些 Ubuntu 的 Ipopt 包可能已经传递依赖 MUMPS，但部署时仍应通过 Ipopt 启动日志确认 MUMPS 可用。

如果从源码构建 Ipopt，仓库旧安装说明还列出 `gfortran liblapack-dev libblas-dev libmumps-dev`。该路线的 Ipopt 版本、编译选项和安装前缀需由部署者选择，当前工程不锁定它们。

## 4.2 CUDA

安装与目标 NVIDIA 驱动、GPU 和 Ubuntu 版本兼容的 CUDA Toolkit，并确认 `nvcc` 与 CMake 可见。仓库没有指定 CUDA 版本，不应盲目照搬其他机器版本。至少检查：

```bash
nvidia-smi
nvcc --version
cmake --version
```

如果机器无 NVIDIA GPU/CUDA，仍可构建三个 Ipopt 优化器；但当前 CMake 把 CUDA 和 CPU-only RRT* 放在同一个 `terrain_friction_segment` project 中，因此不能直接只配置 RRT* 而绕过缺失的 CUDA toolchain。

## 4.3 Python environment

推荐使用项目外或仓库中未提交的虚拟环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy matplotlib
```

只有运行 Python NLP prototypes 时才安装 CasADi：

```bash
python3 -m pip install casadi
```

CuPy 的 package variant 必须与本机 CUDA 匹配；它仅用于性能测试，离线 C++ pipeline 不需要，本文不指定安装包名或版本。

# 5. Build

不要复用仓库内已有 `build*`；它们可能来自不同机器。以下从仓库根目录使用互相独立的新构建目录：

```bash
cmake -S trajopt_cpp -B _build_deploy/trajopt_cpp -DCMAKE_BUILD_TYPE=Release
cmake --build _build_deploy/trajopt_cpp --parallel

cmake -S terrain_friction_segment -B _build_deploy/terrain_friction_segment \
  -DCMAKE_BUILD_TYPE=Release
cmake --build _build_deploy/terrain_friction_segment --parallel

cmake -S seg_trajopt_cpp -B _build_deploy/seg_trajopt_cpp -DCMAKE_BUILD_TYPE=Release
cmake --build _build_deploy/seg_trajopt_cpp --parallel

cmake -S lift_leg_trajopt_cpp -B _build_deploy/lift_leg_trajopt_cpp \
  -DCMAKE_BUILD_TYPE=Release
cmake --build _build_deploy/lift_leg_trajopt_cpp --parallel
```

CUDA CMake 默认生成 `sm_70`, `sm_80`, `sm_86`。部署 GPU 不匹配时，在 configure 阶段显式指定其 compute capability，例如：

```bash
cmake -S terrain_friction_segment -B _build_deploy/terrain_friction_segment \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=<target_arch>
```

`<target_arch>` 必须依据实际 GPU/toolkit 确认，不能从本仓库推断。默认 CUDA 编译使用 `-O3 --use_fast_math`，这可能使 GPU 判定与严格 IEEE 参考计算产生微小边界差异。

构建成功后预期目标为：

```text
_build_deploy/trajopt_cpp/trajopt_cpp
_build_deploy/terrain_friction_segment/terrain_friction_segment
_build_deploy/terrain_friction_segment/rrt_star_region
_build_deploy/seg_trajopt_cpp/seg_trajopt_cpp
_build_deploy/lift_leg_trajopt_cpp/lift_leg_trajopt_cpp
```

多配置生成器可能在中间增加 `Release/`。

# 6. Input Data Preparation

## 6.1 Required RBF JSON contract

当前 C++ readers 只要求：

```json
{
  "cur_rbf_grid": [[x0, y0], [x1, y1]],
  "rbf_weight": [w0, w1]
}
```

- `cur_rbf_grid`：`N x 2` RBF centers；
- `rbf_weight`：长度 N 的 scalar weights；
- centers 与 weights 必须按索引一一对应；
- `sigma` 不从 JSON 读取，当前两个 C++ `JSONReader` 均硬编码 `0.14`。

当前代码按

```text
z = sum_i w_i * exp(-||[x,y]-center_i||^2 / (2*sigma^2))
```

解释数据。实机部署前必须把以下内容变成明确的数据契约：frame ID、center/weight 的长度单位、z 正方向、重力方向、时间戳、地图有效区域、sigma 的单位与来源、snapshot 版本。当前 JSON reader 不验证这些元数据。

## 6.2 `test/terrain_res`

`test/terrain_res/rbf.json` 可作为离线 smoke/baseline 输入；同目录还有已有 segmentation、RRT path 和 labeled waypoints，可直接驱动 segment optimizers。它们缺少生成 commit、命令和环境记录，不应作为实机地图或当前版本验证证据。

# 7. End-to-End Offline Deployment

下面按当前实际 CLI 给出顺序。为避免覆盖历史输出，应使用新的 output directory；但 RRT* 的硬编码输出限制使完全隔离的数据链尚不成立。

## Step 1: RBF terrain input

准备例如：

```text
/data/run_001/rbf.json
```

先离线检查键、数组长度、有限数、坐标范围和预期 frame。当前 executable 只检查解析与数组长度，不验证 frame/单位。

## Step 2: Terrain friction segmentation

实际 CLI：

```bash
terrain_friction_segment [rbf.json] [output_dir] [mu]
```

部署示例：

```bash
RUN=/data/trajopt/run_001
mkdir -p "$RUN/terrain" "$RUN/roll" "$RUN/lift" "$RUN/logs"

MPLBACKEND=Agg \
_build_deploy/terrain_friction_segment/terrain_friction_segment \
  /data/run_001/rbf.json "$RUN/terrain" 0.3 \
  >"$RUN/logs/segmentation.log" 2>&1
```

`mu` 默认 `0.3`，输入会 clamp 到 `[0.1,1.0]`。网格固定 `512 x 512`，范围为 center bounds 外扩 `0.3`。程序还会自动尝试调用 Python visualization；没有 `--no-viz` CLI。

## Step 3: Region-aware RRT*

实际 CLI：

```bash
rrt_star_region [data.txt] [start_x start_y goal_x goal_y] [path_out.txt]
```

示例：

```bash
_build_deploy/terrain_friction_segment/rrt_star_region \
  "$RUN/terrain/terrain_segmentation_data.txt" \
  0.6 1.2 2.4 -3.0 "$RUN/terrain/rrt_star_path.txt" \
  >"$RUN/logs/rrt.log" 2>&1
```

**当前硬编码限制：** 只有 `data.txt`、start/goal 和 `path_out.txt` 可覆盖。`waypoints_segmented.txt` 与 `rrt_star_path.png` 始终使用 `/home/yizhe/trajopt/test/terrain_res`，程序还会尝试 `xdg-open`/`open` PNG。这会破坏隔离部署，甚至覆盖仓库测试文件。在修复 output-dir CLI 前，不应在生产工作副本运行该步骤；可在一次性隔离副本中运行并显式收集硬编码路径下的产物。

RRT* 的 `max_iter=2000`、`step_length=0.15`、`goal_radius=0.2`、`rewire_radius=0.5`、`edge_samples=20`、`lift_cost_weight=5` 和 seed `0` 均硬编码在 `rrt_star_region.hpp/.cpp`，没有 CLI。

## Step 4: Waypoint segmentation

没有单独 executable。RRT* 输出端为每个 path point 查询最近 grid label，写：

```text
x y label
```

其中 `0=lift`, `1=rollable`。seg/lift 的 `WaypointReader` 再识别连续相同标签。当前文件没有显式 segment ID、时间戳或 frame ID。

## Step 5: Rollable/lift segment optimizer

实际 CLI：

```bash
seg_trajopt_cpp [data_dir] [output_dir]
lift_leg_trajopt_cpp [data_dir] [output_dir]
```

`data_dir` 必须含 `rbf.json` 和 `waypoints_segmented.txt`；roll optimizer 还尝试读取 `terrain_segmentation_data.txt`。因为 RRT 输出路由问题，部署者目前必须先把这些文件集中到同一隔离目录。

```bash
_build_deploy/seg_trajopt_cpp/seg_trajopt_cpp \
  "$RUN/terrain" "$RUN/roll" \
  >"$RUN/logs/roll_ipopt.log" 2>&1

_build_deploy/lift_leg_trajopt_cpp/lift_leg_trajopt_cpp \
  "$RUN/terrain" "$RUN/lift" \
  >"$RUN/logs/lift_ipopt.log" 2>&1
```

roll 只优化第一段 label `1`；lift 优先优化第一段 label `0`，找不到时回退第一段 label `1`。两者不会自动遍历整条路线，也不会拼接输出。必须保存并检查 Ipopt status；`final_trajectory.txt` 存在不代表求解成功。

独立通用优化器的实际 CLI 是：

```bash
_build_deploy/trajopt_cpp/trajopt_cpp /data/run_001/rbf.json
```

它使用源码内固定起终状态和参数，不消费 RRT waypoints，不是上述 segment pipeline 的汇总器。

## Step 6: Visualization/export

```bash
python3 seg_trajopt_cpp/scripts/plot_seg_trajectory.py "$RUN/roll" --no-show

python3 lift_leg_trajopt_cpp/scripts/plot_seg_trajectory.py \
  "$RUN/lift" --rbf "$RUN/terrain/rbf.json" --no-show
```

当前 export 是采样几何状态，不包含时间、速度、加速度、控制模式或执行命令。机器人可用格式需要第 11 节的 adapter。

# 8. Outputs

| Output | Current format and meaning |
|---|---|
| `rbf.json` | segmentation executable 对输入 JSON 的覆盖式副本，供后续模块使用。 |
| `segmentation.txt` | line 1: `nx ny xmin xmax ymin ymax`; line 2: `nx*ny` row-major labels。 |
| `segmentation.png` | green rollable / red lift visualization。 |
| `segmentation_comparison.png` | terrain/segmentation comparison image。 |
| `terrain_segmentation_data.txt` | line 1 grid metadata；line 2 `nx*ny` heights；line 3 `nx*ny` labels。 |
| `rrt_star_path.txt` | no header；each line `x y`，start to goal。 |
| `waypoints_segmented.txt` | no header；each line `x y label`，`0=lift`, `1=rollable`。 |
| `rrt_star_path.png` | region map and coarse path image。 |
| `initial_trajectory.txt` | header `s xb yb zb psi xl yl zl xr yr zr`；segment initial guess samples。 |
| `final_trajectory.txt` | same columns；Ipopt final iterate samples。不是 time trajectory。 |
| `terrain_grid.txt` | seg output；line 1 grid metadata，then one height per line on a `100 x 100` grid。 |
| stdout/stderr | objective terms、部分 constraint diagnostics、solve times 和 Ipopt status；必须重定向保存。 |

通用 `trajopt_cpp` 另外使用：

- `/tmp/rbf_terrain_data.txt`；
- `/tmp/trajectory_data.txt`；
- `/tmp/force_data_<pid>.txt`，force visualization 后删除。

这些是临时可视化协议，不是可靠部署 artifact。通用程序当前不会像 segment optimizers 一样把 initial/final trajectory 写入指定 output directory。

# 9. Headless Deployment

在 remote server、SSH、container-like 环境或无 `DISPLAY` 时：

1. 对 seg/lift plot 脚本使用已经实现的 `--no-show`。
2. 可设置 `MPLBACKEND=Agg`，使 Matplotlib 使用非交互 backend：

   ```bash
   export MPLBACKEND=Agg
   ```

3. `visualize_terrain_segmentation.py` 只有 `--auto-close [seconds]`，没有 `--no-show`；可手动运行：

   ```bash
   MPLBACKEND=Agg python3 terrain_friction_segment/scripts/visualize_terrain_segmentation.py \
     "$RUN/terrain/terrain_segmentation_data.txt" --auto-close 1
   ```

4. segmentation executable 会自动调用该 Python 脚本，当前无禁用选项。必须用 headless backend，并用日志/超时监控防止阻塞。
5. RRT* 总会尝试 `xdg-open` 或 `open` PNG，当前无禁用选项；无桌面程序时命令通常失败后继续，但不应依赖这一行为。
6. `trajopt_cpp` 的 C++ visualizers 会通过 `std::system` 启动 Python，默认 `auto_close=false`，也没有总 `--no-viz` CLI；它不适合作为无人值守部署入口。
7. 不要把 `plt.show()` 脚本放入 unattended pipeline，除非明确使用 Agg/auto-close 并验证进程会退出。

# 10. Runtime Configuration

当前没有统一 config file。下表均需在部署 manifest 中记录；标为 hard-coded 的参数变更需要修改源码并重新构建。

| Parameter | Current value/source | Configuration status |
|---|---|---|
| RBF `sigma` | `0.14`, both C++ `json_reader.cpp` | Hard-coded；JSON 不覆盖 |
| segmentation `mu` | default `0.3`, clamp `[0.1,1.0]` | CLI third argument |
| full friction `mu` | `0.6`, `trajopt_cpp/main.cpp` | Hard-coded |
| lift initial-contact slope `mu` | `0.55`, lift main | Hard-coded；只用于初值投影 |
| mass / gravity | `1.0 kg`, `9.81 m/s^2`, full main | Hard-coded |
| nominal wheel spacing | `0.5 m` | Hard-coded in all optimizer mains |
| full base-wheel nominal | `0.8 m` | Hard-coded objective parameter |
| base terrain clearance | `0.8 m` | Hard-coded all optimizer mains |
| lift stance/lift parameters | `0.005 m` stance, `0.24 m` base peak, `0.03 m` no-contact margin | Hard-coded lift main |
| Chebyshev K / samples | full `5/25`; roll `5/15`; lift `8/21` | Hard-coded |
| lift z segment estimate | positive net rise / `0.20 m` | Hard-coded algorithm |
| Ipopt tolerance | `1e-4` for full/seg/lift | Hard-coded |
| Ipopt max iterations | full `1000`; seg `200`; lift `220` | Hard-coded |
| Ipopt Hessian | limited-memory | Hard-coded option |
| Ipopt linear solver | MUMPS | Hard-coded option |
| objective weights | constructors and three optimizer `main.cpp` files | Hard-coded；详见 `ALGORITHM.md` |
| RRT parameters | `rrt_star_region::Params` defaults | Mostly hard-coded；CLI only start/goal/path input |

Python prototypes have separate globals and may diverge from C++。不要把 Python 参数修改视为 C++ runtime configuration。

# 11. Real-Robot Integration

## 11.1 Current repository capability

| Interface/capability | Current status |
|---|---|
| ROS 1/ROS 2 publisher | **Not implemented in this repository** |
| ROS 1/ROS 2 subscriber | **Not implemented in this repository** |
| ROS launch/parameter files | **Not implemented in this repository** |
| Standard/custom trajectory message | **Not implemented in this repository** |
| Robot controller interface | **Not implemented in this repository** |
| Live robot state feedback | **Not implemented in this repository** |
| Online RBF-LIO subscription | **Not implemented in this repository** |
| Real-time replanning | **Not implemented in this repository** |
| Command execution/acknowledgement | **Not implemented in this repository** |
| Emergency stop/fallback state machine | **Not implemented in this repository** |
| Time parameterization | **Not implemented in this repository** |

## 11.2 Recommended integration architecture

```text
RBF-LIO / State Estimator
  |  terrain snapshot + frame_id + stamp + covariance/validity
  v
Planner Input Adapter
  |  validated RBF + current robot state + goal + robot model
  v
Segmentation -> RRT* -> All-Segment Optimizer -> Stitcher
  |  continuous geometric path + contact schedule
  v
Time Parameterizer and Dynamic Feasibility Checker
  |  t, pose/state, velocity, acceleration, contact/force references
  v
Robot-Specific Trajectory Adapter
  |  controller message/service/action
  v
Low-Level Controller
  ^  state feedback, tracking error, contact status, fault/stop
```

建议 adapter schema 至少携带：version、frame ID、source map timestamp、trajectory timestamp、absolute/relative time、base pose、wheel/contact target、velocity、acceleration、contact mode、validity horizon、limits、planner status 和 checksum。

实机 planner 应以当前 robot state 初始化，而不是只使用源码或 waypoint 的理想起点；controller 接受轨迹前应进行 start-state tolerance handshake。以上均为推荐架构，不是当前实现。

# 12. Coordinate Frames

源码能确认的只有所有 RBF centers、base 和 wheels 在同一个数值笛卡尔坐标系中计算：

- terrain 是该平面上的图面 `z=f(x,y)`，向上法向 `[-df/dx,-df/dy,1]`；
- base 用点 `(xb,yb,zb)` 和 yaw `psi` 表示，没有 roll/pitch；
- 左右轮用接触/轮位置点 `(xl,yl,zl)`、`(xr,yr,zr)` 表示；
- yaw 前向向量是 `(cos psi, sin psi)`，横向向量是 `(-sin psi, cos psi)`。

源码/文件没有 frame ID。“world”可作为部署概念，但不能从现有 JSON 确认它究竟是 map、odom、local terrain 还是某个 RBF-LIO frame。也没有 base-to-wheel TF/extrinsic、轮轴方向、轮半径、关节 frame 或 RBF snapshot 到 robot localization frame 的变换定义。

实机部署前必须明确：

- RBF frame 与机器人 localization frame；
- handedness、x/y heading 和 z-up/gravity convention；
- yaw 零方向和正方向；
- base reference point 与真实 base link 的关系；
- wheel point 表示轮心还是接触点；
- snapshot 时间同步和 transform 有效期。

# 13. Units

以下单位由源码常量、注释与物理公式使用方式支持，但 JSON 本身没有单位字段：

| Quantity | Unit |
|---|---|
| position, terrain height, clearance, distance | meter (m) |
| yaw angle | radian (rad) |
| contact force | newton (N) |
| mass | kilogram (kg) |
| gravity / acceleration parameter | meter per second squared (m/s^2) |
| contact stiffness | newton per meter (N/m) |
| friction coefficient | dimensionless |
| RBF sigma | meter, assuming center coordinates are meters |
| normalized path parameter `s` | dimensionless |

当前轨迹没有 seconds 字段，不能从相邻 `s` sample 推断速度单位。

# 14. Safety Checks Before Robot Execution

下表是 deployment requirements，不代表当前代码全部实现。

| Check | Current coverage | Required before robot execution |
|---|---|---|
| Trajectory x/y bounds | sampled hard constraints exist | dense continuous validation plus map-validity margin |
| Velocity limits | Not implemented | time parameterization and per-DOF limits |
| Acceleration/jerk limits | Not implemented | dynamic profile and controller limits |
| Wheel spacing/base-wheel geometry | sampled constraints exist | dense validation against real kinematic model |
| Base terrain clearance | endpoint constraints + soft objective | hard minimum clearance over full timed trajectory |
| Terrain penetration | partial full/lift constraints | dense check with uncertainty margin |
| Contact state | lift has internal masks but does not export them | explicit synchronized contact schedule and sensing feedback |
| Force limits | partial normal bound/penetration proxy | actuator/contact limits and validated dynamics |
| Friction feasibility | segmentation + full soft penalty | hard/robust check using estimated friction uncertainty |
| Segment continuity | only lift wheel-z C0 internal links | base/wheel/contact continuity across every route segment, preferably C1/C2 as required |
| Start-state consistency | ideal endpoint constraint only | compare live state to first command; reject or replan on mismatch |
| Collision/self-collision | only simplified geometry | full robot/environment collision validation |
| Solver success | status printed | require accepted status plus independent max-violation threshold |
| Stop/fallback behavior | Not implemented | watchdog, timeout, controlled stop, safe stance and planner/controller fault handling |
| Map age/frame consistency | Not implemented | reject stale or transform-inconsistent RBF snapshots |

没有完成这些检查时，`final_trajectory.txt` 只能视为研究用候选几何轨迹。

# 15. Known Deployment Limitations

## Critical

1. 没有时间参数化、机器人通信层、controller adapter、状态反馈或安全执行状态机；当前输出不能直接作为实机命令。
2. 没有统一 runner、全 segment 遍历和自动 stitching；当前 seg/lift 各处理第一段匹配区域。
3. frame、单位 schema、RBF snapshot 时间戳和 robot/RBF transform 没有代码级契约。
4. RRT* 的 waypoint/PNG 输出目录硬编码为 `/home/yizhe/trajopt/test/terrain_res`，无法通过 CLI 隔离部署输出。

## Important

1. 多个入口、Python launcher 和 `/tmp` 输出包含 `/home/yizhe/...` 或固定路径。
2. `sigma`、机器人几何、质量、权重、K/Ns、Ipopt options 和大部分 planner 参数硬编码，配置无法审计式注入。
3. `trajopt_cpp` 不消费 RRT waypoints，且不提供稳定 output-dir export；它与 segment pipeline 是平行实现。
4. trajectory constraints 只在有限 samples 检查，没有部署级 dense continuous certification。
5. seg/lift 的 final iterate 可能在非成功 Ipopt status 下仍被写出/视为可用；必须外部判定。
6. force/friction/balance 在通用求解器中多为 soft penalty；lift 只使用穿透力 proxy。
7. Python prototypes 与 C++ 实现在变量数、采样数、目标/约束和运行入口上并非统一 runtime。

## Convenience

1. 没有 root-level superbuild、依赖锁、container image 或安装 target。
2. 某些 executable 自动启动 Matplotlib 或 `xdg-open`，缺少统一 `--headless/--no-viz`。
3. 输出文本没有 schema version、metadata sidecar、checksum 或结构化 metrics。
4. 历史 build/results 缺少 provenance，不适合直接复制到新机器。

# 16. Deployment Checklist

## Fresh machine to offline trajectory generated successfully

- [ ] 确认 Ubuntu、CPU/GPU、NVIDIA driver 和 CUDA compatibility。
- [ ] 安装 C++17 toolchain、CMake、Ipopt/MUMPS 和需要的 Python 环境。
- [ ] 记录所有工具与库版本；不要复用旧 build cache。
- [ ] 从四个独立 build directories 进行 Release configure/build。
- [ ] 在 headless 环境设置 `MPLBACKEND=Agg`，确认没有 GUI 阻塞。
- [ ] 校验 RBF JSON keys、数组长度、有限值、frame、单位、timestamp 和 sigma 约定。
- [ ] 使用新的 run directory，保存输入 checksum、命令和 stdout/stderr。
- [ ] 运行 CUDA segmentation，检查退出码、label 比例、范围和输出尺寸。
- [ ] 在隔离环境处理 RRT* 硬编码输出，确认 start/goal、path 与 waypoint labels。
- [ ] 将 `rbf.json`、`terrain_segmentation_data.txt`、`waypoints_segmented.txt` 放入同一 data directory。
- [ ] 运行 roll/lift optimizer 到新的 output directories。
- [ ] 保存并检查 Ipopt symbolic/numeric status、objective、time 和 constraint diagnostics。
- [ ] 用 `--no-show` 生成图；独立检查间隙、轮距、穿透和段端状态。
- [ ] 明确输出仅覆盖第一段 roll/lift，而非完整路线。

## Offline output to ready for robot integration

以下 checklist 当前不能仅靠本仓库完成：

- [ ] 定义并实现 versioned terrain snapshot interface。
- [ ] 定义 frame/TF、时间同步、单位和 map validity contract。
- [ ] 实现全部 segment orchestration 和连续 trajectory stitching。
- [ ] 将 `s` path 转成满足速度、加速度和 jerk 限制的 timed trajectory。
- [ ] 使用真实机器人运动学/动力学和 collision model 做 dense validation。
- [ ] 导出 contact schedule、force/reference semantics 和 uncertainty margins。
- [ ] 实现 robot-specific trajectory message/action/service adapter。
- [ ] 接入 live state feedback、start-state handshake、tracking-error monitor 和 replan policy。
- [ ] 实现 watchdog、stop/fallback、timeout 和 controller rejection handling。
- [ ] 完成 simulation/HIL 验证、受控场地测试和正式安全评审。

只有第二组也通过后，才能把轨迹标为“ready for robot integration”；这不等于已经达到生产级实机安全认证。

