# 地形摩擦锥分割 (GPU)

根据 RBF 地形在网格上采样，用 **GPU (CUDA)** 判断每个采样点是否满足摩擦锥约束，将地形分为：
- **可滚动区域**：满足摩擦锥，轮足可直接滚动通过
- **需抬腿区域**：不满足摩擦锥，轮足需抬腿通过

满足条件为：表面法向 `n` 的 `n_z >= 1/sqrt(1+μ²)`（竖直支撑力在锥内）。默认 **μ=0.3**（较严，便于把台阶、陡坡判为需抬腿）；RBF 会平滑地形，用更小的 μ 可补偿。

## 原理总览

项目由两段串联：
- **地形分割**（`terrain_friction_segment`）：RBF + GPU 摩擦锥判定，得到 rollable/lift 二值地图。
- **路径规划**（`rrt_star_region`）：在同一分割图上做带区域代价的 RRT*，得到粗路径和带标签 waypoint。

数据链路如下：
1. 读取 RBF JSON（中心点 + 权重 + sigma）
2. 在规则网格上并行计算法向并做摩擦锥判定（CUDA）
3. 输出分割与可视化数据到 `test/terrain_res/`
4. RRT* 读取 `terrain_segmentation_data.txt`，按区域代价搜索路径
5. 输出 `rrt_star_path.txt` 与 `waypoints_segmented.txt`（细规划输入）

## 分割模块原理（CUDA）

### 1) RBF 地形与法向

地形高度建模为高斯基函数叠加：

- `H(x, y) = Σ w_i * exp(-0.5 * ||p-c_i||^2 / sigma^2)`
- 梯度 `Hx, Hy` 由上式对 `x,y` 求导得到
- 未归一化法向取 `[-Hx, -Hy, 1]`，再归一化得到 `n = [nx, ny, nz]`

对应实现：
- `cuda/rbf_friction_kernel.cu`：`rbf_normal_at(...)`
- `src/main.cpp`：`computeTerrainHeightGrid(...)`（CPU 侧仅用于可视化高度图）

### 2) 摩擦锥判定

对每个网格点独立判断：
- 阈值 `nz_min = 1/sqrt(1+mu^2)`
- 若 `nz >= nz_min` -> `label=1`（rollable）
- 否则 `label=0`（lift）

这一步完全并行，适合 GPU；每个线程处理一个 `(i,j)` 网格点。

### 3) 为什么 mu 越小越“严格”

该判定形式下，`mu` 变小会使 `nz_min` 变大，要求表面更“接近水平”，因此更容易把台阶/陡坡判为 lift。项目默认 `mu=0.3` 用于增强这类结构的识别。

## 路径规划原理（区域代价 RRT*）

### 1) 状态空间与采样

- 在分割图覆盖的 `(x,y)` 连续空间采样
- 10% 概率直接采样目标点（goal bias），其余均匀随机采样
- `steer` 以固定步长 `step_length` 从最近节点向采样点扩展

### 2) 边代价（核心）

一条边不是简单欧氏长度，而是按区域积分近似：
- 把边离散为 `edge_samples` 个采样点
- 每个采样点查询 label
- 单段代价累计为：
  - rollable 区域权重 `1.0`
  - lift 区域权重 `k`（`lift_cost_weight`，默认 5.0）

即：
- `cost = L_rollable + k * L_lift`

这样路径会自动偏向可滚动区域，同时允许必要的短距离 lift 穿越。

### 3) RRT* 选择父节点 + 重连（rewire）

- 新节点先基于最近邻得到初始父节点
- 在 `rewire_radius` 内搜索候选父节点，选择总代价最小者
- 新节点加入后，再尝试重连半径内旧节点，若通过新节点更便宜则更新父节点

这保证了 RRT* 的“渐近最优”特性（迭代足够多时路径持续改进）。

### 4) 速度优化点

`src/rrt_star_region.cpp` 中使用了网格桶索引 `NodeGrid`：
- **nearest** 与 **radius 查询**先在桶内做局部检索，避免全量 O(N) 扫描
- 最近邻找不到时再回退线性扫描，保证鲁棒性

因此在迭代次数较大时，性能明显优于纯线性实现。

## 可视化与输出设计

### 1) 分割可视化

- C++ 直接写 `segmentation.png`、`segmentation_comparison.png`
- 同时写 `terrain_segmentation_data.txt`，并自动调用 Python 画高质量对比图（坐标轴、色条、后处理展示）

### 2) 路径可视化与状态切换

`rrt_star_region` 在 C++ 内：
- 直接输出 `rrt_star_path.png`（不依赖手动脚本）
- 打印路径中 rollable/lift 的切换 waypoint
- 生成 `waypoints_segmented.txt`，每行 `x y label`

`waypoints_segmented.txt` 是面向“更细节规划”的关键接口：下游可按 label 决定滚动控制模式或抬腿步态。

## 依赖

- CMake ≥ 3.18
- C++17 编译器
- CUDA Toolkit（与 trajopt_cpp 使用的 RBF JSON 格式兼容）

## 编译

```bash
cd terrain_friction_segment
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

可选：指定 CUDA 架构，例如 `-DCMAKE_CUDA_ARCHITECTURES=86`（按本机 GPU 修改）。

## 运行

```bash
# 不传参：使用与 trajopt_cpp/main.cpp:64 相同的默认 RBF 路径，输出到 test/terrain_res/
./terrain_friction_segment

# 指定 RBF、输出目录，可选摩擦系数（更小=更严，台阶更易识别）
./terrain_friction_segment /path/to/rbf.json test/terrain_res 0.35
```

- 第 1 个参数：RBF 参数 JSON（默认与 trajopt_cpp 一致）
- 第 2 个参数：**输出目录**（默认 `test/terrain_res`），所有分割与 RBF 副本写入该目录
- 第 3 个参数（可选）：摩擦系数 **μ**，默认 **0.3**

程序会：
- 将当前使用的 **RBF JSON 复制**到 `test/terrain_res/rbf.json`，供后续细规划使用
- 在 `test/terrain_res/` 下写入：`segmentation.txt`、`segmentation.png`、`terrain_segmentation_data.txt` 等
- 自动调用 Python 脚本做地形与分割对比

## 可视化

可视化由 **Python 脚本**完成（matplotlib），便于看清坐标轴、色标和对比：
- 运行可执行文件后会自动尝试调用 `scripts/visualize_terrain_segmentation.py`；若自动调用失败，可手动执行：
  ```bash
  python3 scripts/visualize_terrain_segmentation.py test/terrain_res/terrain_segmentation_data.txt
  ```
- 支持 `--auto-close 5`：显示 5 秒后自动关闭（与 trajopt_cpp 一致）。
- 仅查看分割 .txt 时仍可运行：`python3 scripts/visualize_segmentation.py segmentation.txt`

## 路径规划（RRT* + 区域代价）

**C++ 版（推荐，速度快）**：编译后运行 `rrt_star_region`，路径与带区域标签的 waypoint 写入 **test/terrain_res/**，便于后续细规划。

```bash
# 编译（与 terrain_friction_segment 一起）
cd build && make -j

# 默认：读 test/terrain_res/terrain_segmentation_data.txt，起点 (0.6, 1.2)，终点 (2.4, -3.0)，输出到 test/terrain_res/
./rrt_star_region

# 指定数据文件、起终点、路径输出文件
./rrt_star_region test/terrain_res/terrain_segmentation_data.txt 0.6 1.2 2.4 -3.0 test/terrain_res/rrt_star_path.txt
```

**test/terrain_res 输出**（供细规划使用）：
- `rbf.json`：当前地形 RBF 参数（由 terrain_friction_segment 写入）
- `terrain_segmentation_data.txt`：分割网格与地形高度（由 terrain_friction_segment 写入）
- `rrt_star_path.txt`：路径 waypoint 序列 (x y)
- `waypoints_segmented.txt`：**带区域标签的 waypoint**，每行 `x y label`（0=需抬腿，1=可滚动）
- `rrt_star_path.png`：路径可视化

代价 = 可滚动段长度 + k×抬腿段长度（默认 k=5）。C++ 实现使用网格空间索引加速。

### 区域代价 RRT* 详细原理

#### 1) 问题建模

目标是在二维平面上，从起点 `x_start` 到终点 `x_goal` 找一条最小代价路径。  
与普通最短路不同，本项目把地形语义（rollable/lift）编码进代价函数，而不是只看几何长度。

- 状态：`x = (x, y)`（连续）
- 地图：由分割得到标签场 `L(x)`，其中 `L=1` 为 rollable，`L=0` 为 lift
- 目标：最小化路径泛函

```text
J(path) = ∫ w(L(x(s))) ds
其中：
  w(1)=1
  w(0)=k   (k = lift_cost_weight, 默认 5)
```

直观上，算法会“偏好长一点但更可滚动”的路线，避免大段 lift 区域。

#### 2) 边代价离散化（代码里的实际做法）

在实现里，树上每条边 `u -> v` 的代价不是直接欧氏距离，而是做数值积分近似：

1. 把线段均匀分成 `edge_samples` 段  
2. 在每个小段中点查询 `label`  
3. 小段长度乘以对应权重后累加

```text
c(u,v) ≈ Σ_i  Δs * ( L_i==1 ? 1 : k )
```

这使代价函数对“局部穿越 lift 区域”非常敏感，且实现简单稳定。  
对应代码在 `src/rrt_star_region.cpp` 的 `SegData::edge_cost(...)`。

#### 3) RRT* 主循环（逐步解释）

每次迭代执行：

1. **采样 `x_rand`**  
   - 90%：在地图范围均匀采样  
   - 10%：直接采样终点（goal bias，加速收敛）

2. **找最近节点 `x_nearest`**  
   - 优先用 `NodeGrid` 桶索引局部搜索  
   - 若局部未命中，回退到全局线性扫描（鲁棒）

3. **扩展（steer）得到 `x_new`**  
   - 从 `x_nearest` 朝 `x_rand` 走固定步长 `step_length`

4. **选择最优父节点（ChooseParent）**  
   - 在 `rewire_radius` 邻域内枚举候选 `x_i`
   - 选择使 `cost(x_i)+c(x_i,x_new)` 最小的父节点

5. **插入新节点**

6. **重连（Rewire）**  
   - 对邻域内每个旧节点 `x_j`，若经 `x_new` 更便宜，则修改其父节点为 `x_new`

7. **终点连接判定**  
   - 若 `x_new` 落入 `goal_radius`，尝试连接终点并更新当前最优路径

这就是标准 RRT* 框架，只是把“欧氏边长”替换成“区域加权边代价”。

#### 4) 为什么这种代价适合你的任务

- `k` 大：强烈避开 lift 区，生成更“保守”的滚动主导路径
- `k` 小：允许更多穿越 lift 区，路径更短但抬腿动作增多
- 可直接把 `waypoints_segmented.txt` 交给下游控制器：
  - `label=1` 段：滚动模式
  - `label=0` 段：抬腿/跨越模式

等价于在“路径搜索阶段”提前做了模式分配的先验。

#### 5) 关键参数与行为

- `max_iter`：迭代次数；越大越可能找到更优路径，但耗时增加
- `step_length`：扩展步长；小则细腻但慢，大则快但可能抖动/绕行
- `rewire_radius`：重连半径；大则优化更充分但计算更重
- `edge_samples`：边积分采样数；大则代价估计更准但更慢
- `lift_cost_weight (k)`：rollable/lift 的相对惩罚强度（任务语义核心）

工程上常见调参顺序：先定 `k`（行为风格）-> 再调 `step_length`（几何分辨率）-> 最后调 `max_iter/rewire_radius`（质量-速度折中）。

#### 6) 复杂度与加速

朴素 RRT* 的近邻和半径查询通常接近 `O(N)`，总耗时会随节点数明显上升。  
本实现通过 `NodeGrid` 把大部分查询限制在局部桶中，平均代价显著下降，尤其在 `max_iter` 较大时效果明显。

#### 7) 与“全局最优”的关系

RRT* 理论上具有渐近最优性：迭代足够多时，路径会持续逼近最优。  
在本项目中，最优性的度量是“区域加权代价”而非纯长度，因此得到的是“动作成本意义下更优”的路径。

**Python 版**（备用）：`python3 scripts/rrt_star_region.py [data.txt] [sx sy gx gy]`  
仅画路径：`python3 scripts/plot_rrt_path.py <data.txt> <path.txt>`

## 文件说明

| 路径 | 说明 |
|------|------|
| `src/main.cpp` | 读 JSON、调 GPU、写结果 |
| `src/json_reader.cpp` | 读 RBF 的 cur_rbf_grid、rbf_weight（与 trajopt_cpp 一致） |
| `src/rbf_friction_cuda.cpp` | 主机端：分配/拷贝、调用 kernel、写文件 |
| `cuda/rbf_friction_kernel.cu` | GPU 内核：RBF 法向 + 摩擦锥判断 |
| `scripts/visualize_segmentation.py` | 读取 segmentation.txt 并绘图 |
| `scripts/visualize_terrain_segmentation.py` | 地形与分割对比（由可执行文件调用） |
| `scripts/rrt_star_region.py` | 带区域代价的 RRT*（Python 版） |
| `scripts/plot_rrt_path.py` | 读 C++ 输出的 path 文件并画图 |
| `src/rrt_star_region.cpp` | RRT* 核心（网格空间索引、区域代价） |
| 可执行文件 `rrt_star_region` | C++ 版 RRT*，输出路径并调 Python 画图 |

## 输出格式（test/terrain_res/）

- **rbf.json**：与输入一致的 RBF 参数 JSON（供细规划加载地形）
- **segmentation.txt**：第一行 `nx ny xmin xmax ymin ymax`，第二行 `nx*ny` 个 0/1（0=需抬腿，1=可滚动）
- **terrain_segmentation_data.txt**：含网格与地形高度，供可视化与 RRT* 读取
- **rrt_star_path.txt**：路径 waypoint，每行 `x y`
- **waypoints_segmented.txt**：带区域标签的 waypoint，每行 `x y label`（0=需抬腿，1=可滚动），供细规划使用
