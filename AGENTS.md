# AGENTS.md

本文件适用于仓库根目录及全部子目录。修改某个子工程前，应先阅读该子工程的 `README.md`、`CMakeLists.txt` 和入口文件。仓库中的 `build*`、日志、图片及实验输出可能来自旧版本；理解当前实现时，以源码和当前 CMake 配置为准。

## 项目目的

本项目研究轮腿式机器人/双轮机构在 RBF 高度场上的轨迹规划与非线性优化，主要目标是：

1. 从 RBF JSON 重建地形高度、梯度和表面法向；
2. 根据摩擦条件将地形划分为可滚动区和需抬腿区；
3. 使用带区域代价的 RRT* 生成二维粗路径；
4. 使用 Chebyshev 多项式参数化基座及左右轮轨迹；
5. 使用 Ipopt 求解包含边界、几何、接触、摩擦、平衡和轨迹平滑约束的 NLP；
6. 输出轨迹、接触力、分区数据和可视化结果。

区域标签的固定语义为：`0 = lift/需抬腿`，`1 = rollable/可滚动`。不得在没有迁移全部读写端的情况下修改该约定。

### RBF 地形理论与适用边界

根目录 `RBF_LIO.pdf` 是当前 RBF 地形表示的理论参考。论文把局部地形建模为世界坐标系紧致 ROI 内的单值 2.5D 高度场：

```text
z = f(x, y) = sum_i w_i * exp(-||[x,y] - c_i||^2 / (2*sigma^2))
```

- `c_i` 是水平面内的二维 RBF 中心，`w_i` 是与中心一一对应的标量权重，`sigma` 是高斯核带宽；
- 高度梯度可解析求得，向上的未归一化表面法向约定为 `[-df/dx, -df/dy, 1]`；
- 论文上游通过网格候选、自适应中心选择、岭回归初值和活跃中心递归更新在线构建地形；本仓库没有实现完整的 RBF-LIO 里程计或在线拟合，只消费导出的中心和权重进行分区、规划与优化；
- 论文实验参数为中心网格分辨率 `0.07 m`、`sigma = 0.04 m`，而本仓库 C++ `JSONReader` 当前硬编码 `sigma = 0.14`。二者不一致，不能把论文参数直接覆盖到现有工程，但目前参数先全部以代码为准。
- 模型假设局部地形能写成唯一的 `z=f(x,y)`，并在轮地约束中假设稳定接触且每个轮子只有一个有效接触点。陡直面、悬垂/多层表面、尖锐台阶边缘、轮滑、失联和强振动会削弱该假设；相关约束应允许可靠性判定或降权。

`terrain_friction_segment` 的坡度/摩擦分区和各轨迹优化器是 RBF 地形的下游应用，不等同于论文中的 RBF-LIO 状态估计或完整轮地流形残差。文档和代码注释必须保持这一区分。

## 软件架构

工程包含两条实现路线：

- Python/CasADi 原型：`solve_casadi.py`、`trajopt_sdf.py`；
- C++/Ipopt 主线：`terrain_friction_segment`、`seg_trajopt_cpp`、`lift_leg_trajopt_cpp`、`trajopt_cpp`。

主数据流如下：

```text
RBF JSON
  -> terrain_friction_segment（CUDA 地形分区）
  -> terrain_segmentation_data.txt
  -> rrt_star_region（带区域代价的 RRT*）
  -> waypoints_segmented.txt
  -> seg_trajopt_cpp（可滚动段细化）
     或 lift_leg_trajopt_cpp（抬腿段细化）
  -> Ipopt
  -> initial_trajectory.txt / final_trajectory.txt / 可视化结果
```

`trajopt_cpp` 是独立的通用全轨迹优化器，同时优化 10 组轨迹系数和 6 组离散接触力。`seg_trajopt_cpp` 和 `lift_leg_trajopt_cpp` 会在 CMake 中直接复用 `trajopt_cpp` 的 Chebyshev 基、RBF 地形、JSON 读取和部分目标函数源码，因此必须保持这三个目录的同级布局。

目前细规划器只处理 `waypoints_segmented.txt` 中第一段匹配类型的连续路径；工程尚无自动遍历全部区段、拼接轨迹及处理段间连续性的总调度器。不要假设单次运行会生成完整端到端轨迹。

## 编译方法

### 依赖

- CMake 3.10 以上；CUDA 子工程要求 CMake 3.18 以上；
- 支持 C++17 的编译器；
- Ipopt 开发库，通常通过 `pkg-config` 查找；当前配置默认使用 MUMPS；
- CUDA Toolkit，仅 `terrain_friction_segment` 的分区程序需要；
- Python 3、NumPy、Matplotlib；Python NLP 原型还需要 CasADi，GPU 性能实验需要 CuPy。

仓库没有统一依赖锁文件。不要在未获任务授权时升级依赖、修改系统环境或引入新的包管理方案。

### C++ 子工程

从仓库根目录执行：

```bash
cmake -S trajopt_cpp -B trajopt_cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build trajopt_cpp/build --parallel

cmake -S terrain_friction_segment -B terrain_friction_segment/build -DCMAKE_BUILD_TYPE=Release
cmake --build terrain_friction_segment/build --parallel

cmake -S seg_trajopt_cpp -B seg_trajopt_cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build seg_trajopt_cpp/build --parallel

cmake -S lift_leg_trajopt_cpp -B lift_leg_trajopt_cpp/build_lift_new -DCMAKE_BUILD_TYPE=Release
cmake --build lift_leg_trajopt_cpp/build_lift_new --parallel
```

多配置生成器可能把可执行文件放在 `build/Release/`，Windows 下文件名还会带 `.exe`。不要手工编辑 CMake 生成文件或复用明显来自其他机器的缓存。

## 运行方法

默认路径中仍有 `/home/yizhe/...` 和 `/tmp/...` 等 Linux 绝对路径。运行时应尽量显式传入输入、输出路径；不要依赖默认路径。

### 完整 C++ 数据链

```bash
# 1. RBF 地形摩擦分区；第三个参数为 mu
./terrain_friction_segment/build/terrain_friction_segment \
  test/terrain_res/rbf.json test/terrain_res 0.3

# 2. 带区域代价的 RRT* 粗规划
./terrain_friction_segment/build/rrt_star_region \
  test/terrain_res/terrain_segmentation_data.txt \
  0.6 1.2 2.4 -3.0 test/terrain_res/rrt_star_path.txt

# 3a. 优化第一段可滚动路径
./seg_trajopt_cpp/build/seg_trajopt_cpp \
  test/terrain_res seg_trajopt_cpp/seg_out

# 3b. 优化第一段需抬腿路径
./lift_leg_trajopt_cpp/build_lift_new/lift_leg_trajopt_cpp \
  test/terrain_res lift_leg_trajopt_cpp/lift_out
```

注意：`rrt_star_region` 当前只允许通过参数覆盖分区输入和路径文本输出，带标签 waypoint 与 PNG 的默认输出目录仍在源码中硬编码。修改运行路径相关逻辑时，要同时检查 `path_file`、`waypoints_file` 和 `png_file`。

### 通用全轨迹优化器

```bash
./trajopt_cpp/build/trajopt_cpp test/terrain_res/rbf.json
```

### Python 原型

```bash
python solve_casadi.py
python trajopt_sdf.py
```

这两个脚本当前在 `__main__` 中包含硬编码 RBF 路径。除非任务明确要求，不要把本地机器路径提交到源码；新增运行方式应优先使用命令行参数或配置对象。

### 可视化

```bash
python seg_trajopt_cpp/scripts/plot_seg_trajectory.py seg_trajopt_cpp/seg_out --no-show
python lift_leg_trajopt_cpp/scripts/plot_seg_trajectory.py \
  lift_leg_trajopt_cpp/lift_out --rbf test/terrain_res/rbf.json --no-show
python terrain_friction_segment/scripts/visualize_terrain_segmentation.py \
  test/terrain_res/terrain_segmentation_data.txt
```

自动化或无图形环境中必须使用无窗口选项，避免 `plt.show()`、`xdg-open` 或 `std::system()` 阻塞任务。

## 关键模块说明

### Python 模块

- `core/variable.py`：`VarBlock`、`Variables`，管理变量块、边界、缩放以及求解器向量的打包/解包；
- `core/objective.py`：目标项接口、加权组合及终点、平滑、侧滑、轮距、离地等目标；
- `constraints/`：离地间隙、轮距和防交叉约束；
- `solver/base.py`：`NLPProblem`、`Solver` 协议；
- `solver/casadi_ipopt.py`：把 NumPy 目标/约束回调包装成 CasADi callback；
- `geometry/terrain.py`：简化单 RBF 地形演示，不等同于 C++ 多中心 RBF 实现；
- `trajopt_sdf.py`：含显式轮高和接触力的单文件 CasADi 原型。

模块化 Python 层目前缺少将 `Variables + Objective + Constraints` 组装为具体 `NLPProblem` 的实现。不要把根目录 `main.py` 当成完整求解入口。

### C++ 公共层

- `BasisFunction` / `ChebyshevBasis`：轨迹基函数和一、二阶导数；参数 `s` 固定为 `[0,1]`，内部映射到 Chebyshev 区间 `[-1,1]`；
- `RBFTerrain`：RBF 地形高度、梯度、法向和边界；
- `ConstraintBase`：约束值、Jacobian、上下界和稀疏结构接口；
- `ObjectiveBase`：目标值和梯度累加接口；
- `TrajOptNLP` / `SegTrajOptNLP`：聚合目标与约束并实现 Ipopt `TNLP` 回调；
- `JSONReader`：读取 `cur_rbf_grid` 和 `rbf_weight`，二者的索引必须一一对应。当前 `sigma` 固定为 `0.14`，且与 `RBF_LIO.pdf` 的实验值 `0.04 m` 不同；修改带宽或数据格式前必须确认数据生产端并同步所有 Python/C++ 读写端。

### 规划与细化模块

- `terrain_friction_segment/cuda/rbf_friction_kernel.cu`：GPU 地形法向和摩擦可行性判断；
- `terrain_friction_segment/src/rrt_star_region.cpp`：RRT*、区域边代价和空间桶索引；
- `seg_trajopt_cpp`：8 组 Chebyshev 系数，左右轮 z 由地形确定；
- `lift_leg_trajopt_cpp`：左右轮 z 使用分段系数，包含接触掩码、穿透力近似和摆动轮抬升目标；
- `trajopt_cpp`：10 组轨迹系数加 6 组接触力，包含几何、接触、力/力矩平衡和摩擦目标。

## 不允许修改的位置

除非用户明确要求，以下内容不得修改：

1. `.git/` 及任何 Git 内部文件；
2. `**/build*/`、`**/__pycache__/` 和其中的 CMake 缓存、目标文件、依赖文件、二进制及字节码；
3. `*.o`、`*.d`、`*.cubin`、`*.pyc`、已编译可执行文件；
4. `terrain_friction_segment/external/stb_image_write.h`；
5. `ipopt_log.txt`、`presentation.log`、`result` 以及现有 `*_out*`、`seg_out`、实验图片和调试 CSV；
6. 用户未明确要求更新时的 `test/terrain_res/` 测试夹具及其中的 RBF、分区、路径和 waypoint 数据。
7. 根目录 `RBF_LIO.pdf`；它是理论参考原件，不得编辑、覆盖、压缩或重新导出。

需要验证输出时，应写入新建的临时构建/结果目录，不得覆盖现有实验结果。不要使用清理命令批量删除仓库中的历史构建或输出目录。

## 编码规范

### 通用规则

- 只修改任务需要的文件，保留工作区中的既有用户改动；
- 新增文本使用 UTF-8，避免改变无关文件的编码或换行风格；
- 文件路径必须通过参数、配置或 `std::filesystem` 传递，禁止新增机器专属绝对路径；
- 物理量命名或注释应说明单位；摩擦、质量、重力、轮距、间隙等参数不得以无说明的魔数散落；
- 不得静默改变 RBF JSON、分区文件、轨迹文件和 waypoint 文件的数据格式；
- 修改公共数据契约时，必须同步更新生产者、消费者、文档和测试。
- RBF 数据契约还包括坐标系、长度单位、中心/权重顺序、核函数中的 `1/(2*sigma^2)` 约定、带宽单位以及法向方向；任何一项都不得只在单个模块中修改；
- 新增论文中的轮地约束时，应显式记录接触点、轮轴、地形法向所在坐标系，并为打滑、失联或接触歧义提供检测、禁用或降权机制。

### C++

- 使用 C++17；优先 RAII、`std::vector`、`std::shared_ptr`/`std::unique_ptr` 和 `const`；
- 避免新增裸 `new/delete`、不受控全局状态和 C 风格数组；Ipopt `SmartPtr` 场景除外；
- 轨迹变量顺序属于核心 ABI。修改 `cbx/cby/cbz/cpsi/clx/cly/clz/crx/cry/crz` 或接触力偏移时，必须同步全部约束、目标、初值和提取函数；
- 新增约束必须保证 `getNumConstraints()`、`eval_g()`、`getBounds()`、`eval_jac_g()` 和 `getJacobianStructure()` 的维度及顺序完全一致；
- `ObjectiveBase::eval_grad_f()` 必须向已有梯度累加，不得清零或覆盖其他目标项贡献；
- 优先提供解析梯度/Jacobian，并用有限差分校验；不要依赖 Ipopt 对黑盒函数自动求导；
- 对开方、除法、归一化和指数计算保留明确的数值保护，并检查 NaN/Inf；
- CMake 中公共源码的复用关系应保持显式，不要复制一份公共实现到新的子工程。

### Python

- 可复用模块遵循 PEP 8，公共函数和类尽量添加类型标注及简短 docstring；
- 优先 NumPy 向量化，避免在高频目标/约束回调中创建不必要的大型临时数组；
- 可复用功能放入 `core/`、`constraints/`、`geometry/` 或 `solver/`，不要继续扩大单文件原型中的重复实现；
- 绘图函数必须使用传入对象，不依赖仅在 `__main__` 中存在的全局变量；
- CLI 脚本应支持显式输入/输出参数，并提供适用于 CI 的无窗口模式。

### 验证要求

- 变量布局、目标或约束变化：至少进行梯度/Jacobian 有限差分检查；
- C++ NLP 变化：编译对应子工程，并用小规模输入做 Ipopt 冒烟测试；
- 分区/RRT* 变化：验证标签语义、文件尺寸、起终点和 waypoint 输出；
- 数据链变化：至少验证 `RBF JSON -> segmentation -> RRT* -> seg/lift` 的相关接口；
- 不要把新生成的构建产物、日志、图片或大数据文件提交为源码修改。
