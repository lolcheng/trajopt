# TrajOpt C++ 项目

基于Ipopt的轨迹优化C++实现，避免符号图构建，提升性能。

## 项目结构

```
trajopt_cpp/
├── include/              # 头文件
│   ├── trajopt_nlp.hpp  # 主NLP问题类
│   ├── constraint_base.hpp
│   ├── rbf_terrain.hpp
│   ├── polynomial.hpp
│   └── constraints/     # 各约束模块
├── src/                  # 源文件
│   ├── trajopt_nlp.cpp
│   ├── rbf_terrain.cpp
│   ├── polynomial.cpp
│   └── constraints/
├── build/                # 编译目录
├── CMakeLists.txt
├── main.cpp
├── CONSTRAINTS.md        # 约束总结文档
└── ARCHITECTURE.md       # 架构设计文档
```

## 编译

```bash
cd trajopt_cpp
mkdir -p build
cd build
cmake ..
make
```

## 运行

```bash
./trajopt_cpp
```

## 当前版本

当前版本实现了：
- ✅ RBF地形计算
- ✅ 多项式轨迹表示
- ✅ 边界约束（使用RBF计算起点终点高度，包括轮子初末点约束）
- ✅ 轮子-中心距离约束（硬约束，2×n_sample个不等式约束）
- ✅ 左右轮距离约束（硬约束，n_sample个不等式约束）
- ✅ 防交叉约束（硬约束，n_sample个不等式约束）
- ✅ 从JSON读取RBF参数
- ✅ RBF地形可视化
- ✅ 轨迹可视化（求解成功后）
- ✅ 平滑项（目标函数项，对轨迹二阶差分惩罚）
- ✅ 正则化项（目标函数项，防止系数过大）
- ✅ 轮子-中心距离标称项（目标函数项）
- ✅ 左右轮距离标称项（目标函数项）
- ✅ 轮子-地形距离项（目标函数项）
- ✅ 偏航-路径对齐项（目标函数项）
- ✅ 法向力上限约束（硬约束，50个不等式约束）
- ✅ 接触刚度项（目标函数项）
- ✅ 力平衡项（目标函数项）
- ✅ 力矩平衡项（目标函数项）
- ✅ 摩擦锥违反项（目标函数项）
- ✅ 优化时间统计

## 可视化

### 地形可视化

程序在求解前会自动生成RBF地形数据并尝试可视化。如果自动可视化失败，可以手动运行：

```bash
# 从项目根目录
python3 trajopt_cpp/visualize_terrain.py /tmp/rbf_terrain_data.txt
```

地形可视化脚本会显示：
- 3D表面图：展示RBF地形的3D形状
- 2D等高线图：展示地形的俯视图

### 轨迹可视化

程序在求解成功后会自动生成轨迹数据并尝试可视化。如果自动可视化失败，可以手动运行：

```bash
# 从项目根目录
python3 trajopt_cpp/visualize_trajectory.py /tmp/trajectory_data.txt /tmp/rbf_terrain_data.txt
```

轨迹可视化脚本会显示：
- 3D轨迹图：在3D地形上显示基座轨迹、左右轮轨迹，以及偏航角方向
- 2D俯视图：显示轨迹的俯视图，包含地形等高线和偏航角箭头

### 自动关闭功能

可视化窗口支持自动关闭功能，可以在`main.cpp`中配置：

```cpp
// 可视化设置
const bool VISUALIZATION_AUTO_CLOSE = true;  // 是否自动关闭可视化窗口
const double VISUALIZATION_CLOSE_TIME = 5.0;  // 自动关闭时间（秒）
```

- 当`VISUALIZATION_AUTO_CLOSE = true`时，可视化窗口会在指定时间（默认5秒）后自动关闭
- 当`VISUALIZATION_AUTO_CLOSE = false`时，可视化窗口会保持打开，需要手动关闭

手动运行可视化脚本时，也可以使用`--auto-close`参数：

```bash
# 使用默认5秒自动关闭
python3 trajopt_cpp/visualize_terrain.py /tmp/rbf_terrain_data.txt --auto-close

# 自定义关闭时间（例如10秒）
python3 trajopt_cpp/visualize_terrain.py /tmp/rbf_terrain_data.txt --auto-close 10.0
```

## 实现完整性

### 与Python参考实现对比

C++实现已完整实现了Python参考版本（`trajopt_sdf.py`）中的所有约束和目标函数项：

**变量**:
- ✅ 多项式系数：50个（与Python一致）
- ✅ 接触力：150个（Python为600个，差异仅在于采样点数：C++使用25个，Python使用100个）

**硬约束**:
- ✅ Base边界约束：8个（与Python一致）
- ✅ 轮子初末点约束：12个（C++增强，Python无此约束）
- ✅ 轮子-中心距离约束：50个（Python为200个，差异仅在于采样点数）
- ✅ 左右轮距离约束：25个（Python为100个，差异仅在于采样点数）
- ✅ 防交叉约束：25个（Python为100个，差异仅在于采样点数）
- ✅ 法向力上限约束：50个（Python为200个，差异仅在于采样点数）

**软约束（目标函数项）**:
- ✅ 平滑项
- ✅ 正则化项
- ✅ 轮子-中心距离标称项
- ✅ 左右轮距离标称项
- ✅ 轮子-地形距离项
- ✅ 偏航-路径对齐项
- ✅ 接触刚度项
- ✅ 力平衡项
- ✅ 力矩平衡项
- ✅ 摩擦锥违反项

**主要区别**:
1. **采样点数**: C++使用`n_sample=25`（用户要求降低复杂度，间距约20cm），Python使用`n_sample=100`（间距约5cm）
2. **边界约束增强**: C++增加了轮子初末点约束，确保轮子在起点和终点触地且在base附近，这是对Python版本的改进

**结论**: C++实现功能完整，所有约束和目标函数项均已实现。采样点数的差异是设计选择，不影响功能完整性。

## 依赖

- Ipopt (需要安装并配置)
- CMake >= 3.10
- C++17 编译器
- Python3 + matplotlib + numpy (用于可视化)

