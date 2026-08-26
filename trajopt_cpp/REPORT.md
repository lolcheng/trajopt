# 轨迹优化项目技术报告

## 1. 问题概述

本项目实现了一个基于非线性优化的轨迹规划系统，用于在复杂地形上规划两轮机器人的运动轨迹。系统使用多项式参数化轨迹，通过RBF（径向基函数）表示地形，并考虑多种物理约束和优化目标。

### 1.1 核心问题

给定：
- 起点和终点的位置和姿态
- RBF地形模型（通过JSON文件加载）
- 物理参数（质量、重力、摩擦系数等）

求解：
- 基座（base）轨迹：位置 $(x_b(s), y_b(s), z_b(s))$ 和偏航角 $\psi(s)$
- 左右轮轨迹：位置 $(x_l(s), y_l(s), z_l(s))$ 和 $(x_r(s), y_r(s), z_r(s))$
- 接触力：左右轮在各采样点的接触力 $(f_{Lx}, f_{Ly}, f_{Lz})$ 和 $(f_{Rx}, f_{Ry}, f_{Rz})$

其中 $s \in [0,1]$ 是归一化的路径参数。

### 1.2 优化变量

- **多项式系数**（50个）：
  - 基座：$c_{bx}, c_{by}, c_{bz}, c_{\psi}$（各5个系数）
  - 左轮：$c_{lx}, c_{ly}, c_{lz}$（各5个系数）
  - 右轮：$c_{rx}, c_{ry}, c_{rz}$（各5个系数）

- **接触力**（150个）：
  - 左轮：$f_{Lx}, f_{Ly}, f_{Lz}$（各25个采样点）
  - 右轮：$f_{Rx}, f_{Ry}, f_{Rz}$（各25个采样点）

**总变量数**：200个

### 1.3 轨迹表示

使用5阶多项式基函数表示轨迹：
$$x(s) = \sum_{k=0}^{4} c_k s^k$$

在 $n_{\text{sample}} = 25$ 个采样点上评估约束和目标函数。

## 2. 代码结构

### 2.1 目录结构

```
trajopt_cpp/
├── include/                    # 头文件
│   ├── trajopt_nlp.hpp        # 主NLP问题类
│   ├── constraint_base.hpp    # 约束基类
│   ├── objective_base.hpp     # 目标函数基类
│   ├── rbf_terrain.hpp        # RBF地形类
│   ├── polynomial.hpp         # 多项式工具
│   ├── constraints/           # 约束模块
│   │   ├── boundary_constraint.hpp
│   │   ├── wheel_base_distance_constraint.hpp
│   │   ├── wheels_distance_constraint.hpp
│   │   ├── anti_crossing_constraint.hpp
│   │   ├── normal_force_upper_bound_constraint.hpp
│   │   └── wheel_airborne_constraint.hpp
│   └── objectives/            # 目标函数模块
│       ├── smoothing_objective.hpp
│       ├── regularization_objective.hpp
│       ├── wheel_base_distance_nominal_objective.hpp
│       ├── wheels_distance_nominal_objective.hpp
│       ├── wheel_terrain_distance_objective.hpp
│       ├── wheel_airborne_penalty_objective.hpp
│       ├── path_length_objective.hpp
│       ├── yaw_path_alignment_objective.hpp
│       ├── contact_stiffness_objective.hpp
│       ├── force_balance_objective.hpp
│       ├── torque_balance_objective.hpp
│       └── friction_cone_objective.hpp
├── src/                       # 源文件（对应include结构）
├── main.cpp                   # 主程序入口
└── CMakeLists.txt            # 构建配置
```

### 2.2 核心类设计

#### 2.2.1 ConstraintBase（约束基类）

所有约束继承自 `ConstraintBase`，实现统一接口：

```cpp
class ConstraintBase {
    virtual int getNumConstraints() const = 0;
    virtual void eval_g(const double* x, double* g) = 0;
    virtual void eval_jac_g(const double* x, double* jac_g, int offset = 0) = 0;
    virtual void getBounds(double* g_lb, double* g_ub, int offset = 0) = 0;
    virtual int getJacobianStructure(int* iRow, int* jCol, int offset = 0) = 0;
};
```

#### 2.2.2 ObjectiveBase（目标函数基类）

所有目标函数项继承自 `ObjectiveBase`：

```cpp
class ObjectiveBase {
    virtual double eval_f(const double* x) = 0;
    virtual void eval_grad_f(const double* x, double* grad_f) = 0;
    virtual bool isEnabled() const = 0;
};
```

#### 2.2.3 TrajOptNLP（主NLP问题类）

实现Ipopt的 `TNLP` 接口，管理所有约束和目标函数：

```cpp
class TrajOptNLP : public Ipopt::TNLP {
    void addConstraint(std::shared_ptr<ConstraintBase> constraint);
    void addObjective(std::shared_ptr<ObjectiveBase> objective);
    void setTerrain(std::shared_ptr<RBFTerrain> terrain);
    // Ipopt接口实现...
};
```

### 2.3 变量布局

变量按以下顺序排列：
1. 多项式系数（50个）：`cbx, cby, cbz, cpsi, clx, cly, clz, crx, cry, crz`
2. 接触力（150个）：`fLx, fLy, fLz, fRx, fRy, fRz`

## 3. 约束项详细说明

### 3.1 硬约束（Hard Constraints）

硬约束通过 `ConstraintBase` 接口实现，在优化过程中必须严格满足。

#### 3.1.1 边界约束（BoundaryConstraint）

**定义**：约束起点和终点的位置和姿态。

**约束数量**：20个等式约束
- Base边界：8个（起点和终点的 $x_b, y_b, z_b, \psi$）
- 轮子初末点：12个（左右轮在起点和终点的 $x, y, z$ 位置）

**数学表达式**：

起点 ($s=0$)：
- $x_b(0) = x_0$
- $y_b(0) = y_0$
- $z_b(0) = h(x_0, y_0) + \text{CLEAR}_Z$
- $\psi(0) = \psi_0$
- 轮子位置：基于base位置和yaw角计算，确保触地且在base附近

终点 ($s=1$)：
- $x_b(1) = x_1$
- $y_b(1) = y_1$
- $z_b(1) = h(x_1, y_1) + \text{CLEAR}_Z$
- $\psi(1) = \psi_1$
- 轮子位置：基于base位置和yaw角计算

**约束值**：
$$g_i = \Phi(s) \cdot c - \text{target}$$

其中 $\Phi(s)$ 是多项式基矩阵在 $s=0$ 或 $s=1$ 的值。

**雅可比矩阵**：
$$\frac{\partial g_i}{\partial c_j} = \Phi_{ij}(s)$$

**海塞矩阵**：零矩阵（线性约束）

**代码实现说明**：边界约束是线性约束，海塞矩阵为零。代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.1.2 轮子-中心距离约束（WheelBaseDistanceConstraint）

**定义**：约束轮子与基座中心的3D欧氏距离在合理范围内。

**约束数量**：$2 \times n_{\text{sample}} = 50$ 个不等式约束

**数学表达式**：

左轮：
$$\text{WHEEL\_BASE\_MIN}^2 \leq d_L^2 \leq \text{WHEEL\_BASE\_MAX}^2$$

右轮：
$$\text{WHEEL\_BASE\_MIN}^2 \leq d_R^2 \leq \text{WHEEL\_BASE\_MAX}^2$$

其中：
$$d_L^2 = (x_b - x_l)^2 + (y_b - y_l)^2 + (z_b - z_l)^2$$
$$d_R^2 = (x_b - x_r)^2 + (y_b - y_r)^2 + (z_b - z_r)^2$$

**参数**：
- $\text{WHEEL\_BASE\_MIN} = 0.4$ m
- $\text{WHEEL\_BASE\_MAX} = 1.2$ m

**约束值**：
$$g_i = d^2$$

**约束边界**：
- 下界：$g_{\text{lb}} = \text{WHEEL\_BASE\_MIN}^2$
- 上界：$g_{\text{ub}} = \text{WHEEL\_BASE\_MAX}^2$

**雅可比矩阵**：

对左轮距离平方：
$$\frac{\partial d_L^2}{\partial c_{bx}} = 2(x_b - x_l) \cdot \Phi_{bx}(s)$$
$$\frac{\partial d_L^2}{\partial c_{by}} = 2(y_b - y_l) \cdot \Phi_{by}(s)$$
$$\frac{\partial d_L^2}{\partial c_{bz}} = 2(z_b - z_l) \cdot \Phi_{bz}(s)$$
$$\frac{\partial d_L^2}{\partial c_{lx}} = -2(x_b - x_l) \cdot \Phi_{lx}(s)$$
$$\frac{\partial d_L^2}{\partial c_{ly}} = -2(y_b - y_l) \cdot \Phi_{ly}(s)$$
$$\frac{\partial d_L^2}{\partial c_{lz}} = -2(z_b - z_l) \cdot \Phi_{lz}(s)$$

右轮类似。

**海塞矩阵**：

理论上，距离平方的海塞矩阵为：
$$\frac{\partial^2 d_L^2}{\partial c_{bx}^2} = 2 \Phi_{bx}(s)^2$$
$$\frac{\partial^2 d_L^2}{\partial c_{bx} \partial c_{lx}} = -2 \Phi_{bx}(s) \Phi_{lx}(s)$$
$$\frac{\partial^2 d_L^2}{\partial c_{lx}^2} = 2 \Phi_{lx}(s)^2$$

其他项类似。

**代码实现说明**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.1.3 左右轮距离约束（WheelsDistanceConstraint）

**定义**：约束左右轮之间的3D欧氏距离在合理范围内。

**约束数量**：$n_{\text{sample}} = 25$ 个不等式约束

**数学表达式**：

$$\text{WHEELS\_MIN}^2 \leq d_{LR}^2 \leq \text{WHEELS\_MAX}^2$$

其中：
$$d_{LR}^2 = (x_r - x_l)^2 + (y_r - y_l)^2 + (z_r - z_l)^2$$

**参数**：
- $\text{WHEELS\_MIN} = 0.1$ m
- $\text{WHEELS\_MAX} = 0.8$ m

**约束值**：
$$g_i = d_{LR}^2$$

**约束边界**：
- 下界：$g_{\text{lb}} = \text{WHEELS\_MIN}^2$
- 上界：$g_{\text{ub}} = \text{WHEELS\_MAX}^2$

**雅可比矩阵**：

$$\frac{\partial d_{LR}^2}{\partial c_{lx}} = -2(x_r - x_l) \cdot \Phi_{lx}(s)$$
$$\frac{\partial d_{LR}^2}{\partial c_{ly}} = -2(y_r - y_l) \cdot \Phi_{ly}(s)$$
$$\frac{\partial d_{LR}^2}{\partial c_{lz}} = -2(z_r - z_l) \cdot \Phi_{lz}(s)$$
$$\frac{\partial d_{LR}^2}{\partial c_{rx}} = 2(x_r - x_l) \cdot \Phi_{rx}(s)$$
$$\frac{\partial d_{LR}^2}{\partial c_{ry}} = 2(y_r - y_l) \cdot \Phi_{ry}(s)$$
$$\frac{\partial d_{LR}^2}{\partial c_{rz}} = 2(z_r - z_l) \cdot \Phi_{rz}(s)$$

**海塞矩阵**：类似轮子-中心距离约束。

**代码实现说明**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.1.4 防交叉约束（AntiCrossingConstraint）

**定义**：防止左右轮轨迹在base坐标系下交叉。

**约束数量**：$n_{\text{sample}} = 25$ 个不等式约束

**数学表达式**：

在base坐标系下，左轮的横向坐标应大于右轮的横向坐标：
$$y_{L,\text{base}} \geq y_{R,\text{base}}$$

其中：
- 横向单位向量：$\mathbf{u}_{\text{lat}} = [-\sin\psi, \cos\psi]^T$
- 左轮横向坐标：$y_{L,\text{base}} = (\mathbf{p}_L - \mathbf{p}_b) \cdot \mathbf{u}_{\text{lat}}$
- 右轮横向坐标：$y_{R,\text{base}} = (\mathbf{p}_R - \mathbf{p}_b) \cdot \mathbf{u}_{\text{lat}}$

**约束值**：
$$g_i = y_{L,\text{base}} - y_{R,\text{base}}$$

**约束边界**：
- 下界：$g_{\text{lb}} = 0$
- 上界：$g_{\text{ub}} = +\infty$

**雅可比矩阵**：

详细推导：
$$g = y_{L,\text{base}} - y_{R,\text{base}} = (\mathbf{p}_L - \mathbf{p}_b) \cdot \mathbf{u}_{\text{lat}} - (\mathbf{p}_R - \mathbf{p}_b) \cdot \mathbf{u}_{\text{lat}}$$

对 $c_{bx}$ 的导数：
$$\frac{\partial g}{\partial c_{bx}} = \frac{\partial}{\partial c_{bx}}[(\mathbf{p}_L - \mathbf{p}_b) \cdot \mathbf{u}_{\text{lat}}] - \frac{\partial}{\partial c_{bx}}[(\mathbf{p}_R - \mathbf{p}_b) \cdot \mathbf{u}_{\text{lat}}]$$
$$= -u_{\text{lat},x} \cdot \Phi_{bx}(s) + u_{\text{lat},x} \cdot \Phi_{bx}(s) = 0$$

实际上，由于 $\mathbf{p}_L$ 和 $\mathbf{p}_R$ 不依赖于 $c_{bx}$，只有 $\mathbf{p}_b$ 依赖于 $c_{bx}$，所以：
$$\frac{\partial g}{\partial c_{bx}} = -u_{\text{lat},x} \cdot \Phi_{bx}(s)$$

类似地：
$$\frac{\partial g}{\partial c_{by}} = -u_{\text{lat},y} \cdot \Phi_{by}(s)$$

对 $c_{\psi}$ 的导数（需要考虑 $\mathbf{u}_{\text{lat}}$ 对 $\psi$ 的依赖）：
$$\frac{\partial g}{\partial c_{\psi}} = (\mathbf{p}_L - \mathbf{p}_b) \cdot \frac{\partial \mathbf{u}_{\text{lat}}}{\partial \psi} \cdot \Phi_{\psi}(s) - (\mathbf{p}_R - \mathbf{p}_b) \cdot \frac{\partial \mathbf{u}_{\text{lat}}}{\partial \psi} \cdot \Phi_{\psi}(s)$$

其中：
$$\frac{\partial \mathbf{u}_{\text{lat}}}{\partial \psi} = [-\cos\psi, -\sin\psi]^T$$

所以：
$$\frac{\partial g}{\partial c_{\psi}} = [-(v_{Lx} - v_{Rx}) \cos\psi - (v_{Ly} - v_{Ry}) \sin\psi] \cdot \Phi_{\psi}(s)$$

其中 $v_{Lx} = x_l - x_b$，$v_{Ly} = y_l - y_b$，$v_{Rx} = x_r - x_b$，$v_{Ry} = y_r - y_b$。

对轮子位置的导数：
$$\frac{\partial g}{\partial c_{lx}} = u_{\text{lat},x} \cdot \Phi_{lx}(s)$$
$$\frac{\partial g}{\partial c_{ly}} = u_{\text{lat},y} \cdot \Phi_{ly}(s)$$
$$\frac{\partial g}{\partial c_{rx}} = -u_{\text{lat},x} \cdot \Phi_{rx}(s)$$
$$\frac{\partial g}{\partial c_{ry}} = -u_{\text{lat},y} \cdot \Phi_{ry}(s)$$

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.1.5 法向力上限约束（NormalForceUpperBoundConstraint）

**定义**：约束法向力不超过物理上限。

**约束数量**：$2 \times n_{\text{sample}} = 50$ 个不等式约束

**数学表达式**：

左轮：
$$f_{Ln} \leq f_{N,\max} = m \cdot (g + a_{\max})$$

右轮：
$$f_{Rn} \leq f_{N,\max}$$

其中：
- $f_{Ln} = \mathbf{f}_L \cdot \mathbf{n}_L$（法向力）
- $f_{Rn} = \mathbf{f}_R \cdot \mathbf{n}_R$
- $\mathbf{n}_L, \mathbf{n}_R$ 是地形表面法向量

**参数**：
- $m = 1.0$ kg
- $g = 9.81$ m/s²
- $a_{\max} = 1.0$ m/s²
- $f_{N,\max} = 1.0 \times (9.81 + 1.0) = 10.81$ N

**约束值**：
$$g_i = f_n$$

**约束边界**：
- 下界：$g_{\text{lb}} = -\infty$
- 上界：$g_{\text{ub}} = f_{N,\max}$

**雅可比矩阵**：

对接触力的导数（直接）：
$$\frac{\partial f_n}{\partial f_x} = n_x$$
$$\frac{\partial f_n}{\partial f_y} = n_y$$
$$\frac{\partial f_n}{\partial f_z} = n_z$$

对轮子位置的导数（通过法向量）：
理论上，由于法向量依赖于轮子位置，需要对多项式系数求导：
$$\frac{\partial f_n}{\partial c_{lx}} = f_x \frac{\partial n_x}{\partial x_l} \frac{\partial x_l}{\partial c_{lx}} + f_y \frac{\partial n_y}{\partial x_l} \frac{\partial x_l}{\partial c_{lx}} + f_z \frac{\partial n_z}{\partial x_l} \frac{\partial x_l}{\partial c_{lx}}$$

**代码实现说明**：在 `NormalForceUpperBoundConstraint::eval_jac_g` 中，对轮子位置的梯度被简化为0（第143-145行）。这是因为：
1. 法向量对位置的依赖涉及地形梯度的二阶导数，计算复杂
2. 当地形变化不大时，这个近似是合理的
3. 如果需要更精确的梯度，可以通过数值微分或解析计算地形法向量的梯度

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.1.6 轮子离地距离约束（WheelAirborneConstraint）

**定义**：约束轮子与地形的距离在合理范围内（允许下陷和离地）。

**约束数量**：$2 \times n_{\text{sample}} = 50$ 个不等式约束

**数学表达式**：

左轮：
$$\text{MIN\_GROUND\_DISTANCE} \leq z_l - h(x_l, y_l) \leq \text{MAX\_AIRBORNE\_DISTANCE}$$

右轮：
$$\text{MIN\_GROUND\_DISTANCE} \leq z_r - h(x_r, y_r) \leq \text{MAX\_AIRBORNE\_DISTANCE}$$

**参数**：
- $\text{MIN\_GROUND\_DISTANCE} = -0.05$ m（允许5cm下陷）
- $\text{MAX\_AIRBORNE\_DISTANCE} = 0.15$ m（允许15cm离地）

**约束值**：
$$g_i = z_w - h(x_w, y_w)$$

**约束边界**：
- 下界：$g_{\text{lb}} = \text{MIN\_GROUND\_DISTANCE}$
- 上界：$g_{\text{ub}} = \text{MAX\_AIRBORNE\_DISTANCE}$

**雅可比矩阵**：

约束值 $g = z_w - h(x_w, y_w)$，对多项式系数的导数：

$$\frac{\partial g}{\partial c_{lz}} = \Phi_{lz}(s)$$
$$\frac{\partial g}{\partial c_{lx}} = -\frac{\partial h}{\partial x}(x_l, y_l) \cdot \Phi_{lx}(s)$$
$$\frac{\partial g}{\partial c_{ly}} = -\frac{\partial h}{\partial y}(x_l, y_l) \cdot \Phi_{ly}(s)$$

**代码实现说明**：在 `WheelAirborneConstraint::eval_jac_g` 中，使用 `terrain_->gradient()` 计算地形梯度 $\frac{\partial h}{\partial x}$ 和 $\frac{\partial h}{\partial y}$（第120-121行），然后乘以多项式基矩阵得到完整的雅可比。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

### 3.2 软约束（目标函数项）

软约束通过 `ObjectiveBase` 接口实现，作为目标函数的惩罚项。

#### 3.2.1 平滑项（SmoothingObjective）

**定义**：惩罚轨迹的二阶差分，使轨迹更平滑。

**数学表达式**：

$$f_{\text{smooth}} = \frac{1}{n_{\text{sample}}} \sum_{\text{traj}} \sum_{i=1}^{n-2} (v_{i+2} - 2v_{i+1} + v_i)^2$$

其中 $v$ 可以是 $x_l, y_l, z_l, x_r, y_r, z_r, x_b, y_b, z_b, \psi$。

**梯度**：

$$\frac{\partial f_{\text{smooth}}}{\partial c_k} = \frac{2}{n_{\text{sample}}} \sum_{i=1}^{n-2} (v_{i+2} - 2v_{i+1} + v_i) \cdot (\Phi_{i+2,k} - 2\Phi_{i+1,k} + \Phi_{i,k})$$

**海塞矩阵**：

$$\frac{\partial^2 f_{\text{smooth}}}{\partial c_k \partial c_j} = \frac{2}{n_{\text{sample}}} \sum_{i=1}^{n-2} (\Phi_{i+2,k} - 2\Phi_{i+1,k} + \Phi_{i,k}) \cdot (\Phi_{i+2,j} - 2\Phi_{i+1,j} + \Phi_{i,j})$$

#### 3.2.2 正则化项（RegularizationObjective）

**定义**：惩罚多项式系数的大小，防止过拟合。

**数学表达式**：

$$f_{\text{reg}} = \sum_{k} c_k^2$$

**梯度**：

$$\frac{\partial f_{\text{reg}}}{\partial c_k} = 2c_k$$

**海塞矩阵**：

$$\frac{\partial^2 f_{\text{reg}}}{\partial c_k \partial c_j} = 2\delta_{kj}$$

#### 3.2.3 轮子-中心距离标称项（WheelBaseDistanceNominalObjective）

**定义**：惩罚偏离标称距离。

**数学表达式**：

$$f_{\text{base\_nom}} = \frac{1}{n_{\text{sample}}} \sum_i [(d_L^2 - d_{\text{nom}}^2)^2 + (d_R^2 - d_{\text{nom}}^2)^2]$$

其中 $d_{\text{nom}} = 0.8$ m。

**梯度**：

$$\frac{\partial f_{\text{base\_nom}}}{\partial c_{bx}} = \frac{4}{n_{\text{sample}}} \sum_i (d_L^2 - d_{\text{nom}}^2) \cdot (x_b - x_l) \cdot \Phi_{bx}(s_i)$$

类似地，对所有相关系数都有对应的梯度项。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.4 左右轮距离标称项（WheelsDistanceNominalObjective）

**定义**：惩罚偏离标称左右轮距离。

**数学表达式**：

$$f_{\text{wheels\_nom}} = \frac{1}{n_{\text{sample}}} \sum_i (d_{LR}^2 - d_{\text{nom}}^2)^2$$

其中 $d_{\text{nom}} = 0.5$ m。

**梯度**：

$$\frac{\partial f_{\text{wheels\_nom}}}{\partial c_{lx}} = \frac{4}{n_{\text{sample}}} \sum_i (d_{LR}^2 - d_{\text{nom}}^2) \cdot (-(x_r - x_l)) \cdot \Phi_{lx}(s_i)$$

类似地，对所有相关系数都有对应的梯度项。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.5 轮子-地形距离项（WheelTerrainDistanceObjective）

**定义**：惩罚轮子与地形的距离。

**数学表达式**：

$$f_{\text{terrain}} = \frac{w}{n_{\text{sample}}} \sum_i [(z_l - h(x_l, y_l))^2 + (z_r - h(x_r, y_r))^2]$$

**权重**：$w = 20.0$

**梯度**：

$$\frac{\partial f_{\text{terrain}}}{\partial c_{lz}} = \frac{2w}{n_{\text{sample}}} \sum_i (z_l - h(x_l, y_l)) \cdot \Phi_{lz}(s_i)$$
$$\frac{\partial f_{\text{terrain}}}{\partial c_{lx}} = \frac{2w}{n_{\text{sample}}} \sum_i (z_l - h(x_l, y_l)) \cdot \left(-\frac{\partial h}{\partial x}(x_l, y_l)\right) \cdot \Phi_{lx}(s_i)$$
$$\frac{\partial f_{\text{terrain}}}{\partial c_{ly}} = \frac{2w}{n_{\text{sample}}} \sum_i (z_l - h(x_l, y_l)) \cdot \left(-\frac{\partial h}{\partial y}(x_l, y_l)\right) \cdot \Phi_{ly}(s_i)$$

右轮类似。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.6 轮子离地惩罚项（WheelAirbornePenaltyObjective）

**定义**：强烈惩罚轮子离地。

**数学表达式**：

$$f_{\text{airborne}} = \frac{w}{n_{\text{sample}}} \sum_i \max(0, z_w - h(x_w, y_w) - \text{threshold})^2$$

**参数**：
- $w = 1000.0$
- $\text{threshold} = 0.05$ m

**梯度**：

当 $z_w - h(x_w, y_w) > \text{threshold}$ 时：
$$\frac{\partial f_{\text{airborne}}}{\partial c_{lz}} = \frac{2w}{n_{\text{sample}}} \sum_i (z_l - h(x_l, y_l) - \text{threshold}) \cdot \Phi_{lz}(s_i)$$
$$\frac{\partial f_{\text{airborne}}}{\partial c_{lx}} = \frac{2w}{n_{\text{sample}}} \sum_i (z_l - h(x_l, y_l) - \text{threshold}) \cdot \left(-\frac{\partial h}{\partial x}(x_l, y_l)\right) \cdot \Phi_{lx}(s_i)$$

否则梯度为0。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.7 轨迹长度项（PathLengthObjective）

**定义**：惩罚轨迹长度，使轨迹尽可能短。

**数学表达式**：

$$f_{\text{length}} = w \sum_{i=1}^{n-1} \sqrt{(x_{i+1} - x_i)^2 + (y_{i+1} - y_i)^2}$$

**权重**：$w = 10.0$

**梯度**：

$$\frac{\partial f_{\text{length}}}{\partial c_{bx}} = w \sum_{i=1}^{n-1} \frac{(x_{i+1} - x_i)}{\sqrt{(x_{i+1} - x_i)^2 + (y_{i+1} - y_i)^2 + \epsilon}} \cdot (\Phi_{i+1,bx} - \Phi_{i,bx})$$

其中 $\epsilon$ 是防止除零的小常数。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.8 偏航-路径对齐项（YawPathAlignmentObjective）

**定义**：惩罚偏航角与路径切线方向不对齐。

**数学表达式**：

$$f_{\text{yaw}} = \frac{w}{n_{\text{sample}}} \sum_{i=1}^{n-1} (1 - \cos(\text{alignment}_i))^2$$

其中：
- $\mathbf{t}_i = \frac{[\Delta x, \Delta y]^T}{||[\Delta x, \Delta y]||}$（路径切线）
- $\mathbf{f}_i = [\cos\psi_i, \sin\psi_i]^T$（偏航方向）
- $\cos(\text{alignment}_i) = \mathbf{t}_i \cdot \mathbf{f}_i$

**权重**：$w = 50.0$

**梯度**：

$$\frac{\partial f_{\text{yaw}}}{\partial c_{\psi}} = \frac{2w}{n_{\text{sample}}} \sum_{i=1}^{n-1} (1 - \cos(\text{alignment}_i)) \cdot (t_{x,i} \sin\psi_i - t_{y,i} \cos\psi_i) \cdot \Phi_{\psi}(s_i)$$

对 $c_{bx}$ 和 $c_{by}$ 的梯度涉及路径切线对位置的依赖，计算较复杂。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.9 接触刚度项（ContactStiffnessObjective）

**定义**：确保接触力与穿透量符合弹簧模型。

**数学表达式**：

$$f_{\text{stiff}} = \frac{w}{n_{\text{sample}}} \sum_i [(f_{Ln} - K \cdot \text{pen}_l)^2 + (f_{Rn} - K \cdot \text{pen}_r)^2]$$

其中：
- $K = 500.0$ N/m（接触刚度）
- $\text{pen} = \max(0, h(x, y) - z)$（穿透量）

**权重**：$w = 100.0$

**梯度**：

对法向力的梯度：
$$\frac{\partial f_{\text{stiff}}}{\partial f_{Lx}} = \frac{2w}{n_{\text{sample}}} \sum_i (f_{Ln} - K \cdot \text{pen}_l) \cdot n_{Lx}$$

对轮子位置的梯度（通过穿透量）：
$$\frac{\partial f_{\text{stiff}}}{\partial c_{lz}} = \frac{2wK}{n_{\text{sample}}} \sum_i (f_{Ln} - K \cdot \text{pen}_l) \cdot \max(0, \text{sign}(h - z_l)) \cdot \Phi_{lz}(s_i)$$

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.10 力平衡项（ForceBalanceObjective）

**定义**：惩罚合力不为零（准静态平衡）。

**数学表达式**：

$$f_{\text{force}} = \frac{w}{n_{\text{sample}}} \sum_i ||\mathbf{f}_L + \mathbf{f}_R + [0, 0, -mg]^T||^2$$

**权重**：$w = 200.0$

**梯度**：

$$\frac{\partial f_{\text{force}}}{\partial f_{Lx}} = \frac{2w}{n_{\text{sample}}} (f_{Lx} + f_{Rx})$$
$$\frac{\partial f_{\text{force}}}{\partial f_{Lz}} = \frac{2w}{n_{\text{sample}}} (f_{Lz} + f_{Rz} - mg)$$

**海塞矩阵**：

$$\frac{\partial^2 f_{\text{force}}}{\partial f_{Lx}^2} = \frac{2w}{n_{\text{sample}}}$$

#### 3.2.11 力矩平衡项（TorqueBalanceObjective）

**定义**：惩罚合矩不为零。

**数学表达式**：

$$f_{\text{torque}} = \frac{w}{n_{\text{sample}}} \sum_i ||\boldsymbol{\tau}_L + \boldsymbol{\tau}_R||^2$$

其中：
- $\boldsymbol{\tau}_L = \mathbf{r}_L \times \mathbf{f}_L$
- $\boldsymbol{\tau}_R = \mathbf{r}_R \times \mathbf{f}_R$
- $\mathbf{r}_L = [x_l - x_b, y_l - y_b, z_l - z_b]^T$
- $\mathbf{r}_R = [x_r - x_b, y_r - y_b, z_r - z_b]^T$

**权重**：$w = 200.0$

**梯度**：

对接触力的梯度：
$$\frac{\partial f_{\text{torque}}}{\partial f_{Lx}} = \frac{2w}{n_{\text{sample}}} \sum_i [\tau_{Lx} \cdot (r_{Ly} \cdot 0 - r_{Lz} \cdot 0) + \tau_{Ly} \cdot (r_{Lz} \cdot 1 - r_{Lx} \cdot 0) + \tau_{Lz} \cdot (r_{Lx} \cdot 0 - r_{Ly} \cdot 1)]$$

实际上，由于 $\boldsymbol{\tau}_L = \mathbf{r}_L \times \mathbf{f}_L$，对 $\mathbf{f}_L$ 的梯度涉及叉积的导数，计算较复杂。

对位置的梯度（通过 $\mathbf{r}_L$ 和 $\mathbf{r}_R$）：
$$\frac{\partial f_{\text{torque}}}{\partial c_{lx}} = \frac{2w}{n_{\text{sample}}} \sum_i \boldsymbol{\tau}_{\text{total}} \cdot \frac{\partial \boldsymbol{\tau}_L}{\partial x_l} \cdot \Phi_{lx}(s_i)$$

其中 $\frac{\partial \boldsymbol{\tau}_L}{\partial x_l}$ 涉及叉积对位置的导数。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

#### 3.2.12 摩擦锥违反项（FrictionConeObjective）

**定义**：惩罚摩擦锥违反。

**数学表达式**：

$$f_{\text{fric}} = \frac{w}{n_{\text{sample}}} \sum_i [\max(0, ||\mathbf{f}_t||^2 - \mu^2 f_n^2)^2 + \max(0, -f_n)^2]$$

其中：
- $\mathbf{f}_t = \mathbf{f} - f_n \mathbf{n}$（切向力）
- $f_n = \mathbf{f} \cdot \mathbf{n}$（法向力）
- $\mu = 0.6$（摩擦系数）

**权重**：$w = 500.0$

**梯度**：

对接触力的梯度（当违反摩擦锥时）：
$$\frac{\partial f_{\text{fric}}}{\partial f_{Lx}} = \frac{2w}{n_{\text{sample}}} \sum_i [2(f_{tx} - \mu^2 f_n n_x) \cdot (f_{tx} - \mu^2 f_n n_x) + 2\max(0, -f_n) \cdot (-n_x)]$$

其中 $f_{tx} = f_x - f_n n_x$。

对轮子位置的梯度（通过法向量）：
由于法向量依赖于位置，梯度计算涉及地形梯度的二阶导数，计算复杂。代码中可能使用数值微分或简化处理。

**海塞矩阵**：代码中未实现解析海塞矩阵，使用 L-BFGS 近似（见第5节说明）。

## 4. 约束统计

### 4.1 硬约束统计

- **等式约束**：20个（边界约束）
- **不等式约束**：150个
  - 轮子-中心距离：50个
  - 左右轮距离：25个
  - 防交叉：25个
  - 法向力上限：50个
  - 轮子离地距离：50个（实际使用中可能被禁用）

**总约束数**：170个

### 4.2 软约束统计

**目标函数项数**：12个

### 4.3 变量统计

- **多项式系数**：50个
- **接触力**：150个
- **总变量数**：200个

## 5. 海塞矩阵实现说明

### 5.1 海塞矩阵近似方法

**重要说明**：代码中**未实现解析海塞矩阵**，而是使用 **L-BFGS（Limited-memory BFGS）近似方法**。

在 `main.cpp` 中（第473行）设置了：
```cpp
app->Options()->SetStringValue("hessian_approximation", "limited-memory");
```

在 `TrajOptNLP::eval_h` 中（第308-315行）直接返回 `true`，不提供海塞矩阵值：
```cpp
bool TrajOptNLP::eval_h(...) {
    // 当前版本不使用海塞矩阵
    return true;
}
```

### 5.2 L-BFGS 方法说明

L-BFGS 是一种准牛顿方法，通过维护一个有限内存的近似海塞矩阵来加速优化：
- **优点**：
  - 不需要计算和存储完整的海塞矩阵（对于200个变量，完整海塞矩阵有200×200=40000个元素）
  - 内存占用小（只存储最近几次迭代的梯度信息）
  - 对于大规模问题计算效率高
- **缺点**：
  - 收敛速度可能略慢于使用精确海塞矩阵的方法
  - 对于高度非线性的问题，近似可能不够精确

### 5.3 为什么使用 L-BFGS

1. **问题规模**：200个变量，170个约束，完整海塞矩阵计算和存储成本高
2. **实现复杂度**：许多约束和目标函数的海塞矩阵涉及地形梯度的二阶导数，计算复杂
3. **数值稳定性**：L-BFGS 方法在数值上更稳定，特别是对于约束优化问题

### 5.4 如果未来需要精确海塞矩阵

如果需要实现解析海塞矩阵以提高收敛速度，需要：
1. 实现 `TrajOptNLP::eval_h` 方法
2. 为每个约束和目标函数添加海塞矩阵计算方法
3. 特别注意地形相关项的二阶导数计算（RBF 地形的二阶导数）

## 6. 总结

本项目实现了一个完整的轨迹优化系统，包含：

1. **模块化设计**：约束和目标函数都是独立的类，易于扩展和维护
2. **完整的物理建模**：考虑力平衡、力矩平衡、摩擦锥等物理约束
3. **灵活的地形表示**：使用RBF表示复杂地形
4. **高效的数值优化**：使用Ipopt求解器，支持稀疏矩阵和L-BFGS海塞近似
5. **实用的简化处理**：
   - 法向力上限约束中对轮子位置的梯度使用简化（设为0）
   - 海塞矩阵使用L-BFGS近似而非解析计算

系统已通过代码清理，删除了未使用的约束文件（`force_balance_constraint`, `friction_cone_constraint`, `torque_balance_constraint`），这些功能通过对应的目标函数实现，提供了更好的数值稳定性。

