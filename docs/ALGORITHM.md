# 1. Problem Formulation

本文只描述当前源码能够确认的算法与数学模型。除明确标为 **Design / Intended formulation** 的内容外，公式均对应当前实现；无法从代码确认的物理含义标为 **UNKNOWN**。

给定由高斯 RBF 表示的连续 2.5D 地形

\[
z=f(x,y),
\]

以及二维起终点、机器人几何参数和摩擦参数，工程按两级方法求解：先在平面状态空间中生成粗路径，再以 Chebyshev 系数参数化基座和双轮轨迹并求解非线性规划。通用求解器还显式优化左右接触力。

期望输出包括：

- coarse path：\(p_0,\ldots,p_M\)，\(p_i=(x_i,y_i)\)；
- base trajectory：\(q_b(s)=[x_b,y_b,z_b,\psi]^T\)；
- left/right wheel trajectory：\(q_l(s),q_r(s)\)；
- 在 `trajopt_cpp` 中的离散接触力 \(f_l(s_i),f_r(s_i)\)。

**Design / Intended formulation.** 一条完整路线应经过地形分区、区域感知粗规划、路径分段、分别求解 roll/lift 段并连续拼接，同时满足几何、地形接触、摩擦与平衡条件。

**Current implementation.** 这些步骤由多个独立可执行程序组成；roll/lift 求解器各自只读取第一段目标标签，未实现全段遍历与拼接。几何条件多为硬约束；通用求解器的接触刚度、力平衡、力矩平衡和摩擦锥是软目标，不是严格可行性条件。

# 2. RBF Terrain Representation

## 2.1 Current C++ model

`trajopt_cpp/src/rbf_terrain.cpp` 和 CUDA kernel 使用相同核函数。中心 \(c_i=(c_{xi},c_{yi})\)、权重 \(w_i\)、带宽 \(\sigma>0\) 定义

\[
\phi_i(x,y)=\exp\!\left[-\frac{(x-c_{xi})^2+(y-c_{yi})^2}{2\sigma^2}\right],
\qquad
f(x,y)=\sum_{i=1}^{N}w_i\phi_i(x,y).
\]

解析梯度为

\[
f_x=\frac{\partial f}{\partial x}
=\sum_i w_i\phi_i\left(-\frac{x-c_{xi}}{\sigma^2}\right),
\qquad
f_y=\frac{\partial f}{\partial y}
=\sum_i w_i\phi_i\left(-\frac{y-c_{yi}}{\sigma^2}\right).
\]

未归一化向上法向和单位法向为

\[
n_u=[-f_x,-f_y,1]^T,
\qquad
\hat n=\frac{n_u}{\sqrt{f_x^2+f_y^2+1}}.
\]

代码在法向范数判断中使用 \(10^{-8}\)；理论上该范数至少为 1，fallback `[0,0,1]` 只是防御分支。

`JSONReader` 只读取：

- `cur_rbf_grid`：\(N\times2\) centers；
- `rbf_weight`：长度 \(N\) 的 weights。

二者长度必须一致。C++ reader 不从 JSON 读取 \(\sigma\)，而是在 `trajopt_cpp/src/json_reader.cpp` 和 `terrain_friction_segment/src/json_reader.cpp` 中固定赋值 `0.14`。`solve_casadi.py` 与 `trajopt_sdf.py` 的运行入口也设置 `0.14`；`geometry/terrain.py` 是单中心演示模型，示例 `ell=0.35`，不属于 C++ 数据链。`RBF_LIO.pdf` 报告的论文实验值为 `0.04 m`，不能据此推断本仓库输入应自动采用该值。

本仓库没有实现 RBF-LIO 的点云处理、中心在线选择、权重递归更新或 LIO 状态估计。

# 3. Friction-Aware Terrain Segmentation

对规则网格点 \((x_j,y_j)\)，CUDA kernel 计算 \(f_x,f_y,\hat n\)。令坡角

\[
\alpha=\arctan\sqrt{f_x^2+f_y^2},
\qquad
\hat n_z=\frac{1}{\sqrt{1+f_x^2+f_y^2}}.
\]

源码把竖直支撑力落在库仑摩擦锥内写为

\[
\hat n_z\ge \frac{1}{\sqrt{1+\mu^2}},
\]

等价于

\[
\sqrt{f_x^2+f_y^2}\le\mu
\quad\Longleftrightarrow\quad
\alpha\le\arctan\mu.
\]

分类为

\[
L(x,y)=
\begin{cases}
1,&\hat n_z\ge 1/\sqrt{1+\mu^2}\quad(rollable),\\
0,&\text{otherwise}\quad(lift).
\end{cases}
\]

因此更大的 \(\mu\) 降低 \(n_z\) 阈值，使更多陡坡被标为 rollable；更小的 \(\mu\) 更严格。CLI 默认 `mu=0.3`，并将输入截断到 `[0.1,1.0]`。当前网格为 `512 x 512`，边界是 RBF centers 包围盒向外扩 `0.3 m`。每个 CUDA thread 遍历全部 centers，复杂度约为 \(O(n_xn_yN)\)。高度网格随后在 CPU 上重新计算，用于文本与图片输出。

这里没有概率分类、粗糙度模型或空间变化摩擦系数；标签是二值坡度可行性结果。

# 4. Region-Aware RRT*

实现位于 `terrain_friction_segment/src/rrt_star_region.cpp`。

## 4.1 State, sampling, and steering

- 状态空间：分区网格边界内的二维矩形 \(\mathcal X=[x_{min},x_{max}]\times[y_{min},y_{max}]\)。
- 采样：固定随机种子 `0`；概率 `0.1` 直接采 goal，否则在矩形中均匀采样。
- nearest：先用 `64 x 64` bucket 空间索引，在局部桶中找距离不超过 `(3*step_length)^2` 的最近节点；若没有候选，再全树扫描。
- steering：若样本距离不超过 `step_length=0.15`，直接到样本，否则沿直线前进 `0.15`。
- near nodes：固定 `rewire_radius=0.5`，不是随节点数缩小的半径。

## 4.2 Validity and edge cost

当前唯一显式状态有效性检查是新节点落在矩形边界内。源码没有障碍物碰撞检测，也没有禁止 label `0`；两类标签都可通过。

边 \(a\to b\) 被均匀分成 \(M_e=20\) 个小段，在每段中点查询最近网格标签。设 \(\Delta\ell=\|b-a\|/M_e\)，实际代价为

\[
c(a,b)=\sum_{m=0}^{M_e-1}\Delta\ell\,
\rho\!\left(L\left(a+\frac{m+1/2}{M_e}(b-a)\right)\right),
\quad
\rho(1)=1,\quad \rho(0)=k_{lift},
\]

其中默认 \(k_{lift}=5\)。累计节点代价是父节点累计代价加边代价。

## 4.3 Parent selection, rewiring, and goal

新节点先以 nearest 为父，再遍历固定半径邻居，选择使

\[
J_{root}(i)+c(i,new)
\]

最小者。插入后，如果经新节点到邻居成本更小，则更新该邻居的 parent 和 cost。新节点距 goal 不超过 `goal_radius=0.2` 时，额外计算连接 goal 的边；更低成本的 goal 路径通过 parent 指针回溯并反转。

**Design / Intended formulation.** 名称和结构意图是带区域代价的 RRT*。

**Current implementation.** 与经典 RRT* 的重要差异是固定重连半径、固定迭代数 `2000`，并且 rewire 后不递归更新该节点后代的累计 cost。因此不能仅依据名称声称当前实现具有标准 RRT* 的渐近最优保证。

粗路径输出每行 `(x,y)`。`save_waypoints_segmented` 对每个路径节点查询最近网格标签并输出 `(x,y,label)`；它不写显式 segment ID。后续 reader 把连续相同标签视为区段。

# 5. Trajectory Parameterization

## 5.1 Chebyshev basis

令归一化路径参数 \(s\in[0,1]\)，代码映射

\[
\xi=2s-1\in[-1,1].
\]

任一标量轨迹写为

\[
q(s)=\sum_{k=0}^{K-1}c_kT_k(2s-1),
\]

其中 \(T_0=1,T_1=\xi,T_{k}=2\xi T_{k-1}-T_{k-2}\)。利用 \(T_k'(\xi)=kU_{k-1}(\xi)\)，源码的导数关系是

\[
\frac{dq}{ds}=2\sum_kc_kT_k'(\xi),
\qquad
\frac{d^2q}{ds^2}=4\sum_kc_kT_k''(\xi).
\]

端点 \(\xi=\pm1\) 使用解析极限，避免含 \(1-\xi^2\) 的表达式除零。

优化 coefficients 而不是逐 waypoint，意味着轨迹天然处于低维全局函数空间，任一系数同时影响多个采样点；平滑性可由系数或采样差分控制。但这也会使局部约束相互耦合，并增加初值拟合敏感性。当前 `s` 是无时间量纲的路径参数，不是秒；导数也不是物理速度/加速度，除非另行给出时间参数化。

## 5.2 Piecewise wheel height in lift mode

lift 模式把 \([0,1]\) 均分成 \(N_z\) 段。对全局 \(s\)，

\[
j=\min(\lfloor sN_z\rfloor,N_z-1),\qquad
u=\operatorname{clip}(sN_z-j,0,1),
\]

并用第 \(j\) 段系数计算 \(z_w(s)=\sum_kc_{wz,jk}T_k(2u-1)\)。相邻段只加 C0 近似连续约束，没有 C1/C2 连续约束。

# 6. General Trajectory Optimization Problem

`trajopt_cpp` 实现标准形式

\[
\min_x J(x)=\sum_a J_a(x),
\quad\text{s.t.}\quad
g_L\le g(x)\le g_U,\qquad x_L\le x\le x_U.
\]

当前所有原始变量边界都设为近似无界 `[-1e20,1e20]`，实际限制主要来自 `g`。Ipopt callback 对应关系为：

- `eval_f`：累加所有启用目标的标量值；
- `eval_grad_f`：先清零，再累加各目标解析梯度；
- `eval_g`：按注册顺序拼接硬约束值；
- `eval_jac_g`：第一次返回稀疏行列结构，随后返回对应 Jacobian 值；
- `eval_h`：不提供 Hessian，入口选择 `hessian_approximation=limited-memory`。

# 7. Decision Variables

| 优化器 | 当前规模 | 实际变量顺序 |
|---|---:|---|
| `trajopt_cpp` | `K=5`, `Ns=25`；`10K+6Ns=200` | `cbx,cby,cbz,cpsi,clx,cly,clz,crx,cry,crz`，每块 K；随后 `fLx,fLy,fLz,fRx,fRy,fRz`，每块 Ns。 |
| `seg_trajopt_cpp` | `K=5`, `Ns=15`；`8K=40` | `cbx,cby,cbz,cpsi,clx,cly,crx,cry`。`zl=f(xl,yl)`、`zr=f(xr,yr)`，不占变量。 |
| `lift_leg_trajopt_cpp` | `K=8`, `Ns=21`；`8K+2KNz` | `cbx,cby,cbz,cpsi,clx,cly,clz[0:Nz*K],crx,cry,crz[0:Nz*K]`。没有显式接触力。 |

这里 `cb*` 表示 base 系数，`cl*`/`cr*` 表示左右轮系数。通用版本的 \(f_w(s_i)=[f_{wx,i},f_{wy,i},f_{wz,i}]\) 是每个采样点独立变量。

# 8. Objective Functions

记 \(d_{bw}^2=\|p_b-p_w\|^2\)、\(d_{lr}^2=\|p_l-p_r\|^2\)、间隙 \(c_w=z_w-f(x_w,y_w)\)、\([a]_+=\max(0,a)\)。下表的“权重”是传入 objective 后的实际乘数。

## 8.1 `trajopt_cpp`

当前 `main.cpp` 注册并启用了该目录下全部 13 个 objective：

| 目标 | Implemented form、输入和作用 | 当前权重 |
|---|---|---:|
| Smoothing | \(\sum_{q\in10条轨迹}\sum_i(q_{i+2}-2q_{i+1}+q_i)^2\)。惩罚采样值二阶差分；不作用于力。 | `1/Ns = 0.04` |
| Regularization | \(\sum_{10条轨迹}\|c_q\|_2^2\)。抑制系数幅值；不正则化力。 | `1.0` |
| WheelBase nominal | \(\sum_i[(d_{bl,i}^2-0.8^2)^2+(d_{br,i}^2-0.8^2)^2]\)。趋向名义轮-基座距离。 | `1/Ns = 0.04` |
| Wheels nominal | \(\sum_i(d_{lr,i}^2-0.5^2)^2\)。趋向名义轮距。 | `1/Ns = 0.04` |
| WheelTerrainDistance | \(\sum_i(c_{l,i}^2+c_{r,i}^2)\)。同时惩罚离地和穿透。 | `20/Ns = 0.8` |
| Airborne penalty | \(\sum_i([c_{l,i}-0.05]_+^2+[c_{r,i}-0.05]_+^2)\)。减少超过 5 cm 的离地。 | `1000/Ns = 40` |
| Path length | 基座三维弧长 \(\int_0^1\|p_b'(s)\|ds\)，用均匀采样梯形法近似。 | `10` |
| Forward progress | 令 \(d\) 为归一化起点到终点方向，\(r_i=(p_{b,i+1}^{xy}-p_{b,i}^{xy})\cdot d\)，惩罚 \(\sum_i P_\epsilon(-r_i)^2\)，其中 \(P_\epsilon(a)=\tfrac12(a+\sqrt{a^2+10^{-6}})\)。减少回退。 | `300/Ns = 12` |
| Yaw/path alignment | 相邻基座点产生归一化二维切线 \(t_i\)，使用中点 yaw \(\bar\psi_i\)，惩罚 \(\sum_i[1-t_i\cdot(\cos\bar\psi_i,\sin\bar\psi_i)]^2\)。 | 默认 `50/Ns = 2` |
| Contact stiffness | \(\sum_{w,i}[f_{n,w,i}-500[c_{penetration,w,i}]_+]^2\)，其中 penetration \(=f(x_w,y_w)-z_w\)，\(f_n=f\cdot\hat n\)。软拟合弹簧接触。 | `100/Ns = 4` |
| Force balance | \(\sum_i\|f_{l,i}+f_{r,i}+[0,0,-mg]^T\|^2\)，当前 `m=1 kg`,`g=9.81`。准静态合力软平衡。 | `200/Ns = 8` |
| Torque balance | \(\sum_i\|(p_l-p_b)\times f_l+(p_r-p_b)\times f_r\|^2\)。关于 base 点的合矩软平衡。 | `200/Ns = 8` |
| Friction cone | 对每轮惩罚 \([\|f_t\|^2-\mu^2f_n^2]_+^2+[-f_n]_+^2\)，\(f_t=f-f_n\hat n\)，`mu=0.6`。 | `500/Ns = 20` |

所有项的 `enabled_` 在构造时为 true，入口没有关闭任何一项。源码没有“已实现但默认未启用”的通用 objective 文件。

## 8.2 Segment objectives

| 优化器 | 目标 | Implemented form | 当前权重 |
|---|---|---|---:|
| seg | InitialGuessProximity | \(\|x-x_{ref}\|_2^2\)，使细化结果不偏离 waypoint 拟合。 | `0.01` |
| seg | BaseHeightStability | \(\sum_i[z_b-f(x_b,y_b)-0.8]^2\)。 | `2/15` |
| seg | YawPathAlignment | 与通用版相同。 | `10/15` |
| lift | InitialGuessProximity | \(\|x-x_{ref}\|_2^2\)。 | `0.02` |
| lift | BaseHeightStability | \(\sum_i[z_b-f(x_b,y_b)-0.8]^2\)。 | `2/21` |
| lift | YawPathAlignment | 与通用版相同。 | `8/21` |
| lift | SwingWheelLift | \(\sum_i[(c_l-c_l^{ref})^2+(c_r-c_r^{ref})^2]\)，摆动参考为 `0.24*sin^2(pi*s)`，支撑参考为 `0.005`。 | `6/21` |

**Current implementation.** `SwingWheelLiftObjective` 接收一个全段布尔量 `swing_left`，不是逐采样 contact mask。当左右轮都在初始化中出现摆动窗口时，该布尔量保持左轮为摆动轮，因此该目标的参考与多窗口交替 contact mask 并不完全相同。

# 9. Constraints

## 9.1 Shared geometric form

| 约束 | Current implemented form、上下界与物理含义 |
|---|---|
| Boundary | 通用版有 20 个严格等式：base 起终 `x,y,z=f(x,y)+0.8,yaw`；双轮起终平面位置由 \(u_{lat}=(-\sin\psi,\cos\psi)\) 和名义轮距构造，轮高等于目标轮位置处地形。seg 用相同几何残差但允许 base `xy/z/yaw` 容差 `1e-3`、wheel xy 容差 `0.05`。lift 再增加 wheel z 地形残差容差 `0.02` 及端点 z 范围 `[-2,3]`。 |
| Trajectory bounds | 所有采样处 base、左右轮的 `x` 在 `[xmin,xmax]`、`y` 在 `[ymin,ymax]`；不直接限制中间 z/yaw。 |
| Wheel-base distance | 每个采样、每轮约束 \(d_{bw,min}^2\le\|p_b-p_w\|^2\le d_{bw,max}^2\)。通用参数 `[0.4,1.2] m`；seg/lift `[0.2,1.2] m`。 |
| Wheel distance | \(d_{lr,min}^2\le\|p_l-p_r\|^2\le d_{lr,max}^2\)。通用 `[0.1,0.8] m`；seg/lift `[sqrt(0.1),1.0] m`。 |
| Anti-crossing | 仅通用版：\(((p_l^{xy}-p_b^{xy})-(p_r^{xy}-p_b^{xy}))\cdot u_{lat}\ge0\)，保持左右顺序。 |
| Coefficient norm | seg/lift：源码计算 \(g=\|x\|_2^2/n_{var}\le C^2\)，`C=10000`。名称虽为 L2 norm，实际约束的是系数 RMS 不超过 C。 |

seg 的 wheel z 由地形高度计算；lift 的距离约束使用显式分段 wheel z。

## 9.2 Terrain/contact constraints

| 约束 | Current implemented form、上下界与状态 |
|---|---|
| WheelAirborne | 通用版硬约束 \(-0.05\le z_w-f(x_w,y_w)\le0.15\) m。允许 5 cm 穿透和 15 cm 离地。 |
| NormalForceUpperBound | 通用版硬约束 \(f_w\cdot\hat n_w\le m(g+a_{max})=10.81\) N；没有法向力下界。 |
| RollableRegion | seg 在分区文件成功加载时启用。对 base 和两轮的双线性插值标签 \(S(x,y)\)，约束 `threshold-S <= buffer`，即 \(S\ge threshold-buffer\)。当前二者均为 `0.1`，实际为 \(S\ge0\)，对 `[0,1]` 标签不排除 lift 区。 |
| Terrain penetration | lift：\(f(x_w,y_w)-z_w\le0.01\) m，无下界，故抬高不受该约束限制。 |
| Penetration-based force | lift 定义近似力 \(F_w=5000P_\epsilon(f-z_w)\)，\(P_\epsilon(d)=\tfrac12(d+\sqrt{d^2+10^{-6}})\)。contact mask 为 1 时 `5 <= F <= 700 N`；为 0 时 `F <= 8 N` 且无有效下界。它是硬约束在“估算力”上的应用，但该力不是决策变量。 |
| Piecewise z continuity | lift 且 `Nz>1` 时：相邻段端点差 \(|z_{w,j}(1)-z_{w,j+1}(0)|\le0.005\) m；仅 C0。 |

# 10. Contact and Force Model

## 10.1 `trajopt_cpp`

| 机制 | 类型 | 实际含义 |
|---|---|---|
| wheel clearance bounds | hard constraint | 把轮-地高度差限制在 `[-0.05,0.15] m`。 |
| normal-force upper bound | hard constraint | 只限制 \(f_n\) 上界，没有硬性 \(f_n\ge0\)。 |
| contact stiffness | soft objective | 使 \(f_n\) 接近 \(K[f-z]_+\)，`K=500 N/m`；`[.]_+` 实际是非光滑分支，不是平滑函数。 |
| force balance | soft objective | 准静态合力残差平方。没有惯性项 \(ma\)。 |
| torque balance | soft objective | 关于 base 位置的接触合矩残差平方。 |
| friction cone | soft objective | 惩罚锥外力与拉力；不保证最终严格满足。 |

因此该模型是“离散接触力 + 准静态软平衡 + 弹簧接触软拟合”，不是完整刚体动力学或互补接触问题。

## 10.2 `lift_leg_trajopt_cpp`

lift 版没有 \(f_x,f_y,f_z\)。它以穿透 \(d=f-z_w\) 和刚度 `5000 N/m` 计算标量估算力，再依 contact mask 施加区间。摩擦只在**初值构造**中使用 `|grad f| <= 0.55`：支撑 waypoint 若不满足，则在半径 `0.25 m`、`7 x 7` 候选网格内选择坡度最小点。该摩擦条件没有作为 lift NLP 的持续约束。

# 11. Rollable Segment Optimization

输入为 `rbf.json`、`waypoints_segmented.txt`，以及可选成功加载的 `terrain_segmentation_data.txt`。`WaypointReader` 只返回第一段 label `1` 连续路径。

waypoint 参数按索引均匀取 \(s_i=i/(N-1)\)，不是按弧长。入口根据路径差分计算 yaw，以名义轮距 `0.5 m` 沿 \((-\sin\psi,\cos\psi)\) 放置左右轮；base 高度为 `terrain+0.8`，轮高由 terrain 隐式决定。8 组初值系数通过 QR/回代拟合 waypoint 值。

NLP 优化 40 个系数，启用 Boundary、TrajectoryBounds、CoefficientNorm、WheelBaseDistance、WheelsDistance；分区加载成功时再启用 RollableRegion。目标为 InitialGuessProximity、BaseHeightStability 和 YawPathAlignment。输出为带 `s,base,yaw,left,right` 的 initial/final trajectory；wheel z 在提取时重新由 RBF 计算。

**Design / Intended formulation.** rollable 段应留在 label `1` 区域。

**Current implementation.** 默认 rollable 约束化简为插值标签 \(S\ge0\)，因此不会从数学上排除 label `0`。

# 12. Lift-Leg Segment Optimization

入口优先读取第一段 label `0`，若不存在则回退到第一段 label `1`。原始路径按线性插值加密到 `max(12*Nraw,25)` 点，并使用均匀索引参数。

高度段数由

\[
N_z=\max\left(1,\left\lceil\frac{\max(0,h_{end}-h_{start})}{0.20}\right\rceil\right)
\]

确定，而不是从完整高度曲线聚类。每个上升 cycle 在路径进度附近生成两个错开的摆动窗口，左右轮交替；窗口参考高度在线性台阶高度上叠加无接触裕量 `0.03 m` 和 lift profile。局部 peak 为 `max(0.24, step_height+0.06)`。支撑点轮高初始化为 `terrain+0.005 m`，摆动点至少为 `terrain+0.03 m`。接触点还会执行前述摩擦安全平面位置投影。

初值经 QR 拟合：base/yaw/wheel xy 用全局 K=8 系数，wheel z 分段拟合。NLP 使用第 8、9 节所列四个目标和七类约束；没有 rollable-region 约束，也没有显式接触力或力/矩平衡。

与 rollable 优化器的根本区别是：rollable 版令 \(z_w=f(x_w,y_w)\)，而 lift 版把 \(z_w\) 作为分段轨迹变量，通过摆动参考、穿透上界、估算力区间和 C0 连接控制接触模式。

# 13. Complete Planning Pipeline

算法链可写为

\[
\{c_i,w_i,\sigma\}
\xrightarrow{RBF}
(f,\nabla f,\hat n)
\xrightarrow{\|\nabla f\|\le\mu}
L(x,y)
\xrightarrow{RRT^*\;region\;cost}
\{p_i,L_i\}
\xrightarrow{run\;extraction}
\mathcal P_{roll},\mathcal P_{lift}
\xrightarrow{NLP}
q_b,q_l,q_r[,f_l,f_r].
\]

其中 `[contact force]` 只属于独立通用 `trajopt_cpp`；当前 segment pipeline 的 roll 版没有力变量，lift 版只有穿透力近似。当前代码没有把通用 `trajopt_cpp` 自动接在 segment pipeline 末端。

# 14. Algorithm Pseudocode

## Algorithm 1: Terrain friction segmentation

```text
Input: centers c_i, weights w_i, fixed sigma=0.14, mu
Clamp mu to [0.1, 1.0]
Create 512x512 grid over center bounding box plus 0.3 m margin
nz_min <- 1 / sqrt(1 + mu^2)
parallel for each grid point (x,y) on CUDA:
    fx, fy <- sum Gaussian RBF analytic gradients
    n <- normalize([-fx, -fy, 1])
    label <- 1 if n.z >= nz_min else 0
Recompute f(x,y) on CPU for output
Write heights and labels
```

## Algorithm 2: Region-aware RRT*

```text
Input: label grid, start, goal, parameters
Initialize tree with start; RNG seed <- 0
for iter = 1..2000:
    sample goal with probability 0.1, else uniform point in bounds
    nearest <- bucket nearest, with global fallback
    new <- steer(nearest, sample, 0.15)
    reject new only if outside rectangular bounds
    near <- tree nodes within radius 0.5
    parent <- argmin_i cost_root[i] + region_edge_cost(i,new)
    insert new
    for i in near:
        if cost_root[new] + region_edge_cost(new,i) < cost_root[i]:
            parent[i] <- new; update cost_root[i]
            # current code does not propagate update to descendants
    if distance(new,goal) <= 0.2:
        evaluate direct goal edge and retain lower-cost goal path
Backtrack parent pointers from retained goal
Assign nearest-grid label to every returned path node
```

## Algorithm 3: Segment trajectory optimization

```text
Input: RBF, labeled waypoints, requested mode
Read first contiguous waypoint run matching mode
Parameterize points uniformly by waypoint index
Construct base/yaw/left/right geometric reference
if mode == roll:
    set wheel z implicitly to terrain height
    fit 8 global Chebyshev coefficient blocks
    add roll constraints/objectives
else if mode == lift:
    densify path; estimate Nz from positive net height gain
    build alternating swing windows and contact masks
    project stance references toward locally low-slope points
    construct stance/swing wheel-z references
    fit global xy/base/yaw and piecewise wheel-z coefficients
    add lift constraints/objectives
Call Ipopt through TNLP callbacks
Extract sampled initial and final trajectories
```

## Algorithm 4: Intended end-to-end planning pipeline

```text
[MANUAL] Run terrain_friction_segment on RBF JSON
[MANUAL] Run rrt_star_region on generated segmentation data
[MISSING ORCHESTRATOR] Split/iterate all contiguous labeled path segments
for every segment in route:
    [MANUAL, first matching segment only today]
    run seg_trajopt_cpp for label 1 or lift_leg_trajopt_cpp for label 0
[NOT IMPLEMENTED] Check adjacent-segment state/derivative/contact continuity
[NOT IMPLEMENTED] Stitch all segment trajectories into one route
[NOT IMPLEMENTED] Time-parameterize and export to a robot controller
```

# 15. Assumptions

- 地形是单值 2.5D 图面 \(z=f(x,y)\)；悬垂、多层表面和垂直壁面不能唯一表示。
- centers、weights 一一对应，坐标和高度按代码参数使用；JSON 未提供单位 schema，单位来源为 **UNKNOWN**，但常量和文档按米解释。
- 分区假设竖直重力支撑与局部库仑摩擦锥比较可用 \(\|\nabla f\|\le\mu\) 表示；没有轮滑速度或动力学判定。
- 通用力模型采用准静态软平衡；没有 base 线/角加速度惯性项。
- 每轮在采样点有一个局部地形法向和一个等效接触点；接触刚度是线性法向弹簧近似。
- 机构以 base 点、左右轮点和 yaw 表示；roll/lift 的轮宽、轮半径、关节运动学和完整刚体姿态未显式建模。
- `s` 是路径参数，采样间隔不是固定物理时间。

# 16. Numerical Considerations

- **Coefficient scaling.** C++ NLP 没有实现变量 scaling callback；不同单位的轨迹系数和牛顿力变量直接共存。Python `core/variable.py` 支持 scaling，但不作用于 C++。
- **Chebyshev conditioning.** 使用 Chebyshev 基比高阶幂基通常更稳定；当前 K 较小（5 或 8）。初始化以 QR/回代拟合，但 waypoint 的 `s` 按索引均匀，不按弧长。
- **Normalization.** RBF 法向有 `1e-8` 防护；yaw tangent 使用 `sqrt(dx^2+dy^2+1e-6)`；RRT 坐标到网格含小 epsilon。
- **RBF exponential.** 离中心很远时指数可下溢到 0；源码没有截断指数或邻域加速，每次评估遍历全部中心。
- **Ipopt.** 使用 MUMPS、adaptive barrier 和 limited-memory Hessian。目标权重跨度很大，且没有统一量纲归一化。
- **Jacobians.** 大部分约束/目标提供手写梯度/Jacobian。`NormalForceUpperBoundConstraint` 明确把法向对轮位置的导数置 0；`FrictionConeObjective` 也把该位置依赖置 0；`ContactStiffnessObjective` 的位置梯度只写 wheel-z 穿透项，未包含地形 xy 和法向变化。
- **Known derivative mismatch.** `WheelAirbornePenaltyObjective::eval_f` 使用 `clearance-threshold`，但 `eval_grad_f` 使用未减 threshold 的 clearance，二者在 `0 < clearance < 0.05` 等区间不一致。
- **Positive-part functions.** 通用 contact/friction/airborne 中名为 `smoothPos` 的函数实际是分段 `max(0,x)`；lift force 使用真正的 \(P_\epsilon\)。
- **Contact stiffness.** 较大刚度把很小的 penetration 放大为大力，并与 `MAX_PENETRATION`、力上下界共同造成敏感尺度。
- **Initial guess sensitivity.** 低维全局系数、非凸地形约束、离散 contact mask 和大权重罚项使 Ipopt 强依赖 waypoint 拟合。`finalize_solution` 会保存最终 iterate；输出文件存在不等同于求解成功，必须检查 solver status。

# 17. Current Limitations

1. 没有统一端到端 orchestrator、全 segment 遍历、段间连续性约束或 stitching。
2. roll 默认区域硬约束实际不排除 label `0`。
3. lift 的 \(N_z\) 只由正向净高度差决定；下降段或起终等高但中间有台阶时不会据完整地形自动增加高度段。
4. lift 的摆动目标使用单个全段 `swing_left`，而接触力区间使用逐采样交替 mask，两者不完全同构。
5. lift 摩擦安全性只用于初值 waypoint 投影，不是 NLP 约束。
6. 通用接触/平衡/摩擦主要是 soft penalty；不保证严格力学可行。
7. 通用 normal-force 硬约束没有非负下界，且若干力相关解析导数省略地形法向的位置导数。
8. airborne 目标值与解析梯度存在 threshold 使用不一致。
9. RRT* 没有障碍碰撞检查，label `0` 仅增加代价；rewire 不传播后代 cost。
10. `sigma`、起终点、机构参数、权重和大量数值阈值写在入口源码中，没有统一配置来源。
11. 轨迹没有时间参数化、速度/加速度物理界、执行控制或在线重规划接口。
12. 通用 `trajopt_cpp` 与 roll/lift pipeline 是并行实现，不是自动组合的一套完整动力学规划器。

# 18. Symbol Table

| Symbol | Meaning | Unit | Source / current parameter |
|---|---|---|---|
| \((x,y,z)\) | 世界坐标位置 | m（代码使用语境；JSON schema 未声明） | RBF/trajectory code |
| \(c_i\) | 第 i 个 RBF 平面中心 | m | `cur_rbf_grid` |
| \(w_i\) | 第 i 个 RBF 高度权重 | 与 z 相同 | `rbf_weight` |
| \(\sigma\) | Gaussian bandwidth | m | C++ `JSONReader`: `0.14` |
| \(f,f_x,f_y\) | 地形高度及平面梯度 | m；梯度 m/m | `RBFTerrain` |
| \(\hat n\) | 向上单位地形法向 | 1 | `RBFTerrain::normal` |
| \(\mu\) | 摩擦系数 | 1 | segmentation default `0.3`; full force objective `0.6`; lift init `0.55` |
| \(L\) | 区域标签 | 1 | `0=lift`, `1=rollable` |
| \(k_{lift}\) | RRT lift 区长度倍率 | 1 | `5.0` |
| \(s\) | 全局归一化路径参数 | 1 | `[0,1]` |
| \(u\) | lift 分段局部参数 | 1 | `[0,1]` |
| \(T_k,U_k\) | 第一/第二类 Chebyshev 多项式 | 1 | `ChebyshevBasis` |
| \(K\) | 每条轨迹系数数目 | 1 | full/seg `5`; lift `8` |
| \(N_s\) | NLP 采样数 | 1 | full `25`; seg `15`; lift `21` |
| \(N_z\) | lift wheel-z 分段数 | 1 | positive net rise / `0.20` |
| \(p_b,p_l,p_r\) | base、左轮、右轮位置 | m | trajectory variables |
| \(\psi\) | base yaw | rad | `cpsi` |
| \(d_{bw},d_{lr}\) | base-wheel、wheel-wheel 距离 | m | geometry constraints |
| \(c_w=z_w-f(x_w,y_w)\) | wheel terrain clearance | m | terrain objectives/constraints |
| \(f_l,f_r\) | 通用优化器接触力 | N | six discrete force blocks |
| \(f_n,f_t\) | 法向和切向接触力 | N | force objectives |
| \(m,g\) | 质量、重力加速度 | kg, m/s² | full main: `1.0`, `9.81` |
| \(K_c\) | 接触刚度 | N/m | full `500`; lift approximation `5000` |
| \([a]_+\) | positive part `max(0,a)` | follows a | full penalties |
| \(P_\epsilon(a)\) | smooth positive approximation | follows a | lift force, `epsilon=1e-6` |
