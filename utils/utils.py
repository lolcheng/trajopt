import json
import numpy as np
import matplotlib.pyplot as plt

def read_from_json(file_name):
    # 读取 JSON 文件
    with open(file_name, "r") as f:
        data = json.load(f)

    # ------ 读取 timestamp ------
    timestamp = data["timestamp"]
    print("timestamp:", timestamp)

    # ------ 读取 pose ------
    pos = data["pose"]["position"]
    ori = data["pose"]["orientation"]
    print("position:", pos)
    print("orientation:", ori)

    # ------ mean_weight ------
    mean_weight = data["mean_weight"]
    # print("mean weight:", mean_weight)

    # ------ cur_rbf_grid (N x 2 数组) ------
    cur_rbf_grid = data["cur_rbf_grid"]    # list of [x, y]
    # print("cur_rbf_grid:", cur_rbf_grid)

    # ------ rbf_weight (list<double>) ------
    rbf_weight = data["rbf_weight"]
    # print("rbf_weight:", rbf_weight)

    return timestamp, pos, ori, mean_weight, np.array(cur_rbf_grid), np.array(rbf_weight)

def plot_on_terrain(trajs, terrain=None, domain=None):
    if terrain is None:
        raise ValueError("terrain must be provided for plotting RBF terrain")
    (Xb, Yb, Zb, Psib), (Xl, Yl, Zl), (Xr, Yr, Zr) = trajs

    if domain is None:
        # centers 是 (N,2) 的数组
        xmin, ymin = terrain.centers.min(axis=0)
        xmax, ymax = terrain.centers.max(axis=0)
        margin = 0.2     # 适当留一点边界
        xmin -= margin
        xmax += margin
        ymin -= margin
        ymax += margin
    else:
        # 允许 domain=(xmin, xmax, ymin, ymax)
        if len(domain) == 2:
            xmin, xmax = domain
            ymin, ymax = domain
        elif len(domain) == 4:
            xmin, xmax, ymin, ymax = domain
        else:
            raise ValueError("domain must be None, (min,max) or (xmin,xmax,ymin,ymax)")

    xs = np.linspace(xmin, xmax, 161)
    ys = np.linspace(ymin, ymax, 161)
    X, Y = np.meshgrid(xs, ys)
    Z = terrain.h(X, Y)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=2, cstride=2, linewidth=0, alpha=0.8)

    ax.plot(Xb, Yb, Zb, linewidth=2.5, label='base')
    ax.plot(Xl, Yl, Zl, linewidth=1.8, label='left')
    ax.plot(Xr, Yr, Zr, linewidth=1.8, label='right')

    # 标明起点和终点
    ax.scatter([Xb[0]], [Yb[0]], [Zb[0]], color='g', s=60, marker='o', label='start')
    ax.text(Xb[0], Yb[0], Zb[0] + 0.02 * np.ptp(Z), "start", color='g')
    ax.scatter([Xb[-1]], [Yb[-1]], [Zb[-1]], color='r', s=80, marker='X', label='goal')
    ax.text(Xb[-1], Yb[-1], Zb[-1] + 0.02 * np.ptp(Z), "goal", color='r')

    # 用箭头表示偏航角 Psib，放在基座轨迹 (Xb, Yb, Zb) 上，按一定间隔绘制以免过密
    n_pts = len(Xb)
    step = max(1, n_pts // 20)  # 最多显示约20个箭头
    u = np.cos(Psib[::step])
    v = np.sin(Psib[::step])
    w = np.zeros_like(u)
    # 箭头长度相对于 x,y 范围
    arrow_len = 0.08 * max(np.ptp(xs), np.ptp(ys))
    # 在基座轨迹上绘制箭头（不再贴地，而是直接用 Zb 上的高度）
    zs = Zb[::step]
    ax.quiver(Xb[::step], Yb[::step], zs, u, v, w, length=arrow_len, normalize=True, color='k', linewidth=1)

    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title('Power-basis trajectories with boundary constraints on RBF terrain')
    ax.legend(loc='upper left')

    xr, yr, zr = np.ptp(xs), np.ptp(ys), np.ptp(Z)
    ax.set_box_aspect((xr, yr, zr))
    ax.view_init(elev=35, azim=-60)
    plt.tight_layout()

    fig2, axs = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
    axs[0].plot(Xb, label='Xb'); axs[0].set_ylabel('Xb'); axs[0].legend()
    axs[1].plot(Yb, label='Yb'); axs[1].set_ylabel('Yb'); axs[1].legend()
    axs[2].plot(Zb, label='Zb'); axs[2].set_ylabel('Zb'); axs[2].set_xlabel('sample index'); axs[2].legend()
    axs[2].set_ylim(0.0, 1.2)
    fig2.suptitle('Base trajectory components')
    plt.tight_layout()
    plt.show()
