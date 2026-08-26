#!/usr/bin/env python3
"""
可视化分段轨迹优化的初始解与最终解，叠加 RBF 地形。
用法: python3 plot_seg_trajectory.py [输出目录] [--no-show]
读取 输出目录/ 下的 initial_trajectory.txt、final_trajectory.txt、terrain_grid.txt（若有），
绘制 2D (xy) 与 3D (xyz) 轨迹对比图；另存 kinematics 图（对路径参数 s 的导数：|dr/ds|、||d²r/ds²||、dψ/ds）。
默认先弹出窗口可旋转 3D 视角，关闭窗口后保存到 输出目录/seg_trajectory_plot.png 与 seg_trajectory_kinematics.png。
加 --no-show 则直接保存不弹窗（适用于无界面环境）。
"""

import os
import sys
import numpy as np

# 无界面模式（--no-show）用 Agg，否则用默认后端以便 plt.show() 弹窗
if '--no-show' in sys.argv:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def load_terrain_grid(path):
    """读取地形网格文件，返回 (X, Y, Z) 网格数组，或 None。"""
    if not os.path.isfile(path):
        return None
    with open(path, 'r') as f:
        header = f.readline().strip().split()
        if len(header) != 6:
            return None
        nx, ny = int(header[0]), int(header[1])
        xmin, xmax = float(header[2]), float(header[3])
        ymin, ymax = float(header[4]), float(header[5])
        z_flat = []
        for line in f:
            line = line.strip()
            if line:
                z_flat.append(float(line))
    if len(z_flat) != nx * ny:
        return None
    Z = np.array(z_flat).reshape(ny, nx)
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(x, y)
    return X, Y, Z


def _deriv_wrt_s(y, s):
    """dy/ds，非均匀 s 用 np.gradient。"""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    if y.size < 2:
        return np.zeros_like(y)
    return np.gradient(y, s, edge_order=2)


def _speed_accel_3d(traj, xk, yk, zk):
    """
    3D 位置对路径参数 s 的一阶、二阶导数标量：
    |dr/ds| = ||(x',y',z')||,  ||d²r/ds²|| = ||(x'',y'',z'')||
    """
    s = traj['s']
    vx = _deriv_wrt_s(traj[xk], s)
    vy = _deriv_wrt_s(traj[yk], s)
    vz = _deriv_wrt_s(traj[zk], s)
    speed = np.sqrt(vx * vx + vy * vy + vz * vz)
    ax = _deriv_wrt_s(vx, s)
    ay = _deriv_wrt_s(vy, s)
    az = _deriv_wrt_s(vz, s)
    acc_mag = np.sqrt(ax * ax + ay * ay + az * az)
    return speed, acc_mag


def _unwrap_psi(psi):
    """偏航角展开，便于对 s 求导。"""
    p = np.asarray(psi, dtype=float)
    if p.size < 2:
        return p
    return np.unwrap(p)


def plot_kinematics_vs_s(fig, axes, init, final):
    """
    Base / L wheel / R wheel: |dr/ds| and ||d²r/ds²|| vs s; base yaw: dpsi/ds, d²psi/ds².
    axes: (4, 2) ndarray of Axes from subplots.
    init/final: trajectory dict or None.
    """

    def plot_pair(row, title_prefix, traj_init, traj_final, xk, yk, zk, is_psi=False):
        ax_a, ax_b = axes[row, 0], axes[row, 1]
        for tr, name, sty in (
            (traj_init, 'Initial', 'b--'),
            (traj_final, 'Final', 'r-'),
        ):
            if tr is None:
                continue
            s = tr['s']
            if is_psi:
                psi = _unwrap_psi(tr['psi'])
                d1 = _deriv_wrt_s(psi, s)
                d2 = _deriv_wrt_s(d1, s)
                ax_a.plot(s, d1, sty, linewidth=1.5, label=name, alpha=0.9)
                ax_b.plot(s, d2, sty, linewidth=1.5, label=name, alpha=0.9)
            else:
                sp, acc = _speed_accel_3d(tr, xk, yk, zk)
                ax_a.plot(s, sp, sty, linewidth=1.5, label=name, alpha=0.9)
                ax_b.plot(s, acc, sty, linewidth=1.5, label=name, alpha=0.9)
        ax_a.set_title(f'{title_prefix}: |dr/ds|' if not is_psi else f'{title_prefix}: dψ/ds')
        ax_a.set_xlabel('s')
        ax_a.set_ylabel('|dr/ds| (m / unit s)' if not is_psi else 'dψ/ds (rad / unit s)')
        ax_a.grid(True, alpha=0.3)
        ax_a.legend(fontsize=7)
        ax_b.set_title(f'{title_prefix}: ||d²r/ds²||' if not is_psi else f'{title_prefix}: d²ψ/ds²')
        ax_b.set_xlabel('s')
        ax_b.set_ylabel('||d²r/ds²|| (m / unit s²)' if not is_psi else 'd²ψ/ds² (rad / unit s²)')
        ax_b.grid(True, alpha=0.3)
        ax_b.legend(fontsize=7)

    plot_pair(0, 'Base', init, final, 'xb', 'yb', 'zb', is_psi=False)
    plot_pair(1, 'Left wheel', init, final, 'xl', 'yl', 'zl', is_psi=False)
    plot_pair(2, 'Right wheel', init, final, 'xr', 'yr', 'zr', is_psi=False)
    plot_pair(3, 'Base yaw ψ', init, final, None, None, None, is_psi=True)

    fig.suptitle(
        'Kinematics vs path parameter s (numerical d/ds; not physical time unless s ∝ t)',
        fontsize=10,
    )


def load_trajectory(path):
    """读取轨迹文件，返回 dict 含 s, xb, yb, zb, psi, xl, yl, zl, xr, yr, zr 的数组。"""
    if not os.path.isfile(path):
        return None
    data = []
    with open(path, 'r') as f:
        header = f.readline().strip()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 11:
                data.append([float(parts[i]) for i in range(11)])
    if not data:
        return None
    arr = np.array(data)
    return {
        's': arr[:, 0], 'xb': arr[:, 1], 'yb': arr[:, 2], 'zb': arr[:, 3], 'psi': arr[:, 4],
        'xl': arr[:, 5], 'yl': arr[:, 6], 'zl': arr[:, 7],
        'xr': arr[:, 8], 'yr': arr[:, 9], 'zr': arr[:, 10],
    }


def main():
    args = [a for a in sys.argv[1:] if a != '--no-show']
    out_dir = args[0] if args else 'seg_out'
    do_show = '--no-show' not in sys.argv
    if not out_dir.endswith('/'):
        out_dir += '/'
    init_path = out_dir + 'initial_trajectory.txt'
    final_path = out_dir + 'final_trajectory.txt'
    terrain_path = out_dir + 'terrain_grid.txt'
    if not os.path.isfile(init_path) and not os.path.isfile(final_path) and os.path.isdir('build'):
        build_out = 'build/' + out_dir
        if os.path.isfile(build_out + 'initial_trajectory.txt') or os.path.isfile(build_out + 'final_trajectory.txt'):
            out_dir = build_out
            init_path = out_dir + 'initial_trajectory.txt'
            final_path = out_dir + 'final_trajectory.txt'
            terrain_path = out_dir + 'terrain_grid.txt'

    terrain = load_terrain_grid(terrain_path)
    init = load_trajectory(init_path)
    final = load_trajectory(final_path)
    if init is None and final is None:
        print('未找到轨迹文件:', init_path, final_path)
        return 1
    if init is None:
        init = final
    if final is None:
        final = init

    fig = plt.figure(figsize=(14, 6))

    # 2D xy + terrain
    ax1 = fig.add_subplot(1, 2, 1)
    if terrain is not None:
        X, Y, Z = terrain
        cf = ax1.contourf(X, Y, Z, levels=14, cmap='terrain', alpha=0.6)
        plt.colorbar(cf, ax=ax1, label='Terrain z (m)', shrink=0.7)
    if init is not None:
        ax1.plot(init['xb'], init['yb'], 'b-o', markersize=3, label='Initial (base)', alpha=0.9)
        ax1.plot(init['xl'], init['yl'], 'b--', linewidth=1, alpha=0.6)
        ax1.plot(init['xr'], init['yr'], 'b--', linewidth=1, alpha=0.6)
    if final is not None:
        ax1.plot(final['xb'], final['yb'], 'r-s', markersize=3, label='Final (base)', alpha=0.9)
        ax1.plot(final['xl'], final['yl'], 'r-.', linewidth=1, alpha=0.6)
        ax1.plot(final['xr'], final['yr'], 'r-.', linewidth=1, alpha=0.6)
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_title('Trajectory comparison (xy) + RBF terrain')
    ax1.legend(loc='best', fontsize=8)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)

    # 3D: RBF 地形 + 基座 + 左轮、右轮路径
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    if terrain is not None:
        X, Y, Z = terrain
        ax2.plot_surface(X, Y, Z, cmap='terrain', alpha=0.5, rstride=2, cstride=2, label='RBF terrain')
    if init is not None:
        ax2.plot(init['xb'], init['yb'], init['zb'], 'b-o', markersize=2, label='Initial (base)', alpha=0.9)
        ax2.plot(init['xl'], init['yl'], init['zl'], 'c--', linewidth=1.2, label='Initial (L wheel)', alpha=0.8)
        ax2.plot(init['xr'], init['yr'], init['zr'], 'c-.', linewidth=1.2, label='Initial (R wheel)', alpha=0.8)
    if final is not None:
        ax2.plot(final['xb'], final['yb'], final['zb'], 'r-s', markersize=2, label='Final (base)', alpha=0.9)
        ax2.plot(final['xl'], final['yl'], final['zl'], 'm--', linewidth=1.2, label='Final (L wheel)', alpha=0.8)
        ax2.plot(final['xr'], final['yr'], final['zr'], 'm-.', linewidth=1.2, label='Final (R wheel)', alpha=0.8)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.set_zlabel('z (m)')
    ax2.set_title('3D: RBF terrain + base & wheels')
    ax2.legend(loc='best', fontsize=7)

    plt.tight_layout()
    out_png = out_dir + 'seg_trajectory_plot.png'
    fig.savefig(out_png, dpi=120)
    if do_show:
        plt.show()
    plt.close(fig)
    print('已保存:', out_png)

    # Kinematics: velocity/acceleration w.r.t. path parameter s
    fig_k, axes_k = plt.subplots(4, 2, figsize=(12, 14), constrained_layout=True)
    plot_kinematics_vs_s(fig_k, axes_k, init, final)
    out_kin = out_dir + 'seg_trajectory_kinematics.png'
    fig_k.savefig(out_kin, dpi=120)
    if do_show:
        plt.show()
    plt.close(fig_k)
    print('已保存:', out_kin)
    return 0


if __name__ == '__main__':
    sys.exit(main())
