#!/usr/bin/env python3
"""
轨迹可视化脚本
从C++程序生成的数据文件读取并可视化轨迹和地形
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import time

def visualize_trajectory(traj_file='/tmp/trajectory_data.txt',
                         terrain_file='/tmp/rbf_terrain_data.txt',
                         auto_close=False, close_time=5.0):
    """从数据文件读取并可视化轨迹和地形"""
    
    # 读取轨迹数据
    if not os.path.exists(traj_file):
        print(f"Error: Trajectory file {traj_file} not found!")
        return
    
    print(f"Reading trajectory data from: {traj_file}")
    with open(traj_file, 'r') as f:
        n_points = int(f.readline().strip())
        print(f"  Number of points: {n_points}")
        
        Xb, Yb, Zb, Psib = [], [], [], []
        Xl, Yl, Zl = [], [], []
        Xr, Yr, Zr = [], [], []
        airborne_l, airborne_r = [], []  # 轮子离地距离
        
        for i in range(n_points):
            line = f.readline().strip().split()
            Xb.append(float(line[0]))
            Yb.append(float(line[1]))
            Zb.append(float(line[2]))
            Psib.append(float(line[3]))
            Xl.append(float(line[4]))
            Yl.append(float(line[5]))
            Zl.append(float(line[6]))
            Xr.append(float(line[7]))
            Yr.append(float(line[8]))
            Zr.append(float(line[9]))
            # 如果有离地距离数据（新格式）
            if len(line) >= 12:
                airborne_l.append(float(line[10]))
                airborne_r.append(float(line[11]))
            else:
                # 旧格式，计算离地距离（需要地形数据）
                airborne_l.append(0.0)
                airborne_r.append(0.0)
    
    Xb, Yb, Zb, Psib = np.array(Xb), np.array(Yb), np.array(Zb), np.array(Psib)
    Xl, Yl, Zl = np.array(Xl), np.array(Yl), np.array(Zl)
    Xr, Yr, Zr = np.array(Xr), np.array(Yr), np.array(Zr)
    airborne_l = np.array(airborne_l)
    airborne_r = np.array(airborne_r)
    
    # 如果没有离地距离数据，尝试从地形计算（使用简单的最近邻插值）
    if np.all(airborne_l == 0.0) and np.all(airborne_r == 0.0) and terrain_data is not None:
        print("  Computing airborne distances from terrain data...")
        X_terrain, Y_terrain, Z_terrain, _, _, _, _ = terrain_data
        # 使用简单的最近邻插值（避免依赖scipy）
        Hl = np.zeros_like(Zl)
        Hr = np.zeros_like(Zr)
        for i in range(n_points):
            # 找到最近的网格点
            idx_xl = np.argmin(np.abs(X_terrain[0, :] - Xl[i]))
            idx_yl = np.argmin(np.abs(Y_terrain[:, 0] - Yl[i]))
            Hl[i] = Z_terrain[idx_yl, idx_xl]
            idx_xr = np.argmin(np.abs(X_terrain[0, :] - Xr[i]))
            idx_yr = np.argmin(np.abs(Y_terrain[:, 0] - Yr[i]))
            Hr[i] = Z_terrain[idx_yr, idx_xr]
        airborne_l = Zl - Hl
        airborne_r = Zr - Hr
    
    # 读取地形数据
    terrain_data = None
    if os.path.exists(terrain_file):
        print(f"Reading terrain data from: {terrain_file}")
        with open(terrain_file, 'r') as f:
            nx, ny = map(int, f.readline().split())
            xmin, xmax, ymin, ymax = map(float, f.readline().split())
            xs = np.array(list(map(float, f.readline().split())))
            ys = np.array(list(map(float, f.readline().split())))
            Z_terrain = np.zeros((ny, nx))
            for j in range(ny):
                Z_terrain[j, :] = np.array(list(map(float, f.readline().split())))
        
        X_terrain, Y_terrain = np.meshgrid(xs, ys)
        terrain_data = (X_terrain, Y_terrain, Z_terrain, xmin, xmax, ymin, ymax)
        print(f"  Terrain grid: {nx} x {ny}")
    
    # 创建图形（移除统计信息图，只保留3个子图）
    fig = plt.figure(figsize=(18, 6))
    
    # 3D轨迹图
    ax1 = fig.add_subplot(131, projection='3d')
    
    # 绘制地形（如果可用）
    if terrain_data is not None:
        X_terrain, Y_terrain, Z_terrain, _, _, _, _ = terrain_data
        ax1.plot_surface(X_terrain, Y_terrain, Z_terrain, 
                        cmap='viridis', alpha=0.3, linewidth=0, antialiased=False)
    
    # 绘制轨迹
    ax1.plot(Xb, Yb, Zb, 'b-', linewidth=2.5, label='Base trajectory', zorder=10)
    ax1.plot(Xl, Yl, Zl, 'g-', linewidth=1.8, label='Left wheel', zorder=10)
    ax1.plot(Xr, Yr, Zr, 'r-', linewidth=1.8, label='Right wheel', zorder=10)
    
    # 标记起点和终点
    ax1.scatter([Xb[0]], [Yb[0]], [Zb[0]], color='green', s=100, marker='o', 
                label='Start', zorder=20)
    ax1.scatter([Xb[-1]], [Yb[-1]], [Zb[-1]], color='red', s=100, marker='X', 
                label='Goal', zorder=20)
    
    # 绘制偏航角箭头（每隔几个点画一个）
    step = max(1, n_points // 20)
    for i in range(0, n_points, step):
        u = np.cos(Psib[i])
        v = np.sin(Psib[i])
        w = 0.0
        arrow_len = 0.1
        ax1.quiver(Xb[i], Yb[i], Zb[i], u, v, w, 
                  length=arrow_len, normalize=True, color='black', 
                  arrow_length_ratio=0.3, linewidth=1, zorder=15)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.set_title('Trajectory on RBF Terrain (3D)')
    ax1.legend(loc='upper left')
    
    # 2D俯视图
    ax2 = fig.add_subplot(132)
    
    # 绘制地形等高线（如果可用）
    if terrain_data is not None:
        X_terrain, Y_terrain, Z_terrain, _, _, _, _ = terrain_data
        contour = ax2.contourf(X_terrain, Y_terrain, Z_terrain, levels=20, 
                               cmap='viridis', alpha=0.5)
        plt.colorbar(contour, ax=ax2, shrink=0.8)
    
    # 绘制轨迹
    ax2.plot(Xb, Yb, 'b-', linewidth=2.5, label='Base trajectory')
    ax2.plot(Xl, Yl, 'g-', linewidth=1.8, label='Left wheel')
    ax2.plot(Xr, Yr, 'r-', linewidth=1.8, label='Right wheel')
    
    # 标记起点和终点
    ax2.scatter([Xb[0]], [Yb[0]], color='green', s=100, marker='o', 
                label='Start', zorder=20)
    ax2.scatter([Xb[-1]], [Yb[-1]], color='red', s=100, marker='X', 
                label='Goal', zorder=20)
    
    # 绘制偏航角箭头
    for i in range(0, n_points, step):
        u = np.cos(Psib[i])
        v = np.sin(Psib[i])
        arrow_len = 0.15
        ax2.arrow(Xb[i], Yb[i], u * arrow_len, v * arrow_len,
                 head_width=0.05, head_length=0.05, fc='black', ec='black')
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title('Trajectory (Top View)')
    ax2.set_aspect('equal')
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    # 轮子离地距离图
    ax3 = fig.add_subplot(133)
    s = np.linspace(0, 1, n_points)  # 归一化的路径参数
    ax3.plot(s, airborne_l * 100, 'g-', linewidth=2, label='Left wheel (cm)', alpha=0.7)
    ax3.plot(s, airborne_r * 100, 'r-', linewidth=2, label='Right wheel (cm)', alpha=0.7)
    ax3.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5, label='Ground level')
    ax3.axhline(y=5, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='5cm threshold')
    ax3.axhline(y=-5, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    ax3.fill_between(s, 0, 5, alpha=0.1, color='green', label='Allowed range (±5cm)')
    ax3.fill_between(s, -5, 0, alpha=0.1, color='green')
    ax3.set_xlabel('Path parameter s')
    ax3.set_ylabel('Airborne distance (cm)')
    ax3.set_title('Wheel Airborne Distance Over Time')
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3)
    
    # 统计信息可视化已移除（根据用户要求）
    
    plt.tight_layout()
    
    if auto_close:
        print(f"Visualization will auto-close in {close_time} seconds...")
        plt.show(block=False)
        # 先短暂暂停以确保窗口完全渲染
        plt.pause(0.1)
        # 然后等待指定时间
        plt.pause(close_time)
        plt.close()
        print("Visualization closed automatically.")
    else:
        plt.show()
    
    print("Visualization complete!")

if __name__ == '__main__':
    traj_file = '/tmp/trajectory_data.txt'
    terrain_file = '/tmp/rbf_terrain_data.txt'
    auto_close = False
    close_time = 5.0
    
    # 解析命令行参数
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--auto-close':
            auto_close = True
            # 检查下一个参数是否是时间值
            if i + 1 < len(sys.argv):
                try:
                    close_time = float(sys.argv[i + 1])
                    i += 2  # 跳过时间参数
                    continue
                except (ValueError, IndexError):
                    pass
            i += 1
        elif i == 1:
            # 第一个参数是轨迹文件
            traj_file = arg
            i += 1
        elif i == 2:
            # 第二个参数是地形文件
            terrain_file = arg
            i += 1
        else:
            i += 1
    
    visualize_trajectory(traj_file, terrain_file, auto_close, close_time)

