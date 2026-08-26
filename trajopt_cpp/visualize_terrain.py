#!/usr/bin/env python3
"""
RBF地形可视化脚本
从C++程序生成的数据文件读取并可视化RBF地形
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import time

def visualize_terrain(data_file='/tmp/rbf_terrain_data.txt', auto_close=False, close_time=5.0):
    """从数据文件读取并可视化RBF地形"""
    if not os.path.exists(data_file):
        print(f"Error: Data file {data_file} not found!")
        print("Please run the C++ program first to generate the terrain data.")
        return
    
    print(f"Reading terrain data from: {data_file}")
    with open(data_file, 'r') as f:
        # 读取网格尺寸
        nx, ny = map(int, f.readline().split())
        print(f"  Grid size: {nx} x {ny}")
        
        # 读取范围
        xmin, xmax, ymin, ymax = map(float, f.readline().split())
        print(f"  Domain: x=[{xmin:.2f}, {xmax:.2f}], y=[{ymin:.2f}, {ymax:.2f}]")
        
        # 读取x坐标
        xs = np.array(list(map(float, f.readline().split())))
        
        # 读取y坐标
        ys = np.array(list(map(float, f.readline().split())))
        
        # 读取高度数据
        Z = np.zeros((ny, nx))
        for j in range(ny):
            Z[j, :] = np.array(list(map(float, f.readline().split())))
    
    # 创建网格
    X, Y = np.meshgrid(xs, ys)
    
    # 创建3D表面图
    fig = plt.figure(figsize=(12, 8))
    
    # 3D表面图
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, 
                            antialiased=False, alpha=0.9)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z (height)')
    ax1.set_title('RBF Terrain (3D Surface)')
    plt.colorbar(surf, ax=ax1, shrink=0.5)
    
    # 2D等高线图
    ax2 = fig.add_subplot(122)
    contour = ax2.contourf(X, Y, Z, levels=20, cmap='viridis')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title('RBF Terrain (2D Contour)')
    ax2.set_aspect('equal')
    plt.colorbar(contour, ax=ax2)
    
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
    data_file = '/tmp/rbf_terrain_data.txt'
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
            # 第一个参数是数据文件
            data_file = arg
            i += 1
        else:
            i += 1
    
    visualize_terrain(data_file, auto_close, close_time)

