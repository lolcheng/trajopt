#!/usr/bin/env python3
"""
接触力、力平衡、力矩平衡、摩擦锥可视化脚本
从C++程序生成的数据文件读取并可视化
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time

def visualize_forces(data_file='/tmp/force_data.txt', auto_close=False, close_time=5.0):
    """从数据文件读取并可视化接触力、力平衡、力矩平衡、摩擦锥"""
    if not os.path.exists(data_file):
        print(f"Error: Data file {data_file} not found!")
        return
    
    print(f"Reading force data from: {data_file}")
    with open(data_file, 'r') as f:
        n_sample = int(f.readline().strip())
        print(f"  Number of sample points: {n_sample}")
        
        # 读取接触力
        fLx = np.array(list(map(float, f.readline().split())))
        fLy = np.array(list(map(float, f.readline().split())))
        fLz = np.array(list(map(float, f.readline().split())))
        fRx = np.array(list(map(float, f.readline().split())))
        fRy = np.array(list(map(float, f.readline().split())))
        fRz = np.array(list(map(float, f.readline().split())))
        
        # 读取轨迹数据（用于计算力矩）
        Xb = np.array(list(map(float, f.readline().split())))
        Yb = np.array(list(map(float, f.readline().split())))
        Zb = np.array(list(map(float, f.readline().split())))
        Xl = np.array(list(map(float, f.readline().split())))
        Yl = np.array(list(map(float, f.readline().split())))
        Zl = np.array(list(map(float, f.readline().split())))
        Xr = np.array(list(map(float, f.readline().split())))
        Yr = np.array(list(map(float, f.readline().split())))
        Zr = np.array(list(map(float, f.readline().split())))
        
        # 读取地形法向量
        nLx = np.array(list(map(float, f.readline().split())))
        nLy = np.array(list(map(float, f.readline().split())))
        nLz = np.array(list(map(float, f.readline().split())))
        nRx = np.array(list(map(float, f.readline().split())))
        nRy = np.array(list(map(float, f.readline().split())))
        nRz = np.array(list(map(float, f.readline().split())))
        
        # 读取参数
        MASS = float(f.readline().strip())
        GRAVITY = float(f.readline().strip())
        MU_FRICTION = float(f.readline().strip())
    
    # 采样点 s ∈ [0, 1]
    s = np.linspace(0.0, 1.0, n_sample)
    
    # 计算法向力和切向力
    fLn = fLx * nLx + fLy * nLy + fLz * nLz  # 左轮法向力
    fRn = fRx * nRx + fRy * nRy + fRz * nRz  # 右轮法向力
    
    fLtx = fLx - fLn * nLx  # 左轮切向力x
    fLty = fLy - fLn * nLy  # 左轮切向力y
    fLtz = fLz - fLn * nLz  # 左轮切向力z
    fLt_mag = np.sqrt(fLtx**2 + fLty**2 + fLtz**2)  # 左轮切向力大小
    
    fRtx = fRx - fRn * nRx  # 右轮切向力x
    fRty = fRy - fRn * nRy  # 右轮切向力y
    fRtz = fRz - fRn * nRz  # 右轮切向力z
    fRt_mag = np.sqrt(fRtx**2 + fRty**2 + fRtz**2)  # 右轮切向力大小
    
    # 计算力平衡残差
    Fx_res = fLx + fRx
    Fy_res = fLy + fRy
    Fz_res = fLz + fRz - MASS * GRAVITY
    
    # 计算力矩平衡残差
    rLx = Xl - Xb
    rLy = Yl - Yb
    rLz = Zl - Zb
    rRx = Xr - Xb
    rRy = Yr - Yb
    rRz = Zr - Zb
    
    tauLx = rLy * fLz - rLz * fLy
    tauLy = rLz * fLx - rLx * fLz
    tauLz = rLx * fLy - rLy * fLx
    
    tauRx = rRy * fRz - rRz * fRy
    tauRy = rRz * fRx - rRx * fRz
    tauRz = rRx * fRy - rRy * fRx
    
    tau_x = tauLx + tauRx
    tau_y = tauLy + tauRy
    tau_z = tauLz + tauRz
    
    # 计算摩擦锥违反
    mu2 = MU_FRICTION**2
    fric_viol_L = np.maximum(0, fLt_mag**2 - mu2 * fLn**2)
    fric_viol_R = np.maximum(0, fRt_mag**2 - mu2 * fRn**2)
    
    # 创建图形
    fig = plt.figure(figsize=(16, 12))
    
    # 1. 接触力可视化
    ax1 = fig.add_subplot(3, 2, 1)
    ax1.plot(s, fLx, 'r-', label='Left Fx', linewidth=2)
    ax1.plot(s, fLy, 'g-', label='Left Fy', linewidth=2)
    ax1.plot(s, fLz, 'b-', label='Left Fz', linewidth=2)
    ax1.plot(s, fRx, 'r--', label='Right Fx', linewidth=2)
    ax1.plot(s, fRy, 'g--', label='Right Fy', linewidth=2)
    ax1.plot(s, fRz, 'b--', label='Right Fz', linewidth=2)
    ax1.set_xlabel('s')
    ax1.set_ylabel('Force (N)')
    ax1.set_title('Contact Forces')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. 法向力和切向力
    ax2 = fig.add_subplot(3, 2, 2)
    ax2.plot(s, fLn, 'r-', label='Left Normal', linewidth=2)
    ax2.plot(s, fRn, 'r--', label='Right Normal', linewidth=2)
    ax2.plot(s, fLt_mag, 'b-', label='Left Tangential', linewidth=2)
    ax2.plot(s, fRt_mag, 'b--', label='Right Tangential', linewidth=2)
    ax2.set_xlabel('s')
    ax2.set_ylabel('Force (N)')
    ax2.set_title('Normal and Tangential Forces')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. 力平衡残差
    ax3 = fig.add_subplot(3, 2, 3)
    ax3.plot(s, Fx_res, 'r-', label='Fx residual', linewidth=2)
    ax3.plot(s, Fy_res, 'g-', label='Fy residual', linewidth=2)
    ax3.plot(s, Fz_res, 'b-', label='Fz residual', linewidth=2)
    ax3.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax3.set_xlabel('s')
    ax3.set_ylabel('Force Residual (N)')
    ax3.set_title('Force Balance Residuals')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. 力矩平衡残差
    ax4 = fig.add_subplot(3, 2, 4)
    ax4.plot(s, tau_x, 'r-', label='τx residual', linewidth=2)
    ax4.plot(s, tau_y, 'g-', label='τy residual', linewidth=2)
    ax4.plot(s, tau_z, 'b-', label='τz residual', linewidth=2)
    ax4.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax4.set_xlabel('s')
    ax4.set_ylabel('Torque Residual (N·m)')
    ax4.set_title('Torque Balance Residuals')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. 摩擦锥可视化
    ax5 = fig.add_subplot(3, 2, 5)
    # 绘制摩擦锥边界
    fLn_plot = np.linspace(0, np.max([np.max(fLn), np.max(fRn)]) * 1.1, 100)
    fLt_max = MU_FRICTION * fLn_plot
    ax5.plot(fLn_plot, fLt_max, 'k--', label=f'Friction cone (μ={MU_FRICTION})', linewidth=2)
    ax5.plot(fLn_plot, -fLt_max, 'k--', linewidth=2)
    # 绘制实际数据点
    ax5.scatter(fLn, fLt_mag, c='r', s=30, alpha=0.6, label='Left wheel', zorder=5)
    ax5.scatter(fRn, fRt_mag, c='b', s=30, alpha=0.6, marker='^', label='Right wheel', zorder=5)
    ax5.set_xlabel('Normal Force (N)')
    ax5.set_ylabel('Tangential Force Magnitude (N)')
    ax5.set_title('Friction Cone')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    ax5.set_aspect('equal', adjustable='box')
    
    # 6. 摩擦锥违反
    ax6 = fig.add_subplot(3, 2, 6)
    ax6.plot(s, fric_viol_L, 'r-', label='Left violation', linewidth=2)
    ax6.plot(s, fric_viol_R, 'b-', label='Right violation', linewidth=2)
    ax6.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax6.set_xlabel('s')
    ax6.set_ylabel('Friction Cone Violation')
    ax6.set_title('Friction Cone Violations')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if auto_close:
        print(f"Visualization will auto-close in {close_time} seconds...")
        plt.show(block=False)
        plt.pause(0.1)
        plt.pause(close_time)
        plt.close()
        print("Visualization closed automatically.")
    else:
        plt.show()
    
    print("Visualization complete!")

if __name__ == '__main__':
    data_file = '/tmp/force_data.txt'
    auto_close = False
    close_time = 5.0
    
    # 解析命令行参数
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--auto-close':
            auto_close = True
            if i + 1 < len(sys.argv):
                try:
                    close_time = float(sys.argv[i + 1])
                    i += 2
                    continue
                except (ValueError, IndexError):
                    pass
            i += 1
        elif i == 1:
            data_file = arg
            i += 1
        else:
            i += 1
    
    visualize_forces(data_file, auto_close, close_time)

