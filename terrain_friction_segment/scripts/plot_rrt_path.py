#!/usr/bin/env python3
"""Plot segmentation + path from C++ RRT* output. Usage: python3 plot_rrt_path.py data.txt path.txt"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 plot_rrt_path.py <segmentation_data.txt> <path.txt>")
        return 1
    data_file = sys.argv[1]
    path_file = sys.argv[2]
    if not os.path.exists(data_file) or not os.path.exists(path_file):
        print("File not found.")
        return 1

    with open(data_file, 'r') as f:
        line1 = f.readline().strip().split()
        nx, ny = int(line1[0]), int(line1[1])
        xmin, xmax = float(line1[2]), float(line1[3])
        ymin, ymax = float(line1[4]), float(line1[5])
        f.readline()
        line3 = f.readline().strip().split()
        labels = np.array([int(x) for x in line3 if x != ''], dtype=np.int32).reshape(ny, nx)
    extent = [xmin, xmax, ymin, ymax]

    path = []
    with open(path_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                path.append([float(parts[0]), float(parts[1])])
    path = np.array(path) if path else None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    cmap = plt.matplotlib.colors.ListedColormap(['#c0392b', '#27ae60'])
    ax.imshow(labels, extent=extent, origin='lower', cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
    if path is not None and len(path) >= 2:
        ax.plot(path[:, 0], path[:, 1], 'b-', linewidth=2.5, label='RRT* path')
        ax.scatter(path[0, 0], path[0, 1], c='lime', s=120, marker='o', edgecolors='black', linewidths=2, label='Start')
        ax.scatter(path[-1, 0], path[-1, 1], c='orange', s=120, marker='*', edgecolors='black', linewidths=2, label='Goal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title('RRT* path (C++)')
    ax.legend(loc='upper right')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()
    return 0

if __name__ == '__main__':
    sys.exit(main())
