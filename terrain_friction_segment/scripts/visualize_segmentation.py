#!/usr/bin/env python3
"""
Visualize terrain_friction_segment output: green = rollable (friction cone OK), red = lift leg.
Usage: python scripts/visualize_segmentation.py [segmentation.txt]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        seg_path = Path(__file__).resolve().parent.parent / "segmentation.txt"
    else:
        seg_path = Path(sys.argv[1])
    if not seg_path.exists():
        print(f"File not found: {seg_path}")
        sys.exit(1)

    with open(seg_path) as f:
        first = f.readline().strip().split()
        second = f.readline().strip().split()
    nx, ny = int(first[0]), int(first[1])
    xmin, xmax = float(first[2]), float(first[3])
    ymin, ymax = float(first[4]), float(first[5])
    labels = np.array([int(x) for x in second if x], dtype=np.int32)
    if labels.size != nx * ny:
        print(f"Expected {nx*ny} values, got {labels.size}")
        sys.exit(1)
    grid = labels.reshape(ny, nx)

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    cmap = plt.matplotlib.colors.ListedColormap(["#c0392b", "#27ae60"])  # red, green
    im = ax.imshow(
        grid,
        extent=[xmin, xmax, ymin, ymax],
        origin="lower",
        cmap=cmap,
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Friction cone segmentation: green=rollable, red=lift leg")
    cbar = fig.colorbar(im, ax=ax, ticks=[0.25, 0.75], shrink=0.8)
    cbar.ax.set_yticklabels(["Lift leg", "Rollable"])
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
