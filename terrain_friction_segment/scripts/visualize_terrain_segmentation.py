#!/usr/bin/env python3
"""
Terrain vs friction-cone segmentation comparison (called by terrain_friction_segment).
- Panel 1: terrain height
- Panel 2: median-smoothed binary segmentation
- Panel 3: morphologically smoothed (opening + closing) to smooth rollable/lift boundaries.
  Opening (erosion then dilation) removes small islands; closing (dilation then erosion) fills small holes.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    from scipy.ndimage import median_filter, binary_closing, binary_opening
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Median smoothing kernel (odd)
SMOOTH_KERNEL = 5
# Morphological structure size for opening/closing (e.g. 5 = 5x5 square)
MORPH_SIZE = 11
# Use light closing in median step
USE_CLOSING = False


def smooth_segmentation(labels, kernel_size=5, use_closing=True):
    """Median filter + optional binary closing."""
    if kernel_size < 2:
        return labels
    if HAS_SCIPY:
        smoothed = median_filter(labels.astype(np.float64), size=kernel_size)
        out = (smoothed > 0.5).astype(np.int32)
        if use_closing:
            out = binary_closing(out, structure=np.ones((3, 3))).astype(np.int32)
        return out
    ny, nx = labels.shape
    k = kernel_size // 2
    out = np.zeros_like(labels)
    for j in range(ny):
        for i in range(nx):
            j0, j1 = max(0, j - k), min(ny, j + k + 1)
            i0, i1 = max(0, i - k), min(nx, i + k + 1)
            window = labels[j0:j1, i0:i1]
            out[j, i] = 1 if np.sum(window) > window.size / 2 else 0
    return out


def morphological_smooth(labels, size=5):
    """
    Smooth binary 0/1 image with opening then closing.
    Opening (erosion then dilation): removes small islands, smooths boundaries inward.
    Closing (dilation then erosion): fills small holes, smooths boundaries outward.
    """
    if not HAS_SCIPY or size < 2:
        return labels
    structure = np.ones((size, size), dtype=np.int32)
    out = binary_opening(labels.astype(np.uint8), structure=structure)
    out = binary_closing(out, structure=structure)
    return out.astype(np.int32)


def visualize_terrain_segmentation(data_file='/tmp/terrain_segmentation_data.txt',
                                   auto_close=False, close_time=5.0):
    if not os.path.exists(data_file):
        print(f"Error: Data file {data_file} not found!")
        return 1
    print(f"Reading data from: {data_file}")
    with open(data_file, 'r') as f:
        line1 = f.readline().strip().split()
        nx, ny = int(line1[0]), int(line1[1])
        xmin, xmax = float(line1[2]), float(line1[3])
        ymin, ymax = float(line1[4]), float(line1[5])
        line2 = f.readline().strip().split()
        terrain = np.array([float(x) for x in line2], dtype=np.float64)
        line3 = f.readline().strip().split()
        labels_raw = np.array([int(x) for x in line3 if x != ''], dtype=np.int32)
    if terrain.size != nx * ny or labels_raw.size != nx * ny:
        print(f"Error: Expected {nx*ny} values, got terrain {terrain.size}, labels {labels_raw.size}")
        return 1
    terrain = terrain.reshape(ny, nx)
    labels_raw = labels_raw.reshape(ny, nx)

    labels_smooth = smooth_segmentation(labels_raw, kernel_size=SMOOTH_KERNEL, use_closing=USE_CLOSING)
    labels_morph = morphological_smooth(labels_smooth, size=MORPH_SIZE)
    extent = [xmin, xmax, ymin, ymax]

    print(f"Smooth k={SMOOTH_KERNEL}; morph opening+closing size={MORPH_SIZE}")

    n_plots = 3
    fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 6), sharex=True, sharey=True)
    cmap_seg = plt.matplotlib.colors.ListedColormap(['#c0392b', '#27ae60'])

    # Left: terrain
    ax1 = axes[0]
    im1 = ax1.imshow(terrain, extent=extent, origin='lower', cmap='terrain', interpolation='nearest')
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_title('Terrain (RBF height)')
    plt.colorbar(im1, ax=ax1, shrink=0.7, label='Height (m)')

    # Middle: median-smoothed segmentation
    ax2 = axes[1]
    im2 = ax2.imshow(labels_smooth, extent=extent, origin='lower', cmap=cmap_seg, vmin=0, vmax=1, interpolation='nearest')
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.set_title('Segmentation (median smoothed)')
    cbar2 = plt.colorbar(im2, ax=ax2, ticks=[0.25, 0.75], shrink=0.7)
    cbar2.ax.set_yticklabels(['Lift leg', 'Rollable'])

    # Right: morphological smoothing (opening + closing)
    ax3 = axes[2]
    im3 = ax3.imshow(labels_morph, extent=extent, origin='lower', cmap=cmap_seg, vmin=0, vmax=1, interpolation='nearest')
    ax3.set_xlabel('x (m)')
    ax3.set_ylabel('y (m)')
    ax3.set_title('Morphological (opening + closing)')
    cbar3 = plt.colorbar(im3, ax=ax3, ticks=[0.25, 0.75], shrink=0.7)
    cbar3.ax.set_yticklabels(['Lift leg', 'Rollable'])

    plt.tight_layout()
    if auto_close:
        plt.show(block=False)
        plt.pause(close_time)
        plt.close()
        print(f"Visualization closed automatically after {close_time}s.")
    else:
        plt.show()
    print("Visualization complete!")
    return 0


if __name__ == '__main__':
    data_file = '/tmp/terrain_segmentation_data.txt'
    auto_close = False
    close_time = 5.0
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
    sys.exit(visualize_terrain_segmentation(data_file, auto_close, close_time))
