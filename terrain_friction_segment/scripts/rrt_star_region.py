#!/usr/bin/env python3
"""
Minimal RRT* with region cost: read segmentation, plan from start to goal.
Cost = length in rollable (green) + k * length in lift (red). Prefer green.
Usage: python rrt_star_region.py [segmentation_data.txt] [start_x start_y goal_x goal_y]
Default: /tmp/terrain_segmentation_data.txt, start=(0.6, 1.2), goal=(2.4, -3.0)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Cost weight for lift-leg region (k > 1)
LIFT_COST_WEIGHT = 5.0
# RRT* params
MAX_ITER = 2000
STEP_LENGTH = 0.15
GOAL_RADIUS = 0.2
REWIRE_RADIUS = 0.5
# Number of samples along edge for cost integration
EDGE_SAMPLES = 20


def load_segmentation(path):
    """Load segmentation file; return (xmin, xmax, ymin, ymax, labels_2d)."""
    with open(path, 'r') as f:
        line1 = f.readline().strip().split()
        nx, ny = int(line1[0]), int(line1[1])
        xmin, xmax = float(line1[2]), float(line1[3])
        ymin, ymax = float(line1[4]), float(line1[5])
        f.readline()  # terrain line
        line3 = f.readline().strip().split()
        labels = np.array([int(x) for x in line3 if x != ''], dtype=np.int32)
    labels = labels.reshape(ny, nx)
    return xmin, xmax, ymin, ymax, labels


def world_to_grid(x, y, xmin, xmax, ymin, ymax, nx, ny):
    """Map world (x,y) to grid indices (i, j); clamp to valid range."""
    if xmax == xmin:
        i = 0
    else:
        i = (x - xmin) / (xmax - xmin) * (nx - 1)
    if ymax == ymin:
        j = 0
    else:
        j = (y - ymin) / (ymax - ymin) * (ny - 1)
    i = int(round(np.clip(i, 0, nx - 1)))
    j = int(round(np.clip(j, 0, ny - 1)))
    return i, j


def get_label(x, y, xmin, xmax, ymin, ymax, labels):
    """Return 1 (rollable) or 0 (lift) at world (x,y)."""
    ny, nx = labels.shape
    i, j = world_to_grid(x, y, xmin, xmax, ymin, ymax, nx, ny)
    return labels[j, i]


def edge_cost(x1, y1, x2, y2, xmin, xmax, ymin, ymax, labels, k):
    """Cost along segment from (x1,y1) to (x2,y2): length * (fraction_rollable + k*fraction_lift)."""
    n = EDGE_SAMPLES
    cost = 0.0
    dx = (x2 - x1) / n
    dy = (y2 - y1) / n
    seg_len = np.hypot(x2 - x1, y2 - y1) / n
    for i in range(n):
        x = x1 + (i + 0.5) * dx
        y = y1 + (i + 0.5) * dy
        L = get_label(x, y, xmin, xmax, ymin, ymax, labels)
        cost += seg_len * (1.0 if L == 1 else k)
    return cost


def distance(a, b):
    return np.hypot(b[0] - a[0], b[1] - a[1])


def steer(from_pt, to_pt, step):
    """Return point at most step away from from_pt towards to_pt."""
    d = distance(from_pt, to_pt)
    if d <= step:
        return to_pt[0], to_pt[1]
    t = step / d
    return from_pt[0] + t * (to_pt[0] - from_pt[0]), from_pt[1] + t * (to_pt[1] - from_pt[1])


def rrt_star(start, goal, bounds, labels, k=LIFT_COST_WEIGHT, max_iter=MAX_ITER,
             step=STEP_LENGTH, goal_radius=GOAL_RADIUS, rewire_radius=REWIRE_RADIUS):
    """
    RRT* with region cost. bounds = (xmin, xmax, ymin, ymax).
    Tree: list of (x, y), parent index, cost from root.
    """
    xmin, xmax, ymin, ymax = bounds
    nodes = [start]
    parent = [-1]
    cost_from_root = [0.0]

    def cost(a, b):
        return edge_cost(a[0], a[1], b[0], b[1], xmin, xmax, ymin, ymax, labels, k)

    path_to_goal = None
    goal_idx = -1

    for _ in range(max_iter):
        # Sample (bias 10% to goal)
        if np.random.rand() < 0.1:
            rnd = goal
        else:
            rnd = (np.random.uniform(xmin, xmax), np.random.uniform(ymin, ymax))

        # Nearest
        nearest_idx = min(range(len(nodes)), key=lambda i: distance(nodes[i], rnd))
        nearest = nodes[nearest_idx]
        new_pt = steer(nearest, rnd, step)
        c_edge = cost(nearest, new_pt)
        new_cost = cost_from_root[nearest_idx] + c_edge

        # Check bounds
        if not (xmin <= new_pt[0] <= xmax and ymin <= new_pt[1] <= ymax):
            continue

        # Choose parent: nearest or better in rewire radius
        for i in range(len(nodes)):
            if distance(nodes[i], new_pt) <= rewire_radius:
                c = cost(nodes[i], new_pt)
                if cost_from_root[i] + c < new_cost:
                    new_cost = cost_from_root[i] + c
                    nearest_idx = i

        nodes.append(new_pt)
        parent.append(nearest_idx)
        cost_from_root.append(new_cost)

        # Rewire: nodes that could lower cost via new node
        new_idx = len(nodes) - 1
        for i in range(len(nodes) - 1):
            if distance(nodes[i], new_pt) <= rewire_radius:
                c = cost(new_pt, nodes[i])
                if new_cost + c < cost_from_root[i]:
                    parent[i] = new_idx
                    cost_from_root[i] = new_cost + c

        # Goal check
        if distance(new_pt, goal) <= goal_radius:
            c_end = cost(new_pt, goal)
            total = new_cost + c_end
            if goal_idx < 0 or total < cost_from_root[goal_idx]:
                goal_idx = len(nodes)  # index of goal node we are about to add
                nodes.append(goal)
                parent.append(new_idx)
                cost_from_root.append(total)
                # Build path from start to goal
                path = []
                idx = goal_idx
                while idx >= 0:
                    path.append(nodes[idx])
                    idx = parent[idx]
                path.reverse()
                path_to_goal = path

    return path_to_goal, nodes, parent


def main():
    data_file = '/tmp/terrain_segmentation_data.txt'
    start = (0.6, 1.2)
    goal = (2.4, -3.0)
    if len(sys.argv) >= 2:
        data_file = sys.argv[1]
    if len(sys.argv) >= 6:
        start = (float(sys.argv[2]), float(sys.argv[3]))
        goal = (float(sys.argv[4]), float(sys.argv[5]))

    if not os.path.exists(data_file):
        print(f"Error: {data_file} not found. Run terrain_friction_segment first.")
        return 1

    xmin, xmax, ymin, ymax, labels = load_segmentation(data_file)
    ny, nx = labels.shape
    bounds = (xmin, xmax, ymin, ymax)
    extent = [xmin, xmax, ymin, ymax]

    print(f"Segmentation: {nx}x{ny}, bounds x=[{xmin},{xmax}] y=[{ymin},{ymax}]")
    print(f"Start {start}, goal {goal}, lift cost weight k={LIFT_COST_WEIGHT}")
    path, nodes, parent = rrt_star(start, goal, bounds, labels, k=LIFT_COST_WEIGHT)
    if path is None:
        print("No path found (increase MAX_ITER or check bounds).")
        path = []
    else:
        print(f"Path: {len(path)} waypoints")

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    cmap = plt.matplotlib.colors.ListedColormap(['#c0392b', '#27ae60'])
    ax.imshow(labels, extent=extent, origin='lower', cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
    ax.scatter(*start, c='lime', s=120, marker='o', edgecolors='black', linewidths=2, label='Start')
    ax.scatter(*goal, c='orange', s=120, marker='*', edgecolors='black', linewidths=2, label='Goal')
    if path:
        path = np.array(path)
        ax.plot(path[:, 0], path[:, 1], 'b-', linewidth=2.5, label='RRT* path')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title('RRT* with region cost (green=rollable, red=lift)')
    ax.legend(loc='upper right')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()
    return 0


if __name__ == '__main__':
    sys.exit(main())
