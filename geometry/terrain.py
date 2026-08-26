from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

Array = np.ndarray
ScalarOrArray = Union[float, int, Array]

"""
Simplified radial-basis terrain defined on (x,y) ∈ [-1,1]×[-1,1].

Height field:
    h(x,y) = amp * exp( - 0.5 * (x^2 + y^2) / ell^2 )

Its gradient and 3D unit normal are provided. The 3D normal corresponds to the
surface z = h(x,y):
    n3(x,y) = normalize( [-∂h/∂x, -∂h/∂y, 1] )

Notes
-----
- All functions support NumPy broadcasting over x,y shapes.
- Values outside the domain are still computed (smoothly) but you can
  enable hard clipping to avoid accidental usage.
- Set eps to a small positive value to avoid division by zero.
"""

@dataclass
class RBFTerrain:
    amp: float = 0.2          # amplitude (meters)
    ell: float = 0.35         # length scale (meters)
    domain: Tuple[float, float] = (-1.0, 1.0)  # square domain in x and y
    hard_clip: bool = False   # if True, clip x,y into domain before evaluating
    eps: float = 1e-8         # numerical epsilon for normalization

    # --------------------
    # Core height function
    # --------------------
    def h(self, x: ScalarOrArray, y: ScalarOrArray) -> Array:
        """Height h(x,y) of the terrain.
        Parameters
        ----------
        x, y : float or array-like (broadcastable)
        Returns
        -------
        ndarray of the broadcasted shape
        """
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        if self.hard_clip:
            lo, hi = self.domain
            x = np.clip(x, lo, hi)
            y = np.clip(y, lo, hi)
        r2 = x*x + y*y
        g = np.exp(-0.5 * r2 / (self.ell*self.ell))
        return self.amp * g

    # -----------------
    # Planar derivatives
    # -----------------
    def grad(self, x: ScalarOrArray, y: ScalarOrArray) -> Tuple[Array, Array]:
        """Gradient ∇h = [∂h/∂x, ∂h/∂y] at (x,y)."""
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        if self.hard_clip:
            lo, hi = self.domain
            x = np.clip(x, lo, hi)
            y = np.clip(y, lo, hi)
        r2 = x*x + y*y
        g = np.exp(-0.5 * r2 / (self.ell*self.ell))
        common = -(self.amp / (self.ell*self.ell)) * g  # factor in front of [x, y]
        dhdx = common * x
        dhdy = common * y
        return dhdx, dhdy

    # -----------------
    # 3D unit surface normal
    # -----------------
    def normal3d(self, x: ScalarOrArray, y: ScalarOrArray) -> Array:
        """Unit surface normal n̂(x,y) of the graph z=h(x,y).
        Returns an array with shape (..., 3), broadcasting over x,y.
        n̂ = normalize([-∂h/∂x, -∂h/∂y, 1]).
        """
        dhdx, dhdy = self.grad(x, y)
        nx = -dhdx
        ny = -dhdy
        # Construct and normalize
        # Stack will broadcast to a common shape; add last dim=3
        nx, ny = np.asarray(nx), np.asarray(ny)
        one = np.ones_like(nx)
        v = np.stack([nx, ny, one], axis=-1)
        norm = np.linalg.norm(v, axis=-1, keepdims=True)
        return v / (norm + self.eps)

    # -----------------
    # Convenience batch API
    # -----------------
    def eval_all(self, x: ScalarOrArray, y: ScalarOrArray) -> Tuple[Array, Array, Array, Array]:
        """Return (h, dhdx, dhdy, n3) in one call.
        Shapes follow NumPy broadcasting. n3 has an extra trailing axis of size 3.
        """
        z = self.h(x, y)
        dhdx, dhdy = self.grad(x, y)
        n3 = self.normal3d(x, y)
        return z, dhdx, dhdy, n3


def visualize_terrain_3d(terrain, domain=(-1, 1), res=100,
                         normal_stride=8, normal_len=0.15):
    """
    分别可视化地形高度图和法向图
    Parameters
    ----------
    terrain : object
        需要实现 h(x,y) 和 normal3d(x,y)
    domain : tuple(float,float)
        x,y 的范围
    res : int
        表面采样分辨率
    normal_stride : int
        quiver 间隔步长
    normal_len : float
        法向箭头长度
    """
    lo, hi = domain
    xs = np.linspace(lo, hi, res)
    ys = np.linspace(lo, hi, res)
    X, Y = np.meshgrid(xs, ys)
    Z = terrain.h(X, Y)

    # # 1. 高度图
    # fig1 = plt.figure(figsize=(8, 6))
    # ax1 = fig1.add_subplot(111, projection="3d")
    # ax1.plot_surface(X, Y, Z, rstride=2, cstride=2, linewidth=0, alpha=0.9)
    # ax1.set_xlabel("x")
    # ax1.set_ylabel("y")
    # ax1.set_zlabel("z = h(x, y)")
    # ax1.set_title("Terrain Surface")
    # ax1.view_init(elev=35, azim=-60)
    # plt.tight_layout()
    # plt.show()

    # # 2. 法向图
    # Xq, Yq, Zq = X[::normal_stride, ::normal_stride], \
    #              Y[::normal_stride, ::normal_stride], \
    #              Z[::normal_stride, ::normal_stride]
    # N = terrain.normal3d(Xq, Yq)

    # fig2 = plt.figure(figsize=(8, 6))
    # ax2 = fig2.add_subplot(111, projection="3d")
    # ax2.plot_surface(X, Y, Z, rstride=2, cstride=2, linewidth=0, alpha=0.3)
    # ax2.quiver(Xq, Yq, Zq, N[..., 0], N[..., 1], N[..., 2],
    #            length=normal_len, normalize=True, color="k", linewidth=0.5)
    # ax2.set_xlabel("x")
    # ax2.set_ylabel("y")
    # ax2.set_zlabel("z = h(x, y)")
    # ax2.set_title("Terrain Surface Normals")
    # ax2.view_init(elev=35, azim=-60)
    # plt.tight_layout()
    # plt.show()

    # --- surface ---
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=2, cstride=2, linewidth=0, antialiased=True)

    # --- normals (sparser) ---
    step = 8
    Xq = X[::step, ::step]; Yq = Y[::step, ::step]; Zq = Z[::step, ::step]
    Nvec = terr.normal3d(Xq, Yq)
    Nx, Ny, Nz = Nvec[...,0], Nvec[...,1], Nvec[...,2]
    ax.quiver(Xq, Yq, Zq, Nx, Ny, Nz, length=0.08, normalize=True)

    # --- equal aspect box ---
    x_range = np.ptp(xs); y_range = np.ptp(ys); z_range = np.ptp(Z)
    ax.set_box_aspect((x_range, y_range, z_range))

    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title('RBF Surface with Unit Normals (equal aspect)')
    ax.view_init(elev=35, azim=-60)
    plt.tight_layout()
    plt.show()


# ------------------------------------------------------------
# 示例运行
# ------------------------------------------------------------
if __name__ == "__main__":
    terr = RBFTerrain(amp=0.2, ell=0.35)
    visualize_terrain_3d(
        terrain=terr,
        domain=(-1, 1),
        res=150,
        normal_stride=10,
        normal_len=0.15,
    )
