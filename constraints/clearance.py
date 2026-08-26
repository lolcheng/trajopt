from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Optional, Tuple, Dict
import numpy as np

Array = np.ndarray

# Reuse the linear evaluator structures from objective.py if available.
# We duplicate small signatures here to keep this file standalone.

@dataclass
class LinMaps1D:
    N: Array      # (K, n)
    Nd: Array     # (K, n)
    Ndd: Array    # (K, n)
    T: float
    def pos(self, P: Array) -> Array:
        return self.N @ P

@dataclass
class Traj2D:
    X: LinMaps1D
    Y: LinMaps1D
    def pos(self, Px: Array, Py: Array) -> Array:
        return np.stack([self.X.pos(Px), self.Y.pos(Py)], axis=1)  # (K,2)


class ClearanceIneq:
    """Base clearance inequality:  g_k(x) = z_b(k) - h(x_b(k), y_b(k)) - d_min  >= 0

    - If you don't have an explicit base-z decision variable, set `bz_name=None` and
      the term reduces to   g_k(x) = - h(x_b(k), y_b(k)) - d_min  (i.e., z_b ≡ 0),
      which is typically NOT what you want. Prefer providing a z map/block.
    - Terrain is provided as two callables:  h(x, y)  and  grad(x, y) -> (hx, hy).
    - Jacobian includes the coupling via terrain gradient: ∂g/∂Px_b = -hx * N_x, etc.

    Bounds returned by `bounds()` are lb=0, ub=+inf (vector of length K).
    """
    def __init__(self,
                 base_maps_xyz: Tuple[LinMaps1D, LinMaps1D, Optional[LinMaps1D]],
                 base_blocks_xyz: Tuple[str, str, Optional[str]],
                 terrain_h: Callable[[Array, Array], Array],
                 terrain_grad: Optional[Callable[[Array, Array], Tuple[Array, Array]]],
                 d_min: float):
        self.X, self.Y, self.Z = base_maps_xyz
        self.bx, self.by, self.bz = base_blocks_xyz
        self.h_fn = terrain_h
        self.gh_fn = terrain_grad
        self.dmin = float(d_min)
        # cache sizes
        self.K = self.X.N.shape[0]

    # --- Constraint API ---
    def value(self, x: Array, unpack) -> Array:
        u = unpack(x)
        xb = self.X.pos(u[self.bx])
        yb = self.Y.pos(u[self.by])
        zb = self.Z.pos(u[self.bz]) if (self.Z is not None and self.bz is not None) else np.zeros_like(xb)
        h = self.h_fn(xb, yb)
        g = zb - h - self.dmin
        return g.reshape(-1)

    def jac(self, x: Array, unpack) -> Array:
        # Dense Jacobian (K x n). For large problems, consider returning sparse.
        vars_obj = unpack.__self__  # hack: access Variables from bound method
        n = vars_obj.size()
        J = np.zeros((self.K, n), dtype=float)
        u = unpack(x)
        xb = self.X.pos(u[self.bx])
        yb = self.Y.pos(u[self.by])
        # fill z block if available
        if self.Z is not None and self.bz is not None:
            slz = vars_obj.slice(self.bz)
            J[:, slz] += self.Z.N  # ∂g/∂Pz_b = N_z
        # terrain gradient for coupling terms
        if self.gh_fn is not None:
            hx, hy = self.gh_fn(xb, yb)
            hx = np.asarray(hx).reshape(-1)
            hy = np.asarray(hy).reshape(-1)
            slx = vars_obj.slice(self.bx)
            sly = vars_obj.slice(self.by)
            # ∂g/∂Px_b = - hx * N_x;  ∂g/∂Py_b = - hy * N_y
            J[:, slx] += - self.X.N * hx[:, None]
            J[:, sly] += - self.Y.N * hy[:, None]
        return J

    def bounds(self) -> Tuple[Array, Array]:
        lb = np.zeros(self.K)
        ub = np.full(self.K, np.inf)
        return lb, ub
