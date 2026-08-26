from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple
import numpy as np

Array = np.ndarray

@dataclass
class LinMaps1D:
    N: Array
    Nd: Array
    Ndd: Array
    T: float
    def pos(self, P: Array) -> Array:
        return self.N @ P

@dataclass
class Traj2D:
    X: LinMaps1D
    Y: LinMaps1D
    def pos(self, Px: Array, Py: Array) -> Array:
        return np.stack([self.X.pos(Px), self.Y.pos(Py)], axis=1)


class TrackWidthEquality:
    """Track width soft equality as constraint (vector residual with tight bounds).

    Residual per sample:
        r_k(x) = ((p_L - p_R) · b_avg)_k - a1
    Bounds:
        r in [-tol, +tol]  (elementwise)

    If you want an exact equality, set tol=0.
    """
    def __init__(self,
                 left_maps: Traj2D, right_maps: Traj2D,
                 left_blocks: Tuple[str, str], right_blocks: Tuple[str, str],
                 b_avg: Array, a1: float, tol: float = 0.0):
        self.L = left_maps
        self.R = right_maps
        self.Lb = left_blocks
        self.Rb = right_blocks
        self.bavg = np.asarray(b_avg, dtype=float)
        self.a1 = float(a1)
        self.tol = float(tol)
        self.K = self.bavg.shape[0]

    def value(self, x: Array, unpack) -> Array:
        u = unpack(x)
        pL = self.L.pos(u[self.Lb[0]], u[self.Lb[1]])  # (K,2)
        pR = self.R.pos(u[self.Rb[0]], u[self.Rb[1]])
        proj = np.sum((pL - pR) * self.bavg, axis=1)
        r = proj - self.a1
        return r.reshape(-1)

    def jac(self, x: Array, unpack) -> Array:
        vars_obj = unpack.__self__
        n = vars_obj.size()
        J = np.zeros((self.K, n), dtype=float)
        b_x = self.bavg[:, 0:1]
        b_y = self.bavg[:, 1:2]
        # Left contributions:  +b_x * N_x, +b_y * N_y
        slxL = vars_obj.slice(self.Lb[0]); slyL = vars_obj.slice(self.Lb[1])
        J[:, slxL] += self.L.X.N * b_x
        J[:, slyL] += self.L.Y.N * b_y
        # Right contributions: -b_x * N_x, -b_y * N_y
        slxR = vars_obj.slice(self.Rb[0]); slyR = vars_obj.slice(self.Rb[1])
        J[:, slxR] += - self.R.X.N * b_x
        J[:, slyR] += - self.R.Y.N * b_y
        return J

    def bounds(self):
        lb = np.full(self.K, -self.tol)
        ub = np.full(self.K,  self.tol)
        return lb, ub


class TrackNoCrossIneq:
    """Non-crossing inequality:
        g_k(x) = ((p_L - p_R) · b_avg)_k  >=  a_min
    Bounds:  g in [a_min, +inf)

    Use this together with TrackWidthEquality if you want both
    a nominal track width and a strict non-crossing direction.
    """
    def __init__(self,
                 left_maps: Traj2D, right_maps: Traj2D,
                 left_blocks: Tuple[str, str], right_blocks: Tuple[str, str],
                 b_avg: Array, a_min: float):
        self.L = left_maps
        self.R = right_maps
        self.Lb = left_blocks
        self.Rb = right_blocks
        self.bavg = np.asarray(b_avg, dtype=float)
        self.amin = float(a_min)
        self.K = self.bavg.shape[0]

    def value(self, x: Array, unpack) -> Array:
        u = unpack(x)
        pL = self.L.pos(u[self.Lb[0]], u[self.Lb[1]])
        pR = self.R.pos(u[self.Rb[0]], u[self.Rb[1]])
        g = np.sum((pL - pR) * self.bavg, axis=1)  # (K,)
        return g.reshape(-1)

    def jac(self, x: Array, unpack) -> Array:
        vars_obj = unpack.__self__
        n = vars_obj.size()
        J = np.zeros((self.K, n), dtype=float)
        b_x = self.bavg[:, 0:1]
        b_y = self.bavg[:, 1:2]
        # ∂/∂Px_L = +b_x * N_x,  ∂/∂Py_L = +b_y * N_y
        slxL = vars_obj.slice(self.Lb[0]); slyL = vars_obj.slice(self.Lb[1])
        J[:, slxL] += self.L.X.N * b_x
        J[:, slyL] += self.L.Y.N * b_y
        # ∂/∂Px_R = -b_x * N_x,  ∂/∂Py_R = -b_y * N_y
        slxR = vars_obj.slice(self.Rb[0]); slyR = vars_obj.slice(self.Rb[1])
        J[:, slxR] += - self.R.X.N * b_x
        J[:, slyR] += - self.R.Y.N * b_y
        return J

    def bounds(self):
        lb = np.full(self.K, self.amin)
        ub = np.full(self.K, np.inf)
        return lb, ub
