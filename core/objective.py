from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Callable
import numpy as np

# NOTE
# -----
# This module implements objective terms that match the user's PDF (sections 2.1–2.5),
# but in a numerically friendly, differentiable, and solver-agnostic form.
# It works with the Variables class (variables.py) by consuming the packed solver
# vector x (scaled) and an "unpack(x) -> dict" callback to access per-block arrays.
#
# Design goals:
#  - Pure NumPy (no CasADi), so it can plug into SciPy/NLopt/cyipopt or your own SQP.
#  - Closed-form gradients for all quadratic terms (fast + exact).
#  - Hinge/angle terms are smooth or piecewise-smooth with clear gradients.
#  - All per-sample geometry (terrain height/normal, lateral vectors, quadrature weights)
#    are passed in precomputed arrays to keep this module independent of geometry code.
#
# Nomenclature and mapping to PDF:
#  - Terminal goal (2.1):    ω1 * ( ||p_b(s_N) - p_goal||^2 + angErr(ψ_b(s_N), ψ_goal) )
#  - Accel smoothness (2.2): ω2 * ∫ ||p_b''(s)|| ds  +  ω3 * ∫ (||p_L''(s)|| + ||p_R''(s)||) ds
#    → Implemented as weighted sum of squared accelerations with quadrature weights (L2^2).
#  - Lateral slip (2.3):     ω4 * Σ ( (v_L · b_L)^2 + (v_R · b_R)^2 )  using lateral unit vectors b.
#  - Track width (2.4):      ω5 * Σ ( ((p_L - p_R) · b_avg) - a1 )^2
#  - Base height (2.5):      ω6 * ∫ [ z_b - h(x_b, y_b) - d_min ]_+^2  (hinge penalty)
#
# Block names (default, but configurable):
#  - Base:   Px_b, Py_b, (optional) Pz_b, (optional) Psi_b
#  - Left:   Px_L, Py_L
#  - Right:  Px_R, Py_R
#
# Linear evaluation via pre-built matrices:
#   Xb = Nb @ Px_b,   Yb = Nb @ Py_b,   Xdd_b = Ndd_b @ Px_b / T^2, ...
#   and similarly for wheels L/R. (Nd: first derivative, Ndd: second derivative)

Array = np.ndarray

# ----------------------------
# Utility / small linear algebra
# ----------------------------

def wrap_to_pi(angle: Array) -> Array:
    """Wrap angle to (-pi, pi]."""
    a = (angle + np.pi) % (2 * np.pi) - np.pi
    # Keep +pi instead of -pi for continuity near boundary if needed
    return a

# -----------------
# Objective Protocol
# -----------------

class ObjectiveTerm:
    """Interface for an objective term: value and gradient.
    Subclasses must implement value() and grad().
    """
    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        raise NotImplementedError

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        raise NotImplementedError


class WeightedObjective:
    """Weighted sum of objective terms."""
    def __init__(self, terms: Dict[str, Tuple[ObjectiveTerm, float]]):
        self.terms = terms

    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        return float(sum(w * t.value(x, unpack) for t, w in self.terms.values()))

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        g = None
        for t, w in self.terms.values():
            gi = w * t.grad(x, unpack)
            g = gi if g is None else g + gi
        return g if g is not None else np.zeros_like(x)


# ----------------------------
# Linear evaluation containers
# ----------------------------

@dataclass
class LinMaps1D:
    """Linear maps for one coordinate of a trajectory: value, velocity, acceleration.
    All matrices are (K, n_ctrl) for K evaluation points. Output uses these maps as:
        pos = N @ P
        vel = Nd @ P / T
        acc = Ndd @ P / T^2
    """
    N: Array      # (K, n)
    Nd: Array     # (K, n)
    Ndd: Array    # (K, n)
    T: float

    def pos(self, P: Array) -> Array:
        return self.N @ P

    def vel(self, P: Array) -> Array:
        return (self.Nd @ P) / self.T

    def acc(self, P: Array) -> Array:
        return (self.Ndd @ P) / (self.T * self.T)


@dataclass
class Traj2D:
    """2D trajectory maps for x/y coordinates, using separate LinMaps1D for each axis."""
    X: LinMaps1D
    Y: LinMaps1D

    def pos(self, Px: Array, Py: Array) -> Array:
        return np.stack([self.X.pos(Px), self.Y.pos(Py)], axis=1)  # (K,2)

    def vel(self, Px: Array, Py: Array) -> Array:
        return np.stack([self.X.vel(Px), self.Y.vel(Py)], axis=1)

    def acc(self, Px: Array, Py: Array) -> Array:
        return np.stack([self.X.acc(Px), self.Y.acc(Py)], axis=1)


# -----------------
# Objective Terms
# -----------------

class TerminalGoal(ObjectiveTerm):
    """ω1 * ( ||p_b(s_N) - p_goal||^2 + yaw_term )

    Parameters
    ----------
    base_maps : Traj2D for base
    idx_N     : int, index of terminal sample
    p_goal    : (2,) array for goal position [x_g, y_g]
    yaw_row   : Optional[(1,n)] row vector to evaluate ψ_b(s_N) = yaw_row @ Ppsi
    psi_block : Optional[str] name of psi block in Variables
    yaw_goal  : Optional[float] goal yaw angle
    yaw_weight: float weight multiplier for yaw term inside this term (default 1.0)
    """
    def __init__(self,
                 base_maps: Traj2D,
                 idx_N: int,
                 p_goal: Array,
                 yaw_row: Optional[Array] = None,
                 psi_block: Optional[str] = None,
                 yaw_goal: Optional[float] = None,
                 yaw_weight: float = 1.0):
        self.base_maps = base_maps
        self.idx_N = int(idx_N)
        self.p_goal = np.asarray(p_goal, dtype=float).reshape(2)
        self.yaw_row = None if yaw_row is None else np.asarray(yaw_row, dtype=float).reshape(1, -1)
        self.psi_block = psi_block
        self.yaw_goal = float(yaw_goal) if yaw_goal is not None else None
        self.yaw_weight = float(yaw_weight)

    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        u = unpack(x)
        pbN = np.array([
            (self.base_maps.X.N[self.idx_N] @ u['Px_b']),
            (self.base_maps.Y.N[self.idx_N] @ u['Py_b'])
        ])
        val = float(np.dot(pbN - self.p_goal, pbN - self.p_goal))
        # yaw term if provided
        if self.yaw_row is not None and self.psi_block is not None and self.yaw_goal is not None:
            psiN = float(self.yaw_row @ u[self.psi_block])
            dpsi = float(wrap_to_pi(psiN - self.yaw_goal))
            val += self.yaw_weight * (dpsi * dpsi)
        return val

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        u = unpack(x)
        # gradients are formed in solver space, so we need a full-sized grad vector
        g = np.zeros_like(x)
        # base pos part
        Nx = self.base_maps.X.N[self.idx_N]  # (n_x,)
        Ny = self.base_maps.Y.N[self.idx_N]
        pbN = np.array([Nx @ u['Px_b'], Ny @ u['Py_b']])
        diff = pbN - self.p_goal  # (2,)
        # d / d Px_b : 2 * diff_x * Nx
        # d / d Py_b : 2 * diff_y * Ny
        slx = slice(*unpack.__self__.slice('Px_b').indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        sly = slice(*unpack.__self__.slice('Py_b').indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        g[slx] += 2.0 * diff[0] * Nx
        g[sly] += 2.0 * diff[1] * Ny
        # yaw part
        if self.yaw_row is not None and self.psi_block is not None and self.yaw_goal is not None:
            psiN = float(self.yaw_row @ u[self.psi_block])
            dpsi = float(wrap_to_pi(psiN - self.yaw_goal))
            dpsi_dP = self.yaw_row.reshape(-1)  # (n_psi,)
            slp = slice(*unpack.__self__.slice(self.psi_block).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
            g[slp] += 2.0 * self.yaw_weight * dpsi * dpsi_dP
        return g


class AccelSmoothness(ObjectiveTerm):
    """Acceleration L2^2 with quadrature weights.
    For each trajectory (base and/or wheels), penalize sum_k wq[k] * ||a(k)||^2.

    Pass in maps for x/y of the trajectories you want to penalize.
    """
    def __init__(self,
                 maps_and_blocks: Tuple[Tuple[Traj2D, str, str], ...],
                 quad_w: Array):
        self.items = maps_and_blocks  # tuples of (Traj2D, name_x, name_y)
        self.quad_w = quad_w.reshape(-1, 1)  # (K,1)

    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        u = unpack(x)
        val = 0.0
        for maps, name_x, name_y in self.items:
            ax = maps.X.acc(u[name_x])  # (K,)
            ay = maps.Y.acc(u[name_y])
            val += float((self.quad_w[:, 0] * (ax * ax + ay * ay)).sum())
        return val

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        u = unpack(x)
        g = np.zeros_like(x)
        for maps, name_x, name_y in self.items:
            ax = maps.X.acc(u[name_x])  # (K,)
            ay = maps.Y.acc(u[name_y])
            # d/dPx: 2 * sum_k wq * ax * (Ndd_x[k,:]/T^2)
            Jx = maps.X.Ndd / (maps.X.T * maps.X.T)  # (K,n)
            Jy = maps.Y.Ndd / (maps.Y.T * maps.Y.T)
            gx = 2.0 * (Jx * self.quad_w) .T @ ax  # (n,)
            gy = 2.0 * (Jy * self.quad_w) .T @ ay
            slx = slice(*unpack.__self__.slice(name_x).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
            sly = slice(*unpack.__self__.slice(name_y).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
            g[slx] += gx.reshape(-1)
            g[sly] += gy.reshape(-1)
        return g


class LateralSlip(ObjectiveTerm):
    """ω4 * Σ ( (v_L · b_L)^2 + (v_R · b_R)^2 )

    b_L, b_R: (K,2) lateral unit vectors per sample (precomputed)
    We ignore d(b)/d(x) for simplicity (Gauss-Newton style). If you need full
    derivatives, pass b as a function of x and implement its jacobian.
    """
    def __init__(self,
                 left_maps: Traj2D, right_maps: Traj2D,
                 left_blocks: Tuple[str, str], right_blocks: Tuple[str, str],
                 b_left: Array, b_right: Array, quad_w: Optional[Array] = None):
        self.L = left_maps
        self.R = right_maps
        self.L_blocks = left_blocks
        self.R_blocks = right_blocks
        self.bL = b_left  # (K,2)
        self.bR = b_right # (K,2)
        self.quad_w = None if quad_w is None else quad_w.reshape(-1, 1)

    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        u = unpack(x)
        vL = self.L.vel(u[self.L_blocks[0]], u[self.L_blocks[1]])  # (K,2)
        vR = self.R.vel(u[self.R_blocks[0]], u[self.R_blocks[1]])
        sL = np.sum(vL * self.bL, axis=1)  # (K,)
        sR = np.sum(vR * self.bR, axis=1)
        if self.quad_w is None:
            return float((sL * sL + sR * sR).sum())
        else:
            return float((self.quad_w[:, 0] * (sL * sL + sR * sR)).sum())

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        u = unpack(x)
        g = np.zeros_like(x)
        # Left
        vLx = (self.L.X.Nd @ u[self.L_blocks[0]]) / self.L.X.T  # (K,)
        vLy = (self.L.Y.Nd @ u[self.L_blocks[1]]) / self.L.Y.T
        sL = vLx * self.bL[:, 0] + vLy * self.bL[:, 1]          # (K,)
        w = 1.0 if self.quad_w is None else self.quad_w[:, 0]
        # d/dPx_L: 2 * sum_k w * sL_k * (Nd_x[k,:]/T) * bL_x[k]
        JxL = (self.L.X.Nd / self.L.X.T) * self.bL[:, [0]]      # (K,n)
        JyL = (self.L.Y.Nd / self.L.Y.T) * self.bL[:, [1]]
        gxL = 2.0 * (JxL.T @ (w * sL))
        gyL = 2.0 * (JyL.T @ (w * sL))
        slx = slice(*unpack.__self__.slice(self.L_blocks[0]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        sly = slice(*unpack.__self__.slice(self.L_blocks[1]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        g[slx] += gxL.reshape(-1)
        g[sly] += gyL.reshape(-1)
        # Right
        vRx = (self.R.X.Nd @ u[self.R_blocks[0]]) / self.R.X.T
        vRy = (self.R.Y.Nd @ u[self.R_blocks[1]]) / self.R.Y.T
        sR = vRx * self.bR[:, 0] + vRy * self.bR[:, 1]
        JxR = (self.R.X.Nd / self.R.X.T) * self.bR[:, [0]]
        JyR = (self.R.Y.Nd / self.R.Y.T) * self.bR[:, [1]]
        gxR = 2.0 * (JxR.T @ (w * sR))
        gyR = 2.0 * (JyR.T @ (w * sR))
        srx = slice(*unpack.__self__.slice(self.R_blocks[0]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        sry = slice(*unpack.__self__.slice(self.R_blocks[1]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        g[srx] += gxR.reshape(-1)
        g[sry] += gyR.reshape(-1)
        return g


class TrackWidthSoft(ObjectiveTerm):
    """ω5 * Σ ( ((p_L - p_R) · b_avg) - a1 )^2

    b_avg: (K,2) average lateral vectors per sample (unit or not). Ignoring d(b_avg)/d(x).
    """
    def __init__(self,
                 left_maps: Traj2D, right_maps: Traj2D,
                 left_blocks: Tuple[str, str], right_blocks: Tuple[str, str],
                 b_avg: Array, a1: float, quad_w: Optional[Array] = None):
        self.L = left_maps
        self.R = right_maps
        self.L_blocks = left_blocks
        self.R_blocks = right_blocks
        self.bavg = b_avg  # (K,2)
        self.a1 = float(a1)
        self.quad_w = None if quad_w is None else quad_w.reshape(-1, 1)

    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        u = unpack(x)
        pL = self.L.pos(u[self.L_blocks[0]], u[self.L_blocks[1]])  # (K,2)
        pR = self.R.pos(u[self.R_blocks[0]], u[self.R_blocks[1]])
        proj = np.sum((pL - pR) * self.bavg, axis=1)  # (K,)
        r = proj - self.a1
        if self.quad_w is None:
            return float((r * r).sum())
        else:
            return float((self.quad_w[:, 0] * (r * r)).sum())

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        u = unpack(x)
        g = np.zeros_like(x)
        # residual r = ((pL - pR)·bavg) - a1
        # ∂r/∂Px_L = bavg_x * N_x; ∂r/∂Py_L = bavg_y * N_y;  ∂r/∂Px_R = -bavg_x * N_x; ...
        w = 1.0 if self.quad_w is None else self.quad_w[:, 0]
        NxL = self.L.X.N  # (K,n)
        NyL = self.L.Y.N
        NxR = self.R.X.N
        NyR = self.R.Y.N
        b_x = self.bavg[:, 0:1]  # (K,1)
        b_y = self.bavg[:, 1:2]
        # Build Jacobian-vector products efficiently
        # gxL = 2 * Σ w * r_k * b_x[k] * NxL[k,:]
        # etc.
        pL = self.L.pos(u[self.L_blocks[0]], u[self.L_blocks[1]])
        pR = self.R.pos(u[self.R_blocks[0]], u[self.R_blocks[1]])
        r = np.sum((pL - pR) * self.bavg, axis=1)
        r = w * (r - self.a1)  # (K,)
        gxL = 2.0 * (NxL.T @ (r * b_x[:, 0]))
        gyL = 2.0 * (NyL.T @ (r * b_y[:, 0]))
        gxR = -2.0 * (NxR.T @ (r * b_x[:, 0]))
        gyR = -2.0 * (NyR.T @ (r * b_y[:, 0]))
        slxL = slice(*unpack.__self__.slice(self.L_blocks[0]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        slyL = slice(*unpack.__self__.slice(self.L_blocks[1]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        slxR = slice(*unpack.__self__.slice(self.R_blocks[0]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        slyR = slice(*unpack.__self__.slice(self.R_blocks[1]).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
        g[slxL] += gxL.reshape(-1)
        g[slyL] += gyL.reshape(-1)
        g[slxR] += gxR.reshape(-1)
        g[slyR] += gyR.reshape(-1)
        return g


class BaseClearanceHinge(ObjectiveTerm):
    """ω6 * ∫ [ z_b - h(x_b, y_b) - d_min ]_+^2  with quadrature weights.

    We require linear maps for base x/y/z and arrays for terrain height h(x,y) and its gradient
    wrt (x,y) at each sample (for gradient propagation). If you don't have terrain gradient,
    you may set grad_h=None to use only ∂/∂Pz and ignore cross-terms.
    """
    def __init__(self,
                 base_maps_xyz: Tuple[LinMaps1D, LinMaps1D, Optional[LinMaps1D]],
                 base_blocks_xyz: Tuple[str, str, Optional[str]],
                 h_vals: Array,              # (K,) terrain height at base (x,y)
                 grad_h: Optional[Array],    # (K,2) terrain gradient wrt (x,y)
                 d_min: float,
                 quad_w: Array,
                 smooth_hinge_eps: float = 0.0):
        self.X, self.Y, self.Z = base_maps_xyz
        self.bx, self.by, self.bz = base_blocks_xyz
        self.h_vals = h_vals.reshape(-1)
        self.grad_h = None if grad_h is None else grad_h.reshape(-1, 2)
        self.d_min = float(d_min)
        self.quad_w = quad_w.reshape(-1)
        self.eps = float(smooth_hinge_eps)

    def _hinge(self, g: Array) -> Tuple[Array, Array]:
        """Return (phi, dphi/dg) for hinge phi(g) = [g]_+^2.
        If eps>0, use smooth approximation: phi = (softplus(g/eps)*eps)^2.
        """
        if self.eps <= 0.0:
            mask = (g > 0).astype(float)
            phi = (g * mask) ** 2
            dphi = 2.0 * g * mask
            return phi, dphi
        else:
            # softplus(s) = log(1+exp(s)); phi = (eps*softplus(g/eps))^2
            s = g / self.eps
            sp = np.log1p(np.exp(s))
            phi = (self.eps * sp) ** 2
            # d/dg phi = 2 * eps^2 * softplus(s) * sigmoid(s) * (1/eps)
            sig = 1.0 / (1.0 + np.exp(-s))
            dphi = 2.0 * (self.eps * sp) * sig
            return phi, dphi

    def value(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> float:
        u = unpack(x)
        xb = self.X.pos(u[self.bx])
        yb = self.Y.pos(u[self.by])
        zb = self.Z.pos(u[self.bz]) if (self.Z is not None and self.bz is not None) else np.zeros_like(xb)
        g_clear = zb - self.h_vals - self.d_min  # (K,)
        phi, _ = self._hinge(g_clear)
        return float((self.quad_w * phi).sum())

    def grad(self, x: Array, unpack: Callable[[Array], Dict[str, Array]]) -> Array:
        u = unpack(x)
        xb = self.X.pos(u[self.bx])
        yb = self.Y.pos(u[self.by])
        zb = self.Z.pos(u[self.bz]) if (self.Z is not None and self.bz is not None) else np.zeros_like(xb)
        g_clear = zb - self.h_vals - self.d_min
        _, dphi = self._hinge(g_clear)  # (K,)
        w = self.quad_w * dphi  # (K,)
        g = np.zeros_like(x)
        # ∂phi/∂Pz_b = w * ∂g/∂Pz_b = w * N_z
        if self.Z is not None and self.bz is not None:
            Nz = self.Z.N
            gz = (Nz.T @ w)
            slz = slice(*unpack.__self__.slice(self.bz).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
            g[slz] += gz.reshape(-1)
        # terrain coupling: ∂g/∂Px_b = - grad_h_x * N_x, ∂g/∂Py_b = - grad_h_y * N_y
        if self.grad_h is not None:
            Nx = self.X.N
            Ny = self.Y.N
            gx = - (Nx.T @ (w * self.grad_h[:, 0]))
            gy = - (Ny.T @ (w * self.grad_h[:, 1]))
            slx = slice(*unpack.__self__.slice(self.bx).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
            sly = slice(*unpack.__self__.slice(self.by).indices(unpack.__self__.size()))  # type: ignore[attr-defined]
            g[slx] += gx.reshape(-1)
            g[sly] += gy.reshape(-1)
        return g


# -------------------------
# High-level builder helper
# -------------------------

class ObjectiveBuilder:
    """Convenience helper to assemble a WeightedObjective from provided geometry/maps.
    You can add terms one by one and finally call build().
    """
    def __init__(self):
        self.terms: Dict[str, Tuple[ObjectiveTerm, float]] = {}

    def add(self, name: str, term: ObjectiveTerm, weight: float) -> 'ObjectiveBuilder':
        self.terms[name] = (term, float(weight))
        return self

    def build(self) -> WeightedObjective:
        return WeightedObjective(self.terms)
