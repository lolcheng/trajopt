from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Tuple, Union
import numpy as np

ArrayLike = Union[float, int, Iterable[float], np.ndarray]

# ---------------------------
# VarBlock: one logical block
# ---------------------------
@dataclass
class VarBlock:
    """A logical block of decision variables with bounds and scaling.

    Parameters
    ----------
    name : str
        Unique name of the block, e.g. "Px", "Py", "slacks".
    shape : Tuple[int, ...]
        Shape of the block (e.g. (n_ctrl,) or (n_ctrl, 2)). Size = prod(shape).
    lb, ub : ArrayLike
        Lower/upper bounds. Can be scalar or array-like; will be broadcast to size.
    scale : float or ArrayLike, default 1.0
        Scaling applied to *internal* representation. Externally, values are unscaled.
        The solver sees x_scaled = x_unscaled * scale.
    init : Optional[ArrayLike]
        Initial guess (unscaled). If None, will use midpoint of bounds (or zeros if inf).
    dtype : np.dtype
        Storage dtype (float64 recommended for NLP).
    """

    name: str
    shape: Tuple[int, ...]
    lb: ArrayLike = -np.inf
    ub: ArrayLike = np.inf
    scale: ArrayLike = 1.0
    init: Optional[ArrayLike] = None
    dtype: np.dtype = np.float64

    # computed fields
    size: int = field(init=False)
    _lb: np.ndarray = field(init=False, repr=False)
    _ub: np.ndarray = field(init=False, repr=False)
    _scale: np.ndarray = field(init=False, repr=False)
    _init: Optional[np.ndarray] = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        self.size = int(np.prod(self.shape, dtype=int))
        self._lb = _as_array(self.lb, self.size, fill=-np.inf, dtype=self.dtype)
        self._ub = _as_array(self.ub, self.size, fill=np.inf, dtype=self.dtype)
        if np.any(self._lb > self._ub):
            raise ValueError(f"Block {self.name}: lb > ub at some indices.")
        self._scale = _as_array(self.scale, self.size, fill=1.0, dtype=self.dtype)
        if np.any(self._scale == 0):
            raise ValueError(f"Block {self.name}: scale must be non-zero.")
        if self.init is not None:
            # allow scalar or matching-size init by broadcasting via _as_array
            try:
                self._init = _as_array(self.init, self.size, fill=0.0, dtype=self.dtype).copy()
            except ValueError as e:
                raise ValueError(f"Block {self.name}: init size mismatch: {e}")
        else:
            # midpoint where finite, else zero
            finite = np.isfinite(self._lb) & np.isfinite(self._ub)
            mid = np.where(finite, 0.5 * (self._lb + self._ub), 0.0)
            self._init = mid.astype(self.dtype, copy=False)

    # public views (unscaled shapes)
    def lb_unscaled(self) -> np.ndarray:
        return self._lb.reshape(self.shape)

    def ub_unscaled(self) -> np.ndarray:
        return self._ub.reshape(self.shape)

    def scale_vec(self) -> np.ndarray:
        return self._scale.reshape(self.shape)

    def init_unscaled(self) -> np.ndarray:
        return self._init.reshape(self.shape)

    # scaled (solver space)
    def lb_scaled(self) -> np.ndarray:
        return (self._lb * self._scale).reshape(self.shape)

    def ub_scaled(self) -> np.ndarray:
        return (self._ub * self._scale).reshape(self.shape)

    def init_scaled(self) -> np.ndarray:
        return (self._init * self._scale).reshape(self.shape)

    def with_init(self, init: ArrayLike) -> "VarBlock":
        return VarBlock(name=self.name, shape=self.shape, lb=self._lb, ub=self._ub,
                        scale=self._scale, init=init, dtype=self.dtype)

    def freeze(self, value: Optional[ArrayLike] = None) -> "VarBlock":
        """Return a new block with lb == ub == value (or current init if None)."""
        val = self.init_unscaled() if value is None else np.asarray(value, dtype=self.dtype).reshape(self.shape)
        return VarBlock(name=self.name, shape=self.shape, lb=val, ub=val, scale=self._scale, init=val, dtype=self.dtype)


# ---------------------------------
# Variables: pack/unpack & utilities
# ---------------------------------
class Variables:
    """Manage multiple VarBlock objects and map to a single solver vector.

    The internal (solver) vector is the *scaled* concatenation of each block.
    Public APIs accept/return *unscaled* values by default.
    """

    def __init__(self, blocks: Optional[List[VarBlock]] = None) -> None:
        self.blocks: List[VarBlock] = []
        self._index: Dict[str, int] = {}
        self._offsets: np.ndarray = np.array([0], dtype=int)
        if blocks:
            for b in blocks:
                self.add_block(b)

    # ---- block management ----
    def add_block(self, block: VarBlock) -> None:
        if block.name in self._index:
            raise KeyError(f"Duplicate block name: {block.name}")
        self._index[block.name] = len(self.blocks)
        self.blocks.append(block)
        self._recompute_offsets()

    def add(self, name: str, shape: Tuple[int, ...], lb: ArrayLike = -np.inf, ub: ArrayLike = np.inf,
            scale: ArrayLike = 1.0, init: Optional[ArrayLike] = None, dtype: np.dtype = np.float64) -> None:
        self.add_block(VarBlock(name=name, shape=shape, lb=lb, ub=ub, scale=scale, init=init, dtype=dtype))

    def _recompute_offsets(self) -> None:
        sizes = [b.size for b in self.blocks]
        self._offsets = np.concatenate([[0], np.cumsum(sizes)])

    # ---- slicing & views ----
    def slice(self, name: str) -> slice:
        i = self._index[name]
        s = int(self._offsets[i]); e = int(self._offsets[i+1])
        return slice(s, e)

    def shape(self, name: str) -> Tuple[int, ...]:
        return self.blocks[self._index[name]].shape

    def size(self) -> int:
        return int(self._offsets[-1])

    def names(self) -> List[str]:
        return [b.name for b in self.blocks]

    # ---- bounds & scaling (solver space) ----
    def bounds_scaled(self) -> Tuple[np.ndarray, np.ndarray]:
        lbs = [b.lb_scaled().reshape(-1) for b in self.blocks]
        ubs = [b.ub_scaled().reshape(-1) for b in self.blocks]
        return np.concatenate(lbs), np.concatenate(ubs)

    def scale_vector(self) -> np.ndarray:
        return np.concatenate([b.scale_vec().reshape(-1) for b in self.blocks])

    # ---- packing / unpacking ----
    def pack(self, values: Mapping[str, ArrayLike]) -> np.ndarray:
        """Pack *unscaled* block values into a single *scaled* solver vector.
        Missing blocks will use their block.init_unscaled().
        """
        out = np.zeros(self.size(), dtype=self.blocks[0].dtype if self.blocks else np.float64)
        for b in self.blocks:
            v = values.get(b.name, b.init_unscaled())
            v = np.asarray(v, dtype=b.dtype).reshape(b.shape)
            sl = self.slice(b.name)
            out[sl] = (v * b.scale_vec()).reshape(-1)
        return out

    def unpack(self, x_scaled: np.ndarray) -> Dict[str, np.ndarray]:
        """Unpack a *scaled* solver vector to a dict of *unscaled* block arrays with shapes."""
        x_scaled = np.asarray(x_scaled)
        out: Dict[str, np.ndarray] = {}
        for b in self.blocks:
            sl = self.slice(b.name)
            arr = (x_scaled[sl].reshape(-1) / b.scale_vec().reshape(-1)).astype(b.dtype, copy=False)
            out[b.name] = arr.reshape(b.shape)
        return out

    # ---- initial guess ----
    def initial_x(self, overrides: Optional[Mapping[str, ArrayLike]] = None) -> np.ndarray:
        overrides = overrides or {}
        return self.pack(overrides)

    # ---- utilities ----
    def clip_inplace(self, x_scaled: np.ndarray) -> None:
        lb, ub = self.bounds_scaled()
        np.clip(x_scaled, lb, ub, out=x_scaled)

    def is_feasible(self, x_scaled: np.ndarray, tol: float = 0.0) -> bool:
        lb, ub = self.bounds_scaled()
        return bool(np.all(x_scaled >= lb - tol) and np.all(x_scaled <= ub + tol))

    def summary(self) -> str:
        lines = ["Variables summary:"]
        off = 0
        for b in self.blocks:
            lines.append(
                f" - {b.name:12s} shape={b.shape} size={b.size:4d} offset=[{off}:{off+b.size}] scale=min{b._scale.min():.3g}/max{b._scale.max():.3g}"
            )
            off += b.size
        lines.append(f"Total size: {self.size()}")
        return "\n".join(lines)


# --------------
# helper routines
# --------------

def _as_array(x: ArrayLike, size: int, fill: float, dtype: np.dtype) -> np.ndarray:
    if isinstance(x, (int, float)):
        return np.full(size, x, dtype=dtype)
    arr = np.asarray(x, dtype=dtype).reshape(-1)
    if arr.size == 1:
        return np.full(size, float(arr[0]), dtype=dtype)
    if arr.size == size:
        return arr
    raise ValueError(f"Cannot broadcast array of size {arr.size} to size {size}.")


# ---------------------------
# quick self-test (optional)
# ---------------------------
if __name__ == "__main__":
    n_ctrl = 5
    vx = VarBlock(name="Px", shape=(n_ctrl,), lb=-10, ub=10, scale=1.0, init=np.linspace(0, 1, n_ctrl))
    vy = VarBlock(name="Py", shape=(n_ctrl,), lb=-5, ub=5, scale=2.0, init=0.0)
    vs = Variables([vx, vy])
    print(vs.summary())
    x0 = vs.initial_x()
    print("x0 (scaled):", x0)
    parts = vs.unpack(x0)
    print("Px unscaled:", parts["Px"]) 
    print("Py unscaled:", parts["Py"]) 
    assert vs.is_feasible(x0)
