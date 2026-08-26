from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, Dict, Tuple, Optional
import numpy as np

Array = np.ndarray

# ------------------------------------------------------------
# Problem / Solver base interfaces
# ------------------------------------------------------------

class VariablesProtocol(Protocol):
    def unpack(self, x_scaled: Array) -> Dict[str, Array]: ...
    def bounds_scaled(self) -> Tuple[Array, Array]: ...
    def size(self) -> int: ...

class NLPProblem(Protocol):
    """Minimal protocol that a solver expects.

    Implementations typically wrap your Variables, Objective, and Constraints
    into the following methods. The functions must work on the *scaled* solver
    vector x (the one you pass to a solver). If you use trajopt.core.Variables,
    you can forward-pack/unpack internally.
    """
    variables: VariablesProtocol

    def eval_f(self, x: Array) -> float: ...                # scalar objective
    def eval_grad_f(self, x: Array) -> Array: ...            # (n,)
    def eval_g(self, x: Array) -> Array: ...                 # (m,)
    def eval_jac_g(self, x: Array) -> Array: ...             # (m,n) or sparse
    def bounds(self) -> Tuple[Array, Array, Array, Array]: ...  # (xlb, xub, glb, gub)

@dataclass
class SolveResult:
    x: Array                     # optimal x (scaled vector)
    obj: float                   # objective value at x
    status: str                  # solver-dependent status string
    success: bool                # True if solver reports success
    n_iter: Optional[int] = None
    info: Optional[dict] = None  # raw stats from the backend


class Solver(Protocol):
    def solve(self, prob: NLPProblem, x0: Array, options: Optional[Dict] = None) -> SolveResult: ...


# ------------------------------------------------------------
# Utility: check shapes / coercions
# ------------------------------------------------------------

def as_1d(a: Array, name: str) -> Array:
    a = np.asarray(a, dtype=float).reshape(-1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be 1-D")
    return a


def check_problem_shapes(prob: NLPProblem) -> None:
    n = prob.variables.size()
    xlb, xub, glb, gub = prob.bounds()
    xlb = as_1d(xlb, 'xlb'); xub = as_1d(xub, 'xub')
    glb = as_1d(glb, 'glb'); gub = as_1d(gub, 'gub')
    if xlb.size != n or xub.size != n:
        raise ValueError("x bounds size mismatch with number of variables")
    # constraint vector size is runtime-defined; we don't know m here.

