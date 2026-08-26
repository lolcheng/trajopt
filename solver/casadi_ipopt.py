from __future__ import annotations
from typing import Dict, Optional
import numpy as np

try:
    import casadi as ca
except Exception as e:  # pragma: no cover
    raise ImportError("casadi is required for casadi_ipopt solver adapter") from e

from .base import NLPProblem, SolveResult, check_problem_shapes

Array = np.ndarray

# ------------------------------------------------------------------
# CasADi Callback wrappers to expose numpy-based f/∇f and g/J
# ------------------------------------------------------------------

class _ObjCallback(ca.Callback):
    def __init__(self, name: str, prob: NLPProblem, opts: Optional[Dict] = None):
        ca.Callback.__init__(self)
        self.prob = prob
        self.construct(name, {'enable_fd': False, **(opts or {})})

    def get_n_in(self):  # x
        return 1

    def get_n_out(self):  # f
        return 1

    def get_sparsity_in(self, i):
        # x is dense vector of length n
        n = self.prob.variables.size()
        return ca.Sparsity.dense(n, 1)

    def get_sparsity_out(self, i):
        # f is scalar
        return ca.Sparsity.dense(1, 1)

    def eval(self, arg):
        x = np.array(arg[0]).reshape(-1)
        f = float(self.prob.eval_f(x))
        return [ca.DM([f])]

    # Provide analytic gradient to avoid AD through a black-box
    def has_jacobian(self):
        return True

    def get_jacobian(self, *args):
        return _ObjGradCallback('grad_f', self.prob)


class _ObjGradCallback(ca.Callback):
    def __init__(self, name: str, prob: NLPProblem):
        ca.Callback.__init__(self)
        self.prob = prob
        self.construct(name, {'enable_fd': False})

    def get_n_in(self):  # x
        return 1

    def get_n_out(self):  # df/dx (n x 1)
        return 1

    def get_sparsity_in(self, i):
        n = self.prob.variables.size()
        return ca.Sparsity.dense(n, 1)

    def get_sparsity_out(self, i):
        n = self.prob.variables.size()
        return ca.Sparsity.dense(n, 1)

    def eval(self, arg):
        x = np.array(arg[0]).reshape(-1)
        g = np.asarray(self.prob.eval_grad_f(x), dtype=float).reshape(-1, 1)
        return [ca.DM(g)]


class _ConsCallback(ca.Callback):
    def __init__(self, name: str, prob: NLPProblem, m: int, opts: Optional[Dict] = None):
        ca.Callback.__init__(self)
        self.prob = prob
        self.m = int(m)
        self.construct(name, {'enable_fd': False, **(opts or {})})

    def get_n_in(self):  # x
        return 1

    def get_n_out(self):  # g(x) in R^m
        return 1

    def get_sparsity_in(self, i):
        n = self.prob.variables.size()
        return ca.Sparsity.dense(n, 1)

    def get_sparsity_out(self, i):
        return ca.Sparsity.dense(self.m, 1)

    def eval(self, arg):
        x = np.array(arg[0]).reshape(-1)
        g = np.asarray(self.prob.eval_g(x), dtype=float).reshape(self.m, 1)
        return [ca.DM(g)]

    def has_jacobian(self):
        return True

    def get_jacobian(self, *args):
        return _ConsJacCallback('jac_g', self.prob, self.m)


class _ConsJacCallback(ca.Callback):
    def __init__(self, name: str, prob: NLPProblem, m: int):
        ca.Callback.__init__(self)
        self.prob = prob
        self.m = int(m)
        self.construct(name, {'enable_fd': False})

    def get_n_in(self):  # x
        return 1

    def get_n_out(self):  # J(x) in R^{m x n} (as dense column-stacked vector)
        return 1

    def get_sparsity_in(self, i):
        n = self.prob.variables.size()
        return ca.Sparsity.dense(n, 1)

    def get_sparsity_out(self, i):
        # CasADi expects column-stacked jac as (m*n, 1); but when connected as a
        # Function returning an (m x n) DM, it will be vectorized internally.
        m, n = self.m, self.prob.variables.size()
        return ca.Sparsity.dense(m, n)

    def eval(self, arg):
        x = np.array(arg[0]).reshape(-1)
        J = np.asarray(self.prob.eval_jac_g(x), dtype=float)
        if J.ndim != 2:
            J = J.toarray()  # allow scipy.sparse
        m, n = self.m, self.prob.variables.size()
        J = J.reshape(m, n)
        return [ca.DM(J)]


# ------------------------------------------------------------------
# Public solver
# ------------------------------------------------------------------

class CasadiIpoptSolver:
    """Adapter: use CasADi's nlpsol('ipopt', ...) to solve a problem defined
    by numpy-callable objective/gradient and constraints/jacobian.

    This lets you keep a single model (NumPy) and still call IPOPT via CasADi.
    """
    def __init__(self, ipopt_options: Optional[Dict] = None, casadi_options: Optional[Dict] = None):
        self.ipopt_options = ipopt_options or {
            'ipopt.print_level': 5,
            'print_time': 0,
            'ipopt.mu_strategy': 'adaptive',
            'ipopt.hessian_approximation': 'limited-memory',
        }
        self.casadi_options = casadi_options or {}

    def solve(self, prob: NLPProblem, x0: Array, options: Optional[Dict] = None) -> SolveResult:
        check_problem_shapes(prob)
        n = prob.variables.size()
        xlb, xub, glb, gub = prob.bounds()
        x_sym = ca.MX.sym('x', n)

        # Build callbacks
        # constraint dimension m from bounds
        m = int(np.asarray(glb).size)
        f_cb = _ObjCallback('f_cb', prob)
        g_cb = _ConsCallback('g_cb', prob, m)
        F = ca.Function('F', [x_sym], [f_cb(x_sym)])
        G = ca.Function('G', [x_sym], [g_cb(x_sym)])

        nlp = {'x': x_sym, 'f': F(x_sym), 'g': G(x_sym)}

        opts = dict(self.casadi_options)
        opts.update({'ipopt': dict(self.ipopt_options)})
        if options:
            # Allow per-call override / extend
            ip = opts.setdefault('ipopt', {})
            ip.update(options)

        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
        sol = solver(x0=x0,
                     lbx=np.asarray(xlb, dtype=float).reshape(-1),
                     ubx=np.asarray(xub, dtype=float).reshape(-1),
                     lbg=np.asarray(glb, dtype=float).reshape(-1),
                     ubg=np.asarray(gub, dtype=float).reshape(-1))
        x_opt = np.array(sol['x']).reshape(-1)
        f_opt = float(np.array(sol['f']).reshape(()))
        stats = solver.stats()
        success = bool(stats.get('success', False))
        status = stats.get('return_status', 'unknown')
        iters = stats.get('iter_count', None)
        return SolveResult(x=x_opt, obj=f_opt, status=status, success=success, n_iter=iters, info=stats)
