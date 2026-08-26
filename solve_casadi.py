#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solve_terminal_only_casadi_modified.py
Terminal-only demo of anti-crossing wheel path optimization with CasADi.

- Soft wheelbase distance limits (hinge penalties on planar distances)
- Anti-crossing: penalize only reverse orientation via cosine (yaw-dominant, terrain-normal adaptive)
- Regularization & smoothing
- Two-stage homotopy on weights for robustness

Author: ChatGPT (for Yongqi)
"""
import numpy as np
import casadi as ca
import matplotlib.pyplot as plt


# -----------------------------
# User parameters
# -----------------------------
N_COEFF   = 5
K_CONSTR  = 80

# 轮–中心距离范围 [0.4, 1.0]，L_MIN/L_MAX 只用来定义名义 lateral 距离 L0
L_MIN     = 0.40        # min nominal wheel lateral distance
L_MAX     = 1.00        # max nominal wheel lateral distance
R_MIN     = 0.40        # min wheel-to-base 3D distance
R_MAX     = 1.00        # max wheel-to-base 3D distance

AMP       = 0.2         # （原始简单地形已经不用了，可以留着）
ELL       = 0.35

# 中心相对地面高度约 0.8m，允许 [0.6, 1.0]
CLEAR_Z   = 0.80        # base z at s=0 and s=1 (名义高度)
CLEAR_MIN = 0.60        # gap_b >= 0.6
CLEAR_MAX = 1.00        # gap_b <= 1.0

MU = 0.9
NMIN = 1e-3

USE_HARD_CLEARANCE = True

# 中心轨迹：从 (2.5, 1.0) 到 (2.5, -3.5)
# yaw 固定朝 -y 方向（-pi/2），这样左右轮沿 x 方向展开
XB0, YB0, PSI0 = 2.5,  1.0, -np.pi/2
XB1, YB1, PSI1 = 2.5, -3.5, -np.pi/2

L0 = 0.5*(L_MIN + L_MAX)


# --------------- RBF terrain parameters (filled from JSON) ---------------
RBF_CENTERS = None   # (N,2) numpy array
RBF_WEIGHTS = None   # (N,) numpy array
RBF_ELL = None       # scalar length scale


# -----------------------------
# Polynomial basis utilities
# -----------------------------
def phi_poly(s, n):
    """Return Vandermonde-like basis [1, s, s^2, ..., s^{n-1}] for each s (vectorized)."""
    s = np.asarray(s).reshape(-1, 1)
    Phi = np.hstack([s**k for k in range(n)])
    return Phi

def lsq_fit(s, y, n):
    """Least squares polynomial fit of degree n-1 for data (s, y)."""
    s = np.asarray(s).ravel()
    y = np.asarray(y).ravel()
    Phi = phi_poly(s, n)
    # Use pseudo-inverse to get stable coefficients
    c = np.linalg.lstsq(Phi, y, rcond=None)[0]
    return c

# -----------------------------
# Terrain model & gradient (RBF version)
# -----------------------------
def terrain_height_casadi(x, y):
    """
    x, y 可以是 CasADi 标量或向量，这里都当成“按元素运算”，
    对每个 RBF center 累加 w_i * exp(-0.5 * r^2 / ell^2).
    """
    assert RBF_CENTERS is not None, "RBF_CENTERS not set; call from __main__ first."
    ell2 = RBF_ELL ** 2

    h = 0
    for (cx, cy), w in zip(RBF_CENTERS, RBF_WEIGHTS):
        dx = x - cx
        dy = y - cy
        r2 = dx * dx + dy * dy
        h = h + w * ca.exp(-0.5 * r2 / ell2)
    return h


def terrain_grad_casadi(x, y):
    """
    ∂h/∂x = Σ w_i * exp(-0.5*r^2/ell^2) * (-(x-cx)/ell^2)
    ∂h/∂y 类似。
    """
    assert RBF_CENTERS is not None, "RBF_CENTERS not set; call from __main__ first."
    ell2 = RBF_ELL ** 2

    dhdx = 0
    dhdy = 0
    for (cx, cy), w in zip(RBF_CENTERS, RBF_WEIGHTS):
        dx = x - cx
        dy = y - cy
        r2 = dx * dx + dy * dy
        g = ca.exp(-0.5 * r2 / ell2)
        dhdx = dhdx + w * g * (-(dx) / ell2)
        dhdy = dhdy + w * g * (-(dy) / ell2)
    return dhdx, dhdy


class RBFTerrain:
    def __init__(self, centers, weights, ell, domain=None, hard_clip=False, eps=1e-8):
        self.centers = np.asarray(centers, dtype=float)   # (N,2)
        self.weights = np.asarray(weights, dtype=float).ravel()  # (N,)
        self.ell = float(ell)
        self.domain = domain
        self.hard_clip = hard_clip
        self.eps = eps

    def h(self, x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if self.hard_clip and self.domain is not None:
            lo, hi = self.domain
            x = np.clip(x, lo, hi)
            y = np.clip(y, lo, hi)

        H = np.zeros_like(x, dtype=float)
        ell2 = self.ell ** 2
        for (cx, cy), w in zip(self.centers, self.weights):
            dx = x - cx
            dy = y - cy
            r2 = dx * dx + dy * dy
            H += w * np.exp(-0.5 * r2 / ell2)
        return H


# -----------------------------
# Build optimization
# -----------------------------
def build_and_solve(seed=None):
    if seed is not None:
        np.random.seed(seed)

    # Collocation points along s ∈ [0,1]
    Kc = K_CONSTR
    s_c = np.linspace(0.0, 1.0, Kc)
    Phi = phi_poly(s_c, N_COEFF)                # (Kc x N_COEFF)
    # Endpoint evaluation matrices (s=0, s=1)
    Phi0 = phi_poly([0.0], N_COEFF)             # (1 x N_COEFF)
    Phi1 = phi_poly([1.0], N_COEFF)             # (1 x N_COEFF)

    opti = ca.Opti()

    # Polynomial coefficients (variables)
    cbx  = opti.variable(N_COEFF)   # base x(s)
    cby  = opti.variable(N_COEFF)   # base y(s)
    cbz  = opti.variable(N_COEFF)   # base z(s)
    cpsi = opti.variable(N_COEFF)   # base yaw(s)

    clx  = opti.variable(N_COEFF)   # left wheel x(s)
    cly  = opti.variable(N_COEFF)   # left wheel y(s)
    crx  = opti.variable(N_COEFF)   # right wheel x(s)
    cry  = opti.variable(N_COEFF)   # right wheel y(s)

    # Convenient evaluators
    def eval_curve(coeff):
        return ca.mtimes(Phi, coeff)  # (Kc,)

    def eval_curve_0(coeff):
        return ca.mtimes(Phi0, coeff) # (1,)

    def eval_curve_1(coeff):
        return ca.mtimes(Phi1, coeff) # (1,)

    # Evaluate along s_c
    Xb_c  = eval_curve(cbx);  Yb_c  = eval_curve(cby);  Zb_c  = eval_curve(cbz)
    Psib_c= eval_curve(cpsi)

    Xl_c  = eval_curve(clx);  Yl_c  = eval_curve(cly)
    Xr_c  = eval_curve(crx);  Yr_c  = eval_curve(cry)

    # Terrain for wheels (for normals only; height is not constrained here)
    dhx_l, dhy_l = terrain_grad_casadi(Xl_c, Yl_c)
    dhx_r, dhy_r = terrain_grad_casadi(Xr_c, Yr_c)
    h_b          = terrain_height_casadi(Xb_c, Yb_c)    # base terrain height along path

    # -----------------------------
    # Objective assembly (with parameters for homotopy)
    # -----------------------------
    W_CROSS = opti.parameter()  # penalty for anti-crossing (only reverse orientation)
    W_LIM   = opti.parameter()  # penalty for violating [L_MIN, L_MAX]
    W_SMO   = opti.parameter()  # smoothing weight
    W_REG   = opti.parameter()  # coefficient regularization
    W_CLR   = opti.parameter()  # soft clearance penalty

    f = 0.0
    eps = 1e-9
    hinge_pos = lambda z: 0.5*(z + ca.fabs(z))   # max(z,0)

    # --- (A) Soft wheel distance limits in XY (hinge penalties) ---
    vLx = Xb_c - Xl_c; vLy = Yb_c - Yl_c; vLz = Zb_c - terrain_height_casadi(Xl_c, Yl_c)
    vRx = Xb_c - Xr_c; vRy = Yb_c - Yr_c; vRz = Zb_c - terrain_height_casadi(Xr_c, Yr_c)
    dL = ca.sqrt(vLx*vLx + vLy*vLy + vLz*vLz + eps)
    dR = ca.sqrt(vRx*vRx + vRy*vRy + vRz*vRz + eps)
    hinge_pos = lambda z: 0.5*(z + ca.fabs(z))   # max(z,0)

    f += (W_LIM / Kc) * ca.sumsqr(hinge_pos(R_MIN - dL))
    f += (W_LIM / Kc) * ca.sumsqr(hinge_pos(dL - R_MAX))
    f += (W_LIM / Kc) * ca.sumsqr(hinge_pos(R_MIN - dR))
    f += (W_LIM / Kc) * ca.sumsqr(hinge_pos(dR - R_MAX))

    # --- (A2) Left-Right planar separation (anti-overlap) ---
    # 平面距离 dLR = || pL_xy - pR_xy ||
    dLR = ca.sqrt((Xl_c - Xr_c)**2 + (Yl_c - Yr_c)**2 + eps)

    TRACK_MIN = 0.3 * L_MIN      # 允许的最小左右间距（可调，比如0.6*L0）
    TRACK_MAX = L_MAX      # 可选，上界（防止过分张开）
    W_TRACK   = 100.0         # 软惩罚权重（也可做成 opti.parameter()）

    # 软惩罚版本（推荐，数值更稳）：
    # f += (W_TRACK / Kc) * ca.sumsqr(hinge_pos(TRACK_MIN - dLR))
    # f += (W_TRACK / Kc) * ca.sumsqr(hinge_pos(dLR - TRACK_MAX))

    # hard constraint version
    opti.subject_to( (Xl_c - Xr_c)**2 + (Yl_c - Yr_c)**2 >= (TRACK_MIN)**2 )

    # f += (W_TRACK / Kc) * ca.sumsqr(dLR - L0)

    # --- (B) Stable anti-crossing (yaw-dominant lateral, terrain-normal adaptive); only penalize reverse ---
    # yaw lateral unit (in XY): u_lat_yaw = [-sinψ, cosψ]
    ux_yaw = -ca.sin(Psib_c)
    uy_yaw =  ca.cos(Psib_c)

    # terrain normals at wheels (each)
    nx_l = -dhx_l; ny_l = -dhy_l; nz_l = 1.0
    nx_r = -dhx_r; ny_r = -dhy_r; nz_r = 1.0
    nl = ca.sqrt(nx_l*nx_l + ny_l*ny_l + nz_l*nz_l + eps)
    nr = ca.sqrt(nx_r*nx_r + ny_r*ny_r + nz_r*nz_r + eps)
    nx_l, ny_l, nz_l = nx_l/nl, ny_l/nl, nz_l/nl
    nx_r, ny_r, nz_r = nx_r/nr, ny_r/nr, nz_r/nr

    # average normal, normalize
    nx_m = 0.5*(nx_l + nx_r)
    ny_m = 0.5*(ny_l + ny_r)
    nz_m = 0.5*(nz_l + nz_r)
    nm = ca.sqrt(nx_m*nx_m + ny_m*ny_m + nz_m*nz_m + eps)
    nx_m, ny_m, nz_m = nx_m/nm, ny_m/nm, nz_m/nm

    # fwd = [cosψ, sinψ, 0], lateral from normal: v_lat_n = fwd × n_m (take XY and normalize)
    fx = ca.cos(Psib_c); fy = ca.sin(Psib_c)
    vlatx_n =  fy*nz_m
    vlaty_n = -fx*nz_m
    nxy = ca.sqrt(vlatx_n*vlatx_n + vlaty_n*vlaty_n + eps)
    ux_n = vlatx_n / nxy
    uy_n = vlaty_n / nxy

    # adaptive blend: if terrain is sloped (nz small), trust normal more; otherwise just yaw
    beta = hinge_pos(0.6 - nz_m) / 0.6   # nz<0.6 => start to use normal
    beta = ca.fmax(0, ca.fmin(1, beta))

    ux_lat = (1 - beta)*ux_yaw + beta*ux_n
    uy_lat = (1 - beta)*uy_yaw + beta*uy_n
    ulat = ca.sqrt(ux_lat*ux_lat + uy_lat*uy_lat + eps)
    ux_lat /= ulat; uy_lat /= ulat

    # actual right->left unit vector in XY
    vrlx = Xl_c - Xr_c; vrly = Yl_c - Yr_c
    nrl = ca.sqrt(vrlx*vrlx + vrly*vrly + eps)
    vx = vrlx / nrl; vy = vrly / nrl

    cosang = vx*ux_lat + vy*uy_lat               # in [-1,1]
    neg_part = hinge_pos(-cosang)                # = max(-cos, 0) -> only reverse penalized
    f += (W_CROSS / Kc) * ca.sumsqr(neg_part)

    # --- (C) Smoothing (discrete second difference) & coefficient regularization ---
    
    def smooth_penalty(vec):
        # vec: (Kc,1)
        v1 = vec[0:-2]
        v2 = vec[1:-1]
        v3 = vec[2:  ]
        return ca.sumsqr(v3 - 2*v2 + v1)

    f += (W_SMO / Kc) * (
        smooth_penalty(Xl_c) + smooth_penalty(Yl_c) +
        smooth_penalty(Xr_c) + smooth_penalty(Yr_c) +
        smooth_penalty(Xb_c) + smooth_penalty(Yb_c)
    )

    coef_stack = ca.vertcat(cbx, cby, cbz, cpsi, clx, cly, crx, cry)
    f += W_REG * ca.sumsqr(coef_stack)

    # -----------------------------
    # Endpoint boundary conditions (keep light and reasonable)
    # -----------------------------
    # Base anchors at endpoints
    opti.subject_to(eval_curve_0(cbx)  == XB0)
    opti.subject_to(eval_curve_0(cby)  == YB0)
    opti.subject_to(eval_curve_0(cpsi) == PSI0)
    opti.subject_to(eval_curve_0(cbz)  == CLEAR_Z)

    opti.subject_to(eval_curve_1(cbx)  == XB1)
    opti.subject_to(eval_curve_1(cby)  == YB1)
    opti.subject_to(eval_curve_1(cpsi) == PSI1)
    opti.subject_to(eval_curve_1(cbz)  == CLEAR_Z)

    # Wheel endpoints laterally offset from base by L0 using yaw lateral [-sinψ, cosψ]
    # s = 0
    ux0 = -np.sin(PSI0); uy0 = np.cos(PSI0)
    opti.subject_to(eval_curve_0(clx) == XB0 + (-L0/2.0)*ux0)
    opti.subject_to(eval_curve_0(cly) == YB0 + (-L0/2.0)*uy0)
    opti.subject_to(eval_curve_0(crx) == XB0 + ( L0/2.0)*ux0)
    opti.subject_to(eval_curve_0(cry) == YB0 + ( L0/2.0)*uy0)

    # s = 1
    ux1 = -np.sin(PSI1); uy1 = np.cos(PSI1)
    opti.subject_to(eval_curve_1(clx) == XB1 + (-L0/2.0)*ux1)
    opti.subject_to(eval_curve_1(cly) == YB1 + (-L0/2.0)*uy1)
    opti.subject_to(eval_curve_1(crx) == XB1 + ( L0/2.0)*ux1)
    opti.subject_to(eval_curve_1(cry) == YB1 + ( L0/2.0)*uy1)


    # --- (D) Base-to-terrain clearance (with both lower and upper bounds) ---
    gap_b = Zb_c - h_b       # positive if base is above terrain

    if USE_HARD_CLEARANCE:
        # 硬约束：CLEAR_MIN ≤ gap_b ≤ CLEAR_MAX
        opti.subject_to(gap_b - CLEAR_MIN >= 0)
        opti.subject_to(CLEAR_MAX - gap_b >= 0)
    else:
        # 软惩罚：下界和上界都用 hinge
        # 下界：gap_b < CLEAR_MIN → max(CLEAR_MIN - gap_b, 0)
        # 上界：gap_b > CLEAR_MAX → max(gap_b - CLEAR_MAX, 0)
        f += (W_CLR / Kc) * ca.sumsqr(hinge_pos(CLEAR_MIN - gap_b))  # too low
        f += (W_CLR / Kc) * ca.sumsqr(hinge_pos(gap_b - CLEAR_MAX))  # too high
    
    # --- (E) Quasi-static friction cone at wheel contacts ---
    # 左轮：将 vL 分解到法向与切向
    vL_n = vLx*nx_l + vLy*ny_l + vLz*nz_l                  # 标量(逐点)
    tLx  = vLx - vL_n*nx_l
    tLy  = vLy - vL_n*ny_l
    tLz  = vLz - vL_n*nz_l
    tL_sq = tLx*tLx + tLy*tLy + tLz*tLz                    # ||v_t||^2

    # 右轮
    vR_n = vRx*nx_r + vRy*ny_r + vRz*nz_r
    tRx  = vRx - vR_n*nx_r
    tRy  = vRy - vR_n*ny_r
    tRz  = vRz - vR_n*nz_r
    tR_sq = tRx*tRx + tRy*tRy + tRz*tRz

    # 摩擦锥硬约束：||v_t|| <= mu * v_n  <=>  t_sq <= (mu^2) * v_n^2，且 v_n >= NMIN
    opti.subject_to(tL_sq <= (MU*MU) * (vL_n*vL_n))
    opti.subject_to(tR_sq <= (MU*MU) * (vR_n*vR_n))
    opti.subject_to(vL_n >= NMIN)
    opti.subject_to(vR_n >= NMIN)

    # -----------------------------
    # Initial guess (yaw=0 lateral offsets along linear base ref)
    # -----------------------------
    S_init = np.linspace(0, 1, 7)

    # 中心轨迹在 x=2.5 上从 y=1 到 y=-3.5
    xb_ref = np.linspace(XB0, XB1, len(S_init))   # 其实就是常数 2.5
    yb_ref = np.linspace(YB0, YB1, len(S_init))

    # 基座名义高度 CLEAR_Z，对应 gap_b 大约 0.8
    zb_ref  = CLEAR_Z * np.ones_like(xb_ref)

    # yaw 固定为 PSI0（= PSI1 = -pi/2）
    psi_ref = PSI0 * np.ones_like(xb_ref)

    # 左右轮初值：沿 lateral 方向 [-sinψ, cosψ] 展开
    ux_ref = -np.sin(psi_ref)
    uy_ref =  np.cos(psi_ref)

    xl_ref = xb_ref + (-L0/2.0)*ux_ref
    yl_ref = yb_ref + (-L0/2.0)*uy_ref
    xr_ref = xb_ref + ( L0/2.0)*ux_ref
    yr_ref = yb_ref + ( L0/2.0)*uy_ref

    cbx0  = lsq_fit(S_init, xb_ref,  N_COEFF)
    cby0  = lsq_fit(S_init, yb_ref,  N_COEFF)
    cbz0  = lsq_fit(S_init, zb_ref,  N_COEFF)
    cpsi0 = lsq_fit(S_init, psi_ref, N_COEFF)
    clx0  = lsq_fit(S_init, xl_ref,  N_COEFF)
    cly0  = lsq_fit(S_init, yl_ref,  N_COEFF)
    crx0  = lsq_fit(S_init, xr_ref,  N_COEFF)
    cry0  = lsq_fit(S_init, yr_ref,  N_COEFF)

    # Visualize the initial guess before solving
    try:
        print("Visualizing initial guess...")
        Sd_vis = np.linspace(0.0, 1.0, 201)
        Phi_vis = phi_poly(Sd_vis, N_COEFF)
        Xb0_d = Phi_vis @ cbx0
        Yb0_d = Phi_vis @ cby0
        Zb0_d = Phi_vis @ cbz0
        Xl0_d = Phi_vis @ clx0
        Yl0_d = Phi_vis @ cly0
        Xr0_d = Phi_vis @ crx0
        Yr0_d = Phi_vis @ cry0

        traj_init = ((Xb0_d, Yb0_d, Zb0_d),
                    (Xl0_d, Yl0_d),
                    (Xr0_d, Yr0_d))

        # 用真实地形范围
        terrain_vis = RBFTerrain(
            centers=RBF_CENTERS,
            weights=RBF_WEIGHTS,
            ell=RBF_ELL
        )

        # 自动根据 RBF 范围决定 domain（不传 domain 参数）
        plot_on_terrain(traj_init, terrain=terrain_vis)

    except Exception as e:
        print("Initial-guess visualization failed:", e)

    opti.set_initial(cbx,  cbx0)
    opti.set_initial(cby,  cby0)
    opti.set_initial(cbz,  cbz0)
    opti.set_initial(cpsi, cpsi0)
    opti.set_initial(clx,  clx0)
    opti.set_initial(cly,  cly0)
    opti.set_initial(crx,  crx0)
    opti.set_initial(cry,  cry0)

    # -----------------------------
    # Initial constraint diagnostics
    # -----------------------------
    print("Running initial constraint diagnostics ...")

    try:
        g0 = opti.debug.value(opti.g)   # ALWAYS use debug.value()
        if isinstance(g0, np.ndarray) and g0.size > 0:
            print(f"[diag] Initial max |g| = {float(np.max(np.abs(g0))):.3e}")
        else:
            print("[diag] g0 empty or invalid")
    except Exception as e:
        print("[diag] Could not evaluate initial constraints:", e)



    # -----------------------------
    # Solver options
    # -----------------------------
    print("Setting up solver...")
    MAX_IPOPT_ITER = 500
    p_opts = {
        'expand': True,
        'print_time': True,   # CasADi 打印耗时
    }

    # 注意：这里是 Ipopt 原始选项名，没有 'ipopt.' 前缀
    s_opts = {
        'sb': 'no',            # 不静默 (silent barrier)
        'print_level': 5,      # 0~12
        'max_iter': MAX_IPOPT_ITER,
        'tol': 1e-4,
        'acceptable_tol': 5e-4,
        'acceptable_iter': 10,
        'linear_solver': 'mumps',
        'output_file': 'ipopt_log.txt',
    }

    opti.solver('ipopt', p_opts, s_opts)
    print(f"Solver configured: Ipopt, max_iter = {MAX_IPOPT_ITER}")




    # -----------------------------
    # Homotopy: stage 1 & stage 2
    # -----------------------------
    def solve_with_weights(w_cross, w_lim, w_smo, w_reg, w_clr, label=""):
        print(f"\n===== Solving {label} =====")
        print(f"  W_CROSS={w_cross}, W_LIM={w_lim}, W_SMO={w_smo}, W_REG={w_reg}, W_CLR={w_clr}")
        opti.set_value(W_CROSS, w_cross)
        opti.set_value(W_LIM,   w_lim)
        opti.set_value(W_SMO,   w_smo)
        opti.set_value(W_REG,   w_reg)
        opti.set_value(W_CLR,   w_clr)

        try:
            sol = opti.solve()
            f_val = float(sol.value(f))
            print(f"[{label}] SUCCESS. Objective = {f_val:.6e}")
            # 这里 stats 是安全的，因为 solve 成功了
            stats = opti.stats()
            print(f"[{label}] return_status = {stats.get('return_status')}, "
                  f"iter = {stats.get('iter_count')}")
            return sol

        except RuntimeError as e:
            print(f"[{label}] FAILED with RuntimeError: {e}")
            # 这里不能直接 stats()，可能还没 solve 成功
            try:
                stats = opti.stats()
                print("  return_status:", stats.get("return_status"))
                print("  success      :", stats.get("success"))
                print("  iter_count   :", stats.get("iter_count"))
            except Exception as e_stats:
                print("  (no stats available, opti not solved yet)", e_stats)

            # 尝试看一下此刻的约束违反程度
            try:
                g_val = opti.debug.value(opti.g)
                max_g = float(np.max(np.abs(g_val))) if g_val.size > 0 else 0.0
                print("  max |g| at failure point:", max_g)
            except Exception as ee:
                print("  (could not evaluate opti.g at failure point:", ee, ")")

            # 把错误往外抛，让你在上层看到
            raise



    # Stage 1: gentle weights to get feasibility and direction right
    sol1 = solve_with_weights(
        w_cross=80,
        w_lim=3.0,
        w_smo=0.05,
        w_reg=1e-3,
        w_clr=2.0 if not USE_HARD_CLEARANCE else 0.0,
        label="stage 1"
    )

    # Warm-start stage 2 with previous solution (use 2-arg form)
    x1 = sol1.value(opti.x)
    try:
        opti.set_initial(opti.x, x1)
        # also warm-start duals if available
        opti.set_initial(opti.lam_g, sol1.value(opti.lam_g))
        opti.set_initial(opti.lam_x, sol1.value(opti.lam_x))
    except Exception as e:
        # Fall back to just primal variables
        opti.set_initial(opti.x, x1)


    # Stage 2: ramp-up
    sol2 = solve_with_weights(
        w_cross=80,
        w_lim=10.0,
        w_smo=0.10,
        w_reg=1e-3,
        w_clr=8.0 if not USE_HARD_CLEARANCE else 0.0,
        label="stage 2"
    )

    # Extract solution
    cbx_v  = sol2.value(cbx);  cby_v  = sol2.value(cby)
    cbz_v  = sol2.value(cbz);  cpsi_v = sol2.value(cpsi)
    clx_v  = sol2.value(clx);  cly_v  = sol2.value(cly)
    crx_v  = sol2.value(crx);  cry_v  = sol2.value(cry)

    # Evaluate dense path for saving
    Sd = np.linspace(0, 1, 201)
    Phi_d = phi_poly(Sd, N_COEFF)
    def eval_dense(c):
        return Phi_d @ c

    Xb_d = eval_dense(cbx_v);  Yb_d = eval_dense(cby_v)
    Zb_d = eval_dense(cbz_v);  Psib_d = eval_dense(cpsi_v)
    Xl_d = eval_dense(clx_v);  Yl_d = eval_dense(cly_v)
    Xr_d = eval_dense(crx_v);  Yr_d = eval_dense(cry_v)

    # Save to npz
    # np.savez('/mnt/data/optimized_paths.npz',
    #          s=Sd,
    #          Xb=Xb_d, Yb=Yb_d, Zb=Zb_d, Psi=Psib_d,
    #          Xl=Xl_d, Yl=Yl_d, Xr=Xr_d, Yr=Yr_d,
    #          cbx=cbx_v, cby=cby_v, cbz=cbz_v, cpsi=cpsi_v,
    #          clx=clx_v, cly=cly_v, crx=crx_v, cry=cry_v)

    # Print brief summary
    print('Solve OK. Coeff L2 norm:',
          float(np.linalg.norm(np.concatenate([cbx_v, cby_v, cbz_v, cpsi_v, clx_v, cly_v, crx_v, cry_v]))))
    return Sd, (Xb_d, Yb_d, Zb_d), (Xl_d, Yl_d), (Xr_d, Yr_d)


def plot_on_terrain(trajs, terrain=None, domain=None):
    if terrain is None:
        raise ValueError("terrain must be provided for plotting RBF terrain")
    (Xb, Yb, Zb), (Xl, Yl), (Xr, Yr) = trajs

    # ---- 自动根据 RBF center 决定地形范围 ----
    if domain is None:
        # centers 是 (N,2) 的数组
        xmin, ymin = terrain.centers.min(axis=0)
        xmax, ymax = terrain.centers.max(axis=0)
        margin = 0.2     # 适当留一点边界
        xmin -= margin
        xmax += margin
        ymin -= margin
        ymax += margin
    else:
        # 允许 domain=(xmin, xmax, ymin, ymax)
        if len(domain) == 2:
            xmin, xmax = domain
            ymin, ymax = domain
        elif len(domain) == 4:
            xmin, xmax, ymin, ymax = domain
        else:
            raise ValueError("domain must be None, (min,max) or (xmin,xmax,ymin,ymax)")

    xs = np.linspace(xmin, xmax, 161)
    ys = np.linspace(ymin, ymax, 161)
    X, Y = np.meshgrid(xs, ys)
    Z = terrain.h(X, Y)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=2, cstride=2, linewidth=0, alpha=0.8)

    # Wheels' z come from terrain
    Zl = terrain.h(Xl, Yl)
    Zr = terrain.h(Xr, Yr)

    ax.plot(Xb, Yb, Zb, linewidth=2.5, label='base (z poly)')
    ax.plot(Xl, Yl, Zl, linewidth=1.8, label='left (z=h)')
    ax.plot(Xr, Yr, Zr, linewidth=1.8, label='right (z=h)')

    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title('Power-basis trajectories with boundary constraints on RBF terrain')
    ax.legend(loc='upper left')

    xr, yr, zr = np.ptp(xs), np.ptp(ys), np.ptp(Z)
    ax.set_box_aspect((xr, yr, zr))
    ax.view_init(elev=35, azim=-60)
    plt.tight_layout()

    # 下面这几个子图保持不动
    fig2, axs = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
    axs[0].plot(Xb, label='Xb'); axs[0].set_ylabel('Xb'); axs[0].legend()
    axs[1].plot(Yb, label='Yb'); axs[1].set_ylabel('Yb'); axs[1].legend()
    axs[2].plot(Zb, label='Zb'); axs[2].set_ylabel('Zb'); axs[2].set_xlabel('sample index'); axs[2].legend()
    axs[2].set_ylim(0.0, 1.2)
    fig2.suptitle('Base trajectory components')
    plt.tight_layout()
    plt.show()

    
def read_from_json(file_name):
    import json


    # 读取 JSON 文件
    with open(file_name, "r") as f:
        data = json.load(f)

    # ------ 读取 timestamp ------
    timestamp = data["timestamp"]
    print("timestamp:", timestamp)

    # ------ 读取 pose ------
    pos = data["pose"]["position"]
    ori = data["pose"]["orientation"]
    print("position:", pos)
    print("orientation:", ori)

    # ------ mean_weight ------
    mean_weight = data["mean_weight"]
    # print("mean weight:", mean_weight)

    # ------ cur_rbf_grid (N x 2 数组) ------
    cur_rbf_grid = data["cur_rbf_grid"]    # list of [x, y]
    # print("cur_rbf_grid:", cur_rbf_grid)

    # ------ rbf_weight (list<double>) ------
    rbf_weight = data["rbf_weight"]
    # print("rbf_weight:", rbf_weight)

    return timestamp, pos, ori, mean_weight, np.array(cur_rbf_grid), np.array(rbf_weight)

def visualize_rbf_fit(centers, weights, ell=0.14, margin=0.2, nx=240, ny=240, batch=1024):
        import jax
        """
        Visualize an RBF fit given:
         - centers: (N,2) array of [x,y] RBF centers
         - weights: (N,) array of amplitudes
         - ell: length scale of Gaussian RBF
        The implementation automatically adapts to N and computes in chunks to limit memory.
        Uses numpy by default if jax is not available.
        """
        centers = np.asarray(centers, dtype=float)
        weights = np.asarray(weights, dtype=float).ravel()
        if centers.ndim != 2 or centers.shape[1] != 2:
            raise ValueError("centers must be (N,2)")
        if centers.shape[0] != weights.shape[0]:
            print("centers shape:", centers.shape)
            print("weights shape:", weights.shape)
            raise ValueError("centers and weights must have same length")

        # ROI from centers with margin
        xmin, ymin = centers.min(axis=0) - margin
        xmax, ymax = centers.max(axis=0) + margin

        xs = np.linspace(xmin, xmax, nx)
        ys = np.linspace(ymin, ymax, ny)
        X, Y = np.meshgrid(xs, ys)

        H = np.zeros_like(X, dtype=float)
        N = centers.shape[0]

        # Try to use jax if available (vectorized, potentially JIT on accelerator).
        use_jax = False
        try:
            import jax.numpy as jnp
            use_jax = True
            print("Using JAX for RBF evaluation.")
        except Exception:
            use_jax = False
            print("JAX not available; falling back to NumPy for RBF evaluation.")

        if use_jax:
            # JAX path: compute contributions in batches to reduce memory
            X_j = jnp.array(X)
            Y_j = jnp.array(Y)
            centers_j = jnp.array(centers)
            weights_j = jnp.array(weights)

            def batch_sum(i0, i1, acc):
                Cb = centers_j[i0:i1]
                Wb = weights_j[i0:i1]
                dx = X_j[..., None] - Cb[:, 0]
                dy = Y_j[..., None] - Cb[:, 1]
                r2 = dx * dx + dy * dy
                contrib = jnp.exp(-0.5 * r2 / (ell * ell)) * Wb[None, None, :]
                return acc + contrib.sum(axis=-1)

            acc = jnp.zeros_like(X_j)
            for i in range(0, N, batch):
                i1 = min(N, i + batch)
                acc = batch_sum(i, i1, acc)
            H = np.array(acc)
        else:
            # NumPy path: compute contributions in batches
            for i in range(0, N, batch):
                i1 = min(N, i + batch)
                Cb = centers[i:i1]
                Wb = weights[i:i1]
                # dx shape -> (ny, nx, m)
                dx = X[..., None] - Cb[None, :, 0].T  # trick to get broadcasting safe shape
                # to keep code robust use explicit broadcasting
                dx = X[..., None] - Cb[None, :, 0].T  # (ny, nx, m)
                dy = Y[..., None] - Cb[None, :, 1].T
                r2 = dx * dx + dy * dy
                # weight broadcasting (m,) -> (1,1,m)
                contrib = np.exp(-0.5 * r2 / (ell * ell)) * Wb[None, None, :]
                H += contrib.sum(axis=-1)

        # Simple plotting: 3D surface + centers overlay + heatmap
        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax1.plot_surface(X, Y, H, cmap='viridis', linewidth=0, antialiased=False, alpha=0.9)
        # compute z at centers for scatter: evaluate RBF at center positions (cheap)
        z_centers = np.zeros(centers.shape[0], dtype=float)
        for i in range(0, N, batch):
            i1 = min(N, i + batch)
            Cb = centers[i:i1]
            Wb = weights[i:i1]
            # centers[:, 0:1] -> shape (N,1); Cb[:,0][None,:] -> shape (1,m)
            # use 1:2 slice for the second column (avoid empty slice 1:1)
            dx = centers[:, 0:1] - Cb[:, 0][None, :]
            dy = centers[:, 1:2] - Cb[:, 1][None, :]
            # dx,dy shape -> (N, m)
            r2c = dx * dx + dy * dy
            z_centers += (np.exp(-0.5 * r2c / (ell * ell)) * Wb[None, :]).sum(axis=1)
        ax1.scatter(centers[:, 0], centers[:, 1], z_centers, c='r', s=8, depthshade=True)
        ax1.set_xlabel('x'); ax1.set_ylabel('y'); ax1.set_zlabel('z')
        ax1.set_title('RBF fitted terrain (3D surface)')

        ax2 = fig.add_subplot(1, 2, 2)
        im = ax2.imshow(H, extent=(xmin, xmax, ymin, ymax), origin='lower', cmap='viridis')
        ax2.scatter(centers[:, 0], centers[:, 1], c='r', s=6)
        ax2.set_title('RBF height heatmap')
        ax2.set_xlabel('x'); ax2.set_ylabel('y')
        fig.colorbar(im, ax=ax2, shrink=0.8)
        plt.tight_layout()
        plt.show()
        return X, Y, H

   

if __name__ == '__main__':
    # read settings from json
    file_name = "/home/yizhe/trajopt/result/RBF/1753771011_862037182.json"
    timestamp, pos, ori, mean_weight, cur_rbf_grid, rbf_weight = read_from_json(file_name)

    sigma = 0.14

    # 设置全局 RBF 参数，供 CasADi 使用 (module-level variables; no `global` needed here)
    RBF_CENTERS = cur_rbf_grid        # (N,2) numpy
    RBF_WEIGHTS = rbf_weight.ravel()  # (N,)
    RBF_ELL = sigma

    try:
        visualize_rbf_fit(cur_rbf_grid, rbf_weight, ell=sigma, margin=0.25, nx=240, ny=240, batch=1024)
    except Exception as e:
        print("RBF visualization failed:", e)
    # exit()

    Sd, base, left, right = build_and_solve()
    plot_on_terrain((base, left, right),
                terrain=RBFTerrain(centers=RBF_CENTERS,
                                   weights=RBF_WEIGHTS,
                                   ell=RBF_ELL))
