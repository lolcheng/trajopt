import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from utils.utils import read_from_json, plot_on_terrain


# -----------------------------
# User parameters
# -----------------------------
N_COEFF   = 5
N_SAMPLE  = 100

# 两轮距离范围
WHEELS_MIN     = 0.1
WHEELS_MAX     = 0.8
WHEELS_NOMINAL = 0.5

# 中心相对地面高度
CLEAR_Z   = 0.8     
CLEAR_MIN = 0.4
CLEAR_MAX = 1.2

# 轮-中心欧式距离范围
WHEEL_BASE_MIN = 0.4
WHEEL_BASE_MAX = 1.2
WHEEL_BASE_NOMINAL = 0.8

# 起点&终点
XB0, YB0, PSI0 = 2.5,  1.0, -np.pi/2
XB1, YB1, PSI1 = 2.5, -3.5, -np.pi/2

# RBF terrain parameters (filled from JSON)
RBF_CENTERS = None      # (N,2) numpy array
RBF_WEIGHTS = None      # (N,) numpy array
RBF_SIGMA = None        # scalar length scale


MU_FRICTION   = 0.6    # friction coefficient
MASS          = 1.0    
GRAVITY       = 9.81
K_CONTACT     = 500.0  # contact stiffness

MAX_IPOPT_ITER = 1000
W_YAW_ALIGN = 50.0  # weight for yaw-path alignment objective
W_FORCE_BAL   = 10.0   # 合力平衡残差权重
W_TORQUE_BAL  = 10.0   # 合矩平衡残差权重
W_CONTACT_STIFF = 1000.0  # 惩罚接触力与穿透量弹簧模型不符的权重
W_FRIC_CONE = 1000.0    # 摩擦锥软约束权重

# -----------------------------
# MODULE BOOLEANS
# -----------------------------
VISUAL_INITIAL_GUESS = True
VISUAL_OPT = True
USE_SMOOTHING = True
USE_REGULARIZATION = False
USE_WHEEL_DISTANCE = True
USE_WHEEL_BASE_DISTANCE = True
USE_WHEEL_TERRAIN_DISTANCE = False
USE_ANTI_CROSS = True
USE_QUASI_STATIC = True
USE_FRICTION_CONE = True

# -----------------------------
# Polynomial basis utilities
# -----------------------------
def phi_poly(s, n):
    """
    Generate a Vandermonde-style polynomial basis matrix.

    Parameters
    ----------
    s : array_like
        1-D array or column vector of sample points with shape (m,) or (m, 1).
    n : int
        Number of polynomial basis columns. The returned matrix has n columns
        corresponding to [1, s, s**2, ..., s**(n-1)].

    Returns
    -------
        Design matrix where each row corresponds to a sample in `s` and each
        column corresponds to increasing powers of `s` from 0 to n-1.
        The dtype follows NumPy's broadcasting rules for the input `s`.
    """
    s = np.asarray(s).reshape(-1, 1)
    Phi = np.hstack([s**k for k in range(n)])
    return Phi


def lsq_fit(s, y, n):
    """
    Fit a polynomial (in the basis provided by phi_poly) to data (s, y) by
    solving a linear least-squares problem.

    Parameters
    ----------
    s : array_like
        1-D array of independent variable samples. Will be flattened with
        np.asarray(...).ravel().
    y : array_like
        1-D array of dependent variable samples (observations) corresponding to s.
        Will be flattened with np.asarray(...).ravel().
    n : int
        Number of basis functions to use (equivalently, degree + 1). The design
        matrix is constructed by calling phi_poly(s, n), so the returned
        coefficients correspond to the column ordering of phi_poly.

    Returns
    -------
    c : ndarray, shape (n,)
        Least-squares solution for the coefficient vector. Computed using
        numpy.linalg.lstsq on the design matrix Phi = phi_poly(s, n) and the
        observation vector y. The solution is the minimum-norm least-squares
        solution when the system is underdetermined.
    """
    s = np.asarray(s).ravel()
    y = np.asarray(y).ravel()
    Phi = phi_poly(s, n)
    # Use pseudo-inverse to get stable coefficients
    c = np.linalg.lstsq(Phi, y, rcond=None)[0]
    return c

def smooth_pos(z, eps=1e-3):
    # 近似 max(0, z)，但 C¹ 光滑
    return 0.5 * (z + ca.sqrt(z*z + eps))


# -----------------------------
# RBFTerrain Class
# -----------------------------
class RBFTerrain:
    def __init__(self, centers, weights, sigma, domain=None, hard_clip=False, eps=1e-8):
        self.centers = np.asarray(centers, dtype=float)   # (N,2)
        self.weights = np.asarray(weights, dtype=float).ravel()  # (N,)
        self.sigma = float(sigma)
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
        sigma2 = self.sigma ** 2
        for (cx, cy), w in zip(self.centers, self.weights):
            dx = x - cx
            dy = y - cy
            r2 = dx * dx + dy * dy
            H += w * np.exp(-0.5 * r2 / sigma2)
        return H
    
    def terrain_height_casadi(self, x, y):
        """
        Compute the terrain height at (x, y) using CasADi.
        """
        # Create a CasADi-compatible zero with the same shape/type as x
        H = x * 0
        sigma2 = self.sigma ** 2
        # Use the instance data (self.centers, self.weights) instead of global vars
        for (cx, cy), w in zip(self.centers, self.weights):
            dx = x - float(cx)
            dy = y - float(cy)
            r2 = dx * dx + dy * dy
            H = H + float(w) * ca.exp(-0.5 * r2 / sigma2)
        return H
    
    def terrain_grad_casadi(self, x, y):
        """
        Compute the gradient of the terrain height at (x, y) using CasADi.
        Returns the partial derivatives (Hx, Hy).
        """
        Hx = x * 0
        Hy = y * 0
        sigma2 = self.sigma ** 2
        for (cx, cy), w in zip(self.centers, self.weights):
            dx = x - float(cx)
            dy = y - float(cy)
            r2 = dx * dx + dy * dy
            g = float(w) * ca.exp(-0.5 * r2 / sigma2)
            Hx = Hx + g * (-dx / sigma2)
            Hy = Hy + g * (-dy / sigma2)
        return Hx, Hy

# -----------------------------
# Build optimization
# -----------------------------
def build_and_solve(seed=None):
    if seed is not None:
        np.random.seed(seed)
    
    # Collocation points along s ∈ [0,1]
    n_sample = N_SAMPLE
    s_c = np.linspace(0.0, 1.0, n_sample)
    Phi = phi_poly(s_c, N_COEFF)                # (n_sample x N_COEFF)
    # Endpoint evaluation matrices (s=0, s=1)
    Phi0 = phi_poly([0.0], N_COEFF)             # (1 x N_COEFF)
    Phi1 = phi_poly([1.0], N_COEFF)             # (1 x N_COEFF)

    # -----------------------------
    # Solver options
    # -----------------------------
    opti = ca.Opti()
    print("Setting up solver...")
    p_opts = {
        'expand': True,
        'print_time': True,   # CasADi 打印耗时
    }

    s_opts = {
        'sb': 'no',            # 不静默 (silent barrier)
        'print_level': 5,      # 0~12
        'max_iter': MAX_IPOPT_ITER,
        'tol': 1e-4,
        'acceptable_tol': 5e-4,
        'acceptable_iter': 10,
        'linear_solver': 'mumps',
    }

    opti.solver('ipopt', p_opts, s_opts)
    print(f"Solver configured: Ipopt, max_iter = {MAX_IPOPT_ITER}")

    # Polynomial coefficients (variables)
    cbx  = opti.variable(N_COEFF)   # base x(s)
    cby  = opti.variable(N_COEFF)   # base y(s)
    cbz  = opti.variable(N_COEFF)   # base z(s)
    cpsi = opti.variable(N_COEFF)   # base yaw(s)

    clx  = opti.variable(N_COEFF)   # left wheel x(s)
    cly  = opti.variable(N_COEFF)   # left wheel y(s)
    clz  = opti.variable(N_COEFF)   # left wheel z(s)
    crx  = opti.variable(N_COEFF)   # right wheel x(s)
    cry  = opti.variable(N_COEFF)   # right wheel y(s)
    crz  = opti.variable(N_COEFF)   # right wheel z(s)

    fLx = opti.variable(n_sample)   # left wheel contact force x(s)
    fLy = opti.variable(n_sample)   # left wheel contact force y(s)
    fLz = opti.variable(n_sample)   # left wheel contact force z(s)

    fRx = opti.variable(n_sample)   # right wheel contact force x(s)
    fRy = opti.variable(n_sample)   # right wheel contact force y(s)
    fRz = opti.variable(n_sample)   # right wheel contact force z(s)

    # -----------------------------
    # Trajectory evaluation
    # -----------------------------
    def eval_curve(coeff):
        return ca.mtimes(Phi, coeff)  # (n_sample,)

    def eval_curve_0(coeff):
        return ca.mtimes(Phi0, coeff) # (1,)

    def eval_curve_1(coeff):
        return ca.mtimes(Phi1, coeff) # (1,)
    # base trajectory
    Xb_c  = eval_curve(cbx)
    Yb_c  = eval_curve(cby)
    Zb_c  = eval_curve(cbz)
    Psib_c= eval_curve(cpsi)
    # wheel trajectories
    Xl_c  = eval_curve(clx)
    Yl_c  = eval_curve(cly)
    Zl_c  = eval_curve(clz)
    Xr_c  = eval_curve(crx)
    Yr_c  = eval_curve(cry)
    Zr_c  = eval_curve(crz)

    # -----------------------------
    # Terrain instance
    # -----------------------------
    terrain_instance = RBFTerrain(
        centers=RBF_CENTERS,
        weights=RBF_WEIGHTS,
        sigma=RBF_SIGMA
    )

    # -----------------------------
    # Initial guess
    # -----------------------------
    S_init = np.linspace(0, 1, 7)
    # base
    xb_initial = np.linspace(XB0, XB1, len(S_init))
    yb_initial = np.linspace(YB0, YB1, len(S_init))
    ZB0 = terrain_instance.h(XB0, YB0) + CLEAR_Z
    ZB1 = terrain_instance.h(XB1, YB1) + CLEAR_Z
    zb_initial = np.linspace(ZB0, ZB1, len(S_init))

    psi_initial = np.linspace(PSI0, PSI1, len(S_init))
    # wheel
    ux_initial = -np.sin(psi_initial)
    uy_initial =  np.cos(psi_initial)

    xl_initial = xb_initial + (WHEELS_NOMINAL/2) * ux_initial
    yl_initial = yb_initial + (WHEELS_NOMINAL/2) * uy_initial
    zl_initial = zb_initial - CLEAR_Z

    xr_initial = xb_initial - (WHEELS_NOMINAL/2) * ux_initial
    yr_initial = yb_initial - (WHEELS_NOMINAL/2) * uy_initial
    zr_initial = zb_initial - CLEAR_Z

    # base coeffs
    cbx0  = lsq_fit(S_init, xb_initial,  N_COEFF)
    cby0  = lsq_fit(S_init, yb_initial,  N_COEFF)
    cbz0  = lsq_fit(S_init, zb_initial,  N_COEFF)
    cpsi0 = lsq_fit(S_init, psi_initial, N_COEFF)
    # wheel coeffs
    clx0  = lsq_fit(S_init, xl_initial,  N_COEFF)
    cly0  = lsq_fit(S_init, yl_initial,  N_COEFF)
    clz0  = lsq_fit(S_init, zl_initial,  N_COEFF)
    crx0  = lsq_fit(S_init, xr_initial,  N_COEFF)
    cry0  = lsq_fit(S_init, yr_initial,  N_COEFF)
    crz0  = lsq_fit(S_init, zr_initial,  N_COEFF)

    if VISUAL_INITIAL_GUESS:
        print("-----------------------------")
        print("Initial guess polynomial coefficients:")
        print("cbx0:", cbx0)
        print("cby0:", cby0)
        print("cbz0:", cbz0)
        print("cpsi0:", cpsi0)
        print("clx0:", clx0)
        print("cly0:", cly0)
        print("clz0:", clz0)
        print("crx0:", crx0)
        print("cry0:", cry0)
        print("crz0:", crz0)
        print("-----------------------------")
        print()

    # -----------------------------
    # SDF-based surface normals (from RBF gradient)
    # -----------------------------
    dhx_l, dhy_l = terrain_instance.terrain_grad_casadi(Xl_c, Yl_c)
    dhx_r, dhy_r = terrain_instance.terrain_grad_casadi(Xr_c, Yr_c)

    # phi(x,y,z) = z - h(x,y)，梯度为 (-dh/dx, -dh/dy, 1)
    nLx_raw = -dhx_l
    nLy_raw = -dhy_l
    nLz_raw =  1.0

    nRx_raw = -dhx_r
    nRy_raw = -dhy_r
    nRz_raw =  1.0

    # 单位法向
    eps_n = 1e-6
    nL_norm = ca.sqrt(nLx_raw**2 + nLy_raw**2 + nLz_raw**2 + eps_n)
    nR_norm = ca.sqrt(nRx_raw**2 + nRy_raw**2 + nRz_raw**2 + eps_n)

    nLx = nLx_raw / nL_norm
    nLy = nLy_raw / nL_norm
    nLz = nLz_raw / nL_norm

    nRx = nRx_raw / nR_norm
    nRy = nRy_raw / nR_norm
    nRz = nRz_raw / nR_norm


    # -----------------------------
    # Visualize the initial guess before solving
    # -----------------------------
    if VISUAL_INITIAL_GUESS:
        try:
            print("Visualizing initial guess...")
            Sd_vis = np.linspace(0.0, 1.0, 201)
            Phi_vis = phi_poly(Sd_vis, N_COEFF)
            Xb0_d = Phi_vis @ cbx0
            Yb0_d = Phi_vis @ cby0
            Zb0_d = Phi_vis @ cbz0
            Psi0_d = Phi_vis @ cpsi0
            Xl0_d = Phi_vis @ clx0
            Yl0_d = Phi_vis @ cly0
            Zl0_d = Phi_vis @ clz0
            Xr0_d = Phi_vis @ crx0
            Yr0_d = Phi_vis @ cry0
            Zr0_d = Phi_vis @ crz0

            traj_init = ((Xb0_d, Yb0_d, Zb0_d, Psi0_d),
                        (Xl0_d, Yl0_d, Zl0_d),
                        (Xr0_d, Yr0_d, Zr0_d))

            # 自动根据 RBF 范围决定 domain（不传 domain 参数）
            plot_on_terrain(traj_init, terrain=terrain_instance)

        except Exception as e:
            print("Initial-guess visualization failed:", e)
    
    # -----------------------------
    # Set initial guess
    # -----------------------------
    opti.set_initial(cbx,  cbx0)
    opti.set_initial(cby,  cby0)
    opti.set_initial(cbz,  cbz0)
    opti.set_initial(cpsi, cpsi0)
    opti.set_initial(clx,  clx0)
    opti.set_initial(cly,  cly0)
    opti.set_initial(clz,  clz0)
    opti.set_initial(crx,  crx0)
    opti.set_initial(cry,  cry0)
    opti.set_initial(crz,  crz0)

    # -----------------------------
    # objective function
    # -----------------------------
    f = 0.0

    # -----------------------------
    # Boundary conditions
    # -----------------------------
    opti.subject_to( eval_curve_0(cbx)  == XB0 )
    opti.subject_to( eval_curve_0(cby)  == YB0 )
    opti.subject_to( eval_curve_0(cbz)  == ZB0 )
    opti.subject_to( eval_curve_0(cpsi) == PSI0 )
    opti.subject_to( eval_curve_1(cbx)  == XB1 )
    opti.subject_to( eval_curve_1(cby)  == YB1 )
    opti.subject_to( eval_curve_1(cbz)  == ZB1 )
    opti.subject_to( eval_curve_1(cpsi) == PSI1 )


    # -----------------------------
    # Smoothing & regularization objective
    # -----------------------------
    if USE_SMOOTHING:
        def smooth_penalty(vec):
            # vec: (Kc,1)
            v1 = vec[0:-2]
            v2 = vec[1:-1]
            v3 = vec[2:  ]
            return ca.sumsqr(v3 - 2*v2 + v1)
        f += (1.0 / n_sample) * (
            smooth_penalty(Xl_c) + smooth_penalty(Yl_c) + smooth_penalty(Zl_c) +
            smooth_penalty(Xr_c) + smooth_penalty(Yr_c) + smooth_penalty(Zr_c) +
            smooth_penalty(Xb_c) + smooth_penalty(Yb_c) + smooth_penalty(Zb_c) + 
            smooth_penalty(Psib_c)
        )
    if USE_REGULARIZATION:
        coef_stack = ca.vertcat(cbx, cby, cbz, cpsi, 
                                clx, cly, clz,
                                crx, cry, crz)
        f += ca.sumsqr(coef_stack)

    # -----------------------------
    # Wheel-base distance constraints
    # -----------------------------
    if USE_WHEEL_BASE_DISTANCE:
        # squared distances (elementwise over collocation points)
        dist2_L = (Xb_c - Xl_c)**2 + (Yb_c - Yl_c)**2 + (Zb_c - Zl_c)**2
        dist2_R = (Xb_c - Xr_c)**2 + (Yb_c - Yr_c)**2 + (Zb_c - Zr_c)**2

        opti.subject_to(dist2_L >= (WHEEL_BASE_MIN**2))
        opti.subject_to(dist2_L <= (WHEEL_BASE_MAX**2))
        opti.subject_to(dist2_R >= (WHEEL_BASE_MIN**2))
        opti.subject_to(dist2_R <= (WHEEL_BASE_MAX**2))

        f += (1.0 / n_sample) * (ca.sumsqr(dist2_L - WHEEL_BASE_NOMINAL**2) +
                                 ca.sumsqr(dist2_R - WHEEL_BASE_NOMINAL**2))

    # -----------------------------
    # Wheels distance constraints
    # -----------------------------    
    if USE_WHEEL_DISTANCE:
        # squared distances (elementwise over collocation points)
        dist2 = (Xr_c - Xl_c)**2 + (Yr_c - Yl_c)**2 + (Zr_c - Zl_c)**2

        opti.subject_to(dist2 >= (WHEELS_MIN**2))
        opti.subject_to(dist2 <= (WHEELS_MAX**2))

        f += (1.0 / n_sample) * ca.sumsqr(dist2 - WHEELS_NOMINAL**2)

    # -----------------------------
    # Wheels terrain distance objective
    # -----------------------------    
    H_l_c = terrain_instance.terrain_height_casadi(Xl_c, Yl_c)
    H_r_c = terrain_instance.terrain_height_casadi(Xr_c, Yr_c)

    dist_terrain_l = Zl_c - H_l_c
    dist_terrain_r = Zr_c - H_r_c
    if USE_WHEEL_TERRAIN_DISTANCE:
        f += (1.0 / n_sample) * (ca.sumsqr(dist_terrain_l) +
                                 ca.sumsqr(dist_terrain_r))

    # -----------------------------
    # Anti-crossing objective
    # 首先对齐base的偏航角和路径切线方向，这时可以获得base坐标系。
    # 防止轮子轨迹交叉可以通过base坐标系下左轮y坐标始终大于右轮y坐标来实现。
    # -----------------------------
    if USE_ANTI_CROSS:
        # -----------------------------
        # yaw-path alignment
        # -----------------------------
        Xb_next = Xb_c[1:]
        Xb_prev = Xb_c[:-1]
        Yb_next = Yb_c[1:]
        Yb_prev = Yb_c[:-1]

        dx = Xb_next - Xb_prev
        dy = Yb_next - Yb_prev

        Psi_next = Psib_c[1:]
        Psi_prev = Psib_c[:-1]
        Psi_mid  = 0.5 * (Psi_next + Psi_prev)

        cos_psi = ca.cos(Psi_mid)
        sin_psi = ca.sin(Psi_mid)

        # evaluated yaw from base path tangent
        eps_norm = 1e-6
        norm_t   = ca.sqrt(dx*dx + dy*dy + eps_norm)
        tx = dx / norm_t
        ty = dy / norm_t
        # fitted yaw direction
        fx = cos_psi
        fy = sin_psi

        cos_align = tx*fx + ty*fy
        yaw_misalign = 1.0 - cos_align   # the smaller the better
        f += (W_YAW_ALIGN / (n_sample)) * ca.sumsqr(yaw_misalign)

        # -----------------------------
        # Left–right ordering in base frame
        # -----------------------------
        cos_psi_full = ca.cos(Psib_c)
        sin_psi_full = ca.sin(Psib_c)
        u_lat_x = -sin_psi_full 
        u_lat_y =  cos_psi_full

        vLx = Xl_c - Xb_c
        vLy = Yl_c - Yb_c
        vRx = Xr_c - Xb_c
        vRy = Yr_c - Yb_c

        yL_base = vLx * u_lat_x + vLy * u_lat_y
        yR_base = vRx * u_lat_x + vRy * u_lat_y

        opti.subject_to( yL_base >= yR_base )  # hard constraint

    # -----------------------------
    # Quasi-static force & moment balance (soft)
    # -----------------------------
    if USE_QUASI_STATIC:
        # contact force estimates from penetration
        pen_l = smooth_pos(-dist_terrain_l)
        pen_r = smooth_pos(-dist_terrain_r)
        fLn_spring = K_CONTACT * pen_l
        fRn_spring = K_CONTACT * pen_r

        fLn  = fLx * nLx + fLy * nLy + fLz * nLz  # 法向分量
        fRn  = fRx * nRx + fRy * nRy + fRz * nRz
        f += (W_CONTACT_STIFF / n_sample) * (
            ca.sumsqr(fLn - fLn_spring) +
            ca.sumsqr(fRn - fRn_spring)
        )

        # -----------------------------
        # 合力：fL + fR + gravity ≈ 0
        # -----------------------------
        Fx_res = fLx + fRx           # x 方向合力（理想为 0）
        Fy_res = fLy + fRy           # y 方向合力
        Fz_res = fLz + fRz - MASS * GRAVITY  # z 方向：fLz + fRz ≈ mg

        f += (W_FORCE_BAL / n_sample) * (
            ca.sumsqr(Fx_res) + ca.sumsqr(Fy_res) + ca.sumsqr(Fz_res)
        )

        # -----------------------------
        # 合矩：tau = r_L × f_L + r_R × f_R ≈ 0
        # -----------------------------
        rLx = Xl_c - Xb_c
        rLy = Yl_c - Yb_c
        rLz = Zl_c - Zb_c

        rRx = Xr_c - Xb_c
        rRy = Yr_c - Yb_c
        rRz = Zr_c - Zb_c

        # r × f
        tauLx = rLy * fLz - rLz * fLy
        tauLy = rLz * fLx - rLx * fLz
        tauLz = rLx * fLy - rLy * fLx

        tauRx = rRy * fRz - rRz * fRy
        tauRy = rRz * fRx - rRx * fRz
        tauRz = rRx * fRy - rRy * fRx

        tau_x = tauLx + tauRx
        tau_y = tauLy + tauRy
        tau_z = tauLz + tauRz

        f += (W_TORQUE_BAL / n_sample) * (
            ca.sumsqr(tau_x) + ca.sumsqr(tau_y) + ca.sumsqr(tau_z)
        )

    # -----------------------------
    # friction constraints
    # -----------------------------
    if USE_FRICTION_CONE:
        # 左轮
        fLn = fLx * nLx + fLy * nLy + fLz * nLz          # 法向分量
        fLtx = fLx - fLn * nLx                            # 切向分量
        fLty = fLy - fLn * nLy
        fLtz = fLz - fLn * nLz
        fLt_sq = fLtx*fLtx + fLty*fLty + fLtz*fLtz        # ||f_t||^2

        # 右轮
        fRn = fRx * nRx + fRy * nRy + fRz * nRz
        fRtx = fRx - fRn * nRx
        fRty = fRy - fRn * nRy
        fRtz = fRz - fRn * nRz
        fRt_sq = fRtx*fRtx + fRty*fRty + fRtz*fRtz

        mu2 = MU_FRICTION**2

        # 软摩擦锥：viol = max(0, ||f_t||^2 - mu^2 * f_n^2)
        viol_fric_L = smooth_pos(fLt_sq - mu2 * (fLn**2))
        viol_fric_R = smooth_pos(fRt_sq - mu2 * (fRn**2))

        # 不许拉力：f_n >= 0  （软约束）
        viol_normal_L = smooth_pos(-fLn)
        viol_normal_R = smooth_pos(-fRn)

        # 最大拉力： f_l <= mg + ma  （软约束）
        # a = 1.0 m/s²
        fN_max = MASS * (GRAVITY + 1.0)
        # viol_normal_L_max = smooth_pos(fLn - fN_max)
        # viol_normal_R_max = smooth_pos(fRn - fN_max)
        # f += (W_FRIC_CONE / n_sample) * (
        #     ca.sumsqr(viol_normal_L_max) + ca.sumsqr(viol_normal_R_max)
        # )
        opti.subject_to( fLn <= fN_max )
        opti.subject_to( fRn <= fN_max )
        f += (W_FRIC_CONE / n_sample) * (
            ca.sumsqr(viol_fric_L) + ca.sumsqr(viol_fric_R) +
            ca.sumsqr(viol_normal_L) + ca.sumsqr(viol_normal_R)
        )

    # -----------------------------
    # minimize objective
    # -----------------------------
    opti.minimize(f)

    # -----------------------------
    # Solve
    # -----------------------------
    print("Solving NLP...")
    success = False
    try:
        sol = opti.solve()
        success = True
    except RuntimeError as e:
        print("Solver failed, using debug solution. Error:", e)
        sol = opti.debug
    if success:
        cbx_opt  = sol.value(cbx)
        cby_opt  = sol.value(cby)
        cbz_opt  = sol.value(cbz)
        cpsi_opt = sol.value(cpsi)
        clx_opt  = sol.value(clx)
        cly_opt  = sol.value(cly)
        clz_opt  = sol.value(clz)
        crx_opt  = sol.value(crx)
        cry_opt  = sol.value(cry)
        crz_opt  = sol.value(crz)
        if VISUAL_OPT:
            print("-----------------------------")
            print("Optimized polynomial coefficients:")
            print("cbx0:", cbx_opt)
            print("cby0:", cby_opt)
            print("cbz0:", cbz_opt)
            print("cpsi0:", cpsi_opt)
            print("clx0:", clx_opt)
            print("cly0:", cly_opt)
            print("clz0:", clz_opt)
            print("crx0:", crx_opt)
            print("cry0:", cry_opt)
            print("crz0:", crz_opt)
            print("-----------------------------")
            print()

            Sd_vis = np.linspace(0.0, 1.0, 201)
            Phi_vis = phi_poly(Sd_vis, N_COEFF)
            Xb_d = Phi_vis @ cbx_opt
            Yb_d = Phi_vis @ cby_opt
            Zb_d = Phi_vis @ cbz_opt
            Psi_d = Phi_vis @ cpsi_opt
            Xl_d = Phi_vis @ clx_opt
            Yl_d = Phi_vis @ cly_opt
            Zl_d = Phi_vis @ clz_opt
            Xr_d = Phi_vis @ crx_opt
            Yr_d = Phi_vis @ cry_opt
            Zr_d = Phi_vis @ crz_opt

            traj_opt = ((Xb_d, Yb_d, Zb_d, Psi_d),
                        (Xl_d, Yl_d, Zl_d),
                        (Xr_d, Yr_d, Zr_d))

            plot_on_terrain(traj_opt, terrain=terrain_instance)

if __name__ == "__main__":
    # read settings from json
    file_name = "/home/yizhe/trajopt/result/RBF/1753771011_476223707.json"
    timestamp, pos, ori, mean_weight, cur_rbf_grid, rbf_weight = read_from_json(file_name)

    RBF_CENTERS = cur_rbf_grid        # (N,2) numpy
    RBF_WEIGHTS = rbf_weight.ravel()  # (N,)
    RBF_SIGMA = 0.14

    build_and_solve()