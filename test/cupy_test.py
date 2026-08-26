import numpy as np, cupy as cp, time

N, M = 2000, 20000
XY = np.random.randn(N,2)
C  = np.random.randn(M,2)
w  = np.random.randn(M)
l  = 0.3

def rbf_cpu(XY, C, w, l):
    XYb = XY[:,None,:]
    Cb  = C[None,:,:]
    D   = XYb - Cb
    r2  = np.sum(D*D, axis=2)
    phi = np.exp(-r2/(2*l*l))
    h   = phi @ w
    gx  = ((-phi * D[:,:,0] / (l*l)) @ w)
    gy  = ((-phi * D[:,:,1] / (l*l)) @ w)
    return h,gx,gy

def rbf_gpu(XY, C, w, l):
    XYg = cp.asarray(XY, dtype=cp.float32)
    Cg  = cp.asarray(C, dtype=cp.float32)
    wg  = cp.asarray(w, dtype=cp.float32)
    XYb = XYg[:,None,:]; Cb = Cg[None,:,:]
    D   = XYb - Cb
    r2  = cp.sum(D*D, axis=2)
    phi = cp.exp(-r2/(2*l*l))
    h   = phi @ wg
    gx  = ((-phi * D[:,:,0] / (l*l)) @ wg)
    gy  = ((-phi * D[:,:,1] / (l*l)) @ wg)
    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(h), cp.asnumpy(gx), cp.asnumpy(gy)

# 1) 先把逐元素链式操作融合成单个 kernel —— 计算 phi
@cp.fuse()
def fused_phi(XY0, XY1, C0, C1, inv_two_l2):
    # D0 = XY[:,None,0] - C[None,:,0]
    # D1 = XY[:,None,1] - C[None,:,1]
    D0 = XY0[:, None] - C0[None, :]
    D1 = XY1[:, None] - C1[None, :]
    # r2 = D0^2 + D1^2; phi = exp(-r2 * inv_two_l2)
    return cp.exp(-(D0*D0 + D1*D1) * inv_two_l2)

# 2) 同样融合 “-phi * D0 * inv_l2” 这类逐元素式子，供 gx / gy 使用
@cp.fuse()
def fused_phi_dx(phi, XY0, C0, inv_l2):
    D0 = XY0[:, None] - C0[None, :]
    return -phi * D0 * inv_l2

@cp.fuse()
def fused_phi_dy(phi, XY1, C1, inv_l2):
    D1 = XY1[:, None] - C1[None, :]
    return -phi * D1 * inv_l2

def rbf_gpu_fused(XY, C, w, l):
    # 全部用 float32（对 RTX 3070 非常关键）
    XYg = cp.asarray(XY, dtype=cp.float32)   # (N,2)
    Cg  = cp.asarray(C,  dtype=cp.float32)   # (M,2)
    wg  = cp.asarray(w,  dtype=cp.float32)   # (M,)

    l   = cp.float32(l)
    inv_two_l2 = cp.float32(1.0) / (cp.float32(2.0) * l * l)
    inv_l2     = cp.float32(1.0) / (l * l)

    # 逐元素融合核：一次完成 (减 -> 平方 -> 求和 -> exp)
    phi = fused_phi(XYg[:, 0], XYg[:, 1], Cg[:, 0], Cg[:, 1], inv_two_l2)  # (N,M)

    # h = phi @ w  （交给 cuBLAS 的 GEMV，很快）
    h  = phi @ wg  # (N,)

    # 为 gx / gy 准备逐元素融合核（-phi * D / l^2），然后再 @ w
    phi_dx = fused_phi_dx(phi, XYg[:, 0], Cg[:, 0], inv_l2)   # (N,M)
    phi_dy = fused_phi_dy(phi, XYg[:, 1], Cg[:, 1], inv_l2)   # (N,M)

    gx = phi_dx @ wg
    gy = phi_dy @ wg

    cp.cuda.Stream.null.synchronize()   # 同步测时/取数前
    return cp.asnumpy(h), cp.asnumpy(gx), cp.asnumpy(gy)

t0=time.time(); _=rbf_cpu(XY,C,w,l); t1=time.time()
_ = rbf_gpu(XY,C,w,l); t2=time.time()
_ = rbf_gpu_fused(XY,C,w,l); t3=time.time()
print(f"CPU: {t1-t0:.3f}s, GPU: {t2-t1:.3f}s, GPU Fused: {t3-t2:.3f}s")
