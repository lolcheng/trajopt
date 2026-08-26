#ifndef BASIS_FUNCTION_HPP
#define BASIS_FUNCTION_HPP

/**
 * 基函数接口类
 * 所有基函数（多项式基、切比雪夫基等）都应继承此类
 */
class BasisFunction {
public:
    virtual ~BasisFunction() = default;
    
    /**
     * 构造基矩阵
     * 计算在给定采样点上的基函数值
     * 
     * @param s 采样点数组 (n_samples)，通常 s ∈ [0, 1]
     * @param n_samples 采样点数量
     * @param n_coeff 基函数数量（系数数量）
     * @param Phi 输出的基矩阵 (n_samples × n_coeff)，按行优先存储
     *            Phi[i, j] = phi_j(s[i])，其中 phi_j 是第 j 个基函数
     */
    virtual void computeBasisMatrix(const double* s, int n_samples, int n_coeff,
                                    double* Phi) = 0;
    
    /**
     * 评估轨迹在单个点的值: x(s) = sum(coeff[j] * phi_j(s))
     * 
     * @param s 采样点，通常 s ∈ [0, 1]
     * @param coeff 系数数组 (n_coeff)
     * @param n_coeff 系数数量
     * @return 轨迹值 x(s)
     */
    virtual double evalPoint(double s, const double* coeff, int n_coeff) = 0;
    
    /**
     * 评估轨迹: x(s) = Phi * coeff
     * 
     * @param Phi 基矩阵 (n_samples × n_coeff)，按行优先存储
     * @param n_samples 采样点数量
     * @param coeff 系数数组 (n_coeff)
     * @param n_coeff 系数数量
     * @param x 输出的轨迹值 (n_samples)
     */
    virtual void eval(const double* Phi, int n_samples, const double* coeff,
                     int n_coeff, double* x) = 0;
    
    /**
     * 评估轨迹导数: dx/ds = sum(coeff[j] * dphi_j(s)/ds)
     * 
     * @param s 采样点数组 (n_samples)
     * @param n_samples 采样点数量
     * @param coeff 系数数组 (n_coeff)
     * @param n_coeff 系数数量
     * @param dx 输出的导数值 (n_samples)
     */
    virtual void evalDerivative(const double* s, int n_samples,
                               const double* coeff, int n_coeff, double* dx) = 0;
    
    /**
     * 评估轨迹导数在单个点的值: dx/ds(s) = sum(coeff[j] * dphi_j(s)/ds)
     * 
     * @param s 采样点
     * @param coeff 系数数组 (n_coeff)
     * @param n_coeff 系数数量
     * @return 导数值 dx/ds(s)
     */
    virtual double evalDerivativePoint(double s, const double* coeff, int n_coeff) = 0;

    /**
     * 在单个参数 s 处，填充各基函数对 s 的一阶、二阶导数行向量：
     *   dphi_ds[j]   = d phi_j(s) / ds
     *   d2phi_ds2[j] = d^2 phi_j(s) / ds^2
     * 用于端点速度/加速度（对路径参数）等线性约束的雅可比。
     */
    virtual void fillDerivativeBasisAt(double s, int n_coeff,
                                       double* dphi_ds, double* d2phi_ds2) = 0;
};

#endif // BASIS_FUNCTION_HPP

