#ifndef OBJECTIVE_BASE_HPP
#define OBJECTIVE_BASE_HPP

/**
 * 目标函数项基类
 * 所有目标函数项都应继承此类
 */
class ObjectiveBase {
public:
    virtual ~ObjectiveBase() = default;
    
    /**
     * 评估目标项值
     * @param x 优化变量
     * @return 目标项值
     */
    virtual double eval_f(const double* x) = 0;
    
    /**
     * 评估目标项梯度
     * @param x 优化变量
     * @param grad_f 输出的梯度（需要累加，不是覆盖）
     */
    virtual void eval_grad_f(const double* x, double* grad_f) = 0;
    
    /**
     * 是否启用此目标项
     */
    virtual bool isEnabled() const = 0;
};

#endif // OBJECTIVE_BASE_HPP

