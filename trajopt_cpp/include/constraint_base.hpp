#ifndef CONSTRAINT_BASE_HPP
#define CONSTRAINT_BASE_HPP

/**
 * 约束基类
 * 所有约束都继承此类，实现统一的接口
 */
class ConstraintBase {
public:
    virtual ~ConstraintBase() = default;
    
    /**
     * 返回约束数量
     */
    virtual int getNumConstraints() const = 0;
    
    /**
     * 评估约束值 g(x)
     * @param x 优化变量
     * @param g 输出的约束值数组
     */
    virtual void eval_g(const double* x, double* g) = 0;
    
    /**
     * 评估约束雅可比 J = ∂g/∂x
     * @param x 优化变量
     * @param jac_g 输出的雅可比矩阵（稀疏格式，按行优先）
     * @param offset 约束索引偏移（用于多个约束组合）
     */
    virtual void eval_jac_g(const double* x, double* jac_g, int offset = 0) = 0;
    
    /**
     * 获取约束边界 [g_lb, g_ub]
     * @param g_lb 约束下界
     * @param g_ub 约束上界
     * @param offset 约束索引偏移
     */
    virtual void getBounds(double* g_lb, double* g_ub, int offset = 0) = 0;
    
    /**
     * 获取雅可比稀疏结构
     * @param iRow 行索引数组
     * @param jCol 列索引数组
     * @param offset 约束索引偏移
     * @return 非零元素数量
     */
    virtual int getJacobianStructure(int* iRow, int* jCol, int offset = 0) = 0;
    
    /**
     * 是否启用此约束
     */
    virtual bool isEnabled() const { return true; }
};

#endif // CONSTRAINT_BASE_HPP

