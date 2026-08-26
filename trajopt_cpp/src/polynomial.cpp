#include "polynomial.hpp"

// 静态成员定义
std::shared_ptr<BasisFunction> Polynomial::basis_ = std::make_shared<ChebyshevBasis>();

