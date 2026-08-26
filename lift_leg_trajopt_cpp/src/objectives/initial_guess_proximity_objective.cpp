#include "objectives/initial_guess_proximity_objective.hpp"

double InitialGuessProximityObjective::eval_f(const double* x) {
    if (!enabled_) return 0.0;
    double obj = 0.0;
    const int n = static_cast<int>(x_ref_.size());
    for (int i = 0; i < n; ++i) {
        const double d = x[i] - x_ref_[i];
        obj += d * d;
    }
    return weight_ * obj;
}

void InitialGuessProximityObjective::eval_grad_f(const double* x, double* grad_f) {
    if (!enabled_) return;
    const int n = static_cast<int>(x_ref_.size());
    for (int i = 0; i < n; ++i) {
        grad_f[i] += 2.0 * weight_ * (x[i] - x_ref_[i]);
    }
}
