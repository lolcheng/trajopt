#ifndef INITIAL_GUESS_PROXIMITY_OBJECTIVE_HPP
#define INITIAL_GUESS_PROXIMITY_OBJECTIVE_HPP

#include "objective_base.hpp"
#include <vector>

class InitialGuessProximityObjective : public ObjectiveBase {
public:
    InitialGuessProximityObjective(const std::vector<double>& x_ref, double weight)
        : x_ref_(x_ref),
          weight_(weight),
          enabled_(true) {}

    double eval_f(const double* x) override;
    void eval_grad_f(const double* x, double* grad_f) override;
    bool isEnabled() const override { return enabled_; }
    void setEnabled(bool enabled) { enabled_ = enabled; }

private:
    std::vector<double> x_ref_;
    double weight_;
    bool enabled_;
};

#endif
