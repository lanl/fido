#include <cmath>
#include <fmt/core.h>
#include <limits>
#include <nlopt.hpp>
#include <vector>
#include <functional>

double my_func(unsigned n, const double* x, double* grad, void* my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / std::sqrt(x[1]);
    }
    return std::sqrt(x[1]);
}

struct my_constraint_data {
    double a, b;
};

double my_constraint(unsigned n, const double* x, double* grad, void* data)
{
    auto [a, b] = *reinterpret_cast<my_constraint_data*>(data);

    if (grad) {
        grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
        grad[1] = -1.0;
    }
    return (a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1];
}

template <typename T>
double wrap_lambda(unsigned n, const double* x, double* grad, void* data)
{
    T* l = reinterpret_cast<T*>(data);
    return std::invoke(*l, n, x, grad);
}

int main()
{
    nlopt::opt opt(nlopt::LD_MMA, 2);

    std::vector<double> lower_bounds{std::numeric_limits<double>::lowest(), 0};

    opt.set_lower_bounds(lower_bounds);

    int count {};
    auto my_l_func = [&count](unsigned n, const double* x, double* grad) mutable {
        ++count;
        if (grad) {
            grad[0] = 0.0;
            grad[1] = 0.5 / std::sqrt(x[1]);
        }
        return std::sqrt(x[1]);
    };
    
    opt.set_min_objective(wrap_lambda<decltype(my_l_func)>, &my_l_func);

    auto cd = std::vector<my_constraint_data>{{2.0, 0}, {-1, 1}};
    opt.add_inequality_constraint(my_constraint, &cd[0], 1e-8);
    opt.add_inequality_constraint(my_constraint, &cd[1], 1e-8);

    opt.set_xtol_rel(1e-4);

    auto x = std::vector{1.234, 5.678};
    double minf;

    auto result = opt.optimize(x, minf);
    fmt::print("found minimum in {} evaluations at f({}, {}) = {}\n", count, x[0], x[1], minf);
}