#include "main.h"

int main() {
    constexpr double X_prec = 0.8;
    std::vector<double> X{0.2, 0.6, 1.0, 1.4};

    print_nodes(function, X);

    const double lagrange_val = Lagrange_polynomial(function, X, X_prec);

    const double newton_val = Newton_polynomial(function, X, X_prec);

    const double exact = function(X_prec);
    std::cout << std::format("Точное значение f(x*) = ln({:.1f}) = {:.8f}\n", X_prec, exact);

    std::cout << std::format("Абсолютная погрешность (Лагранж): {:.3e}\n",
                             std::abs(lagrange_val - exact));
    std::cout << std::format("Абсолютная погрешность (Ньютон):  {:.3e}\n",
                             std::abs(newton_val - exact));

    std::cout << std::format("Разность |Лагранж - Ньютон|: {:.3e}\n",
                             std::abs(lagrange_val - newton_val));

    return 0;
}