#include "main.h"
#include <cmath>
#include <format>
#include <iostream>
#include <string>
#include <vector>

// MSE по коэффициентам
double compute_mse(const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<double>& coeffs) {
    double mse = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double f = 0.0, xc = 1.0;
        for (size_t j = 0; j < coeffs.size(); ++j) {
            f += coeffs[j] * xc;
            xc *= x[i];
        }
        mse += (f - y[i]) * (f - y[i]);
    }
    return mse / x.size();
}

void print_polynomial(const std::vector<double>& c, const std::string& name) {
    std::cout << name << "(x) = ";
    for (size_t i = 0; i < c.size(); ++i) {
        if (i == 0) {
            std::cout << std::format("{:.4f}", c[i]);
        } else {
            std::cout << (c[i] >= 0 ? " + " : " - ")
                      << std::format("{:.4f}", std::abs(c[i])) << "x";
            if (i > 1) std::cout << "^" << i;
        }
    }
    std::cout << "\n";
}

int main() {
    const std::vector<double> x{0.1, 0.5, 0.9, 1.3, 1.7, 2.1};
    const std::vector<double> y{-2.3026, -0.69315, -0.10536, 0.26236, 0.53063, 0.74194};

    std::cout << " Многочлен 1-й степени\n";
    auto c1 = approximation_polynomial_first_degree(x, y);
    print_polynomial(c1, "P1");
    double mse1 = compute_mse(x, y, c1);
    std::cout << std::format("MSE1 = {:.6f}\n\n", mse1);

    std::cout << " Многочлен 2-й степени\n";
    auto c2 = approximation_polynomial_second_degree(x, y);
    print_polynomial(c2, "P2");
    double mse2 = compute_mse(x, y, c2);
    std::cout << std::format("MSE2 = {:.6f}\n\n", mse2);

    std::cout << " Многочлен 3-й степени\n";
    auto c3 = approximation_polynomial_third_degree(x, y);
    print_polynomial(c3, "P3");
    double mse3 = compute_mse(x, y, c3);
    std::cout << std::format("MSE3 = {:.6f}\n\n", mse3);

    std::cout << " Многочлен 5-й степени\n";
    auto c5 = approximation_polynomial_fifth_degree(x, y);
    print_polynomial(c5, "P5");
    double mse5 = compute_mse(x, y, c5);
    std::cout << std::format("MSE5 = {:.6f}\n\n", mse5);


    std::cout << std::format("MSE1 = {:.6f}\n", mse1);
    std::cout << std::format("MSE2 = {:.6f}\n", mse2);
    std::cout << std::format("MSE3 = {:.6f}\n", mse3);
    std::cout << std::format("MSE5 = {:.6f}\n", mse5);

    /*double best = std::min({mse1, mse2, mse3, mse5});
    if (best == mse1)      std::cout << "Лучшее приближение: P1(x)\n";
    else if (best == mse2) std::cout << "Лучшее приближение: P2(x)\n";
    else if (best == mse3) std::cout << "Лучшее приближение: P3(x)\n";
    else*/                   std::cout << "Сравнивать не корректно с полиномом 5-й степени из-за "
                                        "плохой обусловленности системы 6x6\n";

    return 0;
}