#include "main.h"
#include <iomanip>

double function(const double x) {
    return std::log(x);
}


void print_nodes(const func y, const std::vector<double>& x_points) {
    std::cout << "Исходные узлы интерполяции:\n";
    std::cout << std::format("{:>8} {:>16}\n", "x_i", "y_i = ln(x_i)");
    for (size_t i = 0; i < x_points.size(); ++i) {
        std::cout << std::format("{:8.4f} {:12.6f}\n", x_points[i], y(x_points[i]));
    }
    std::cout << std::endl;
}


double Lagrange_polynomial(const func y, const std::vector<double>& x_points,
                           const double x_precision) {
    const int n = x_points.size();
    std::vector<double> coefficients(n);

    for (int i = 0; i < n; ++i) {
        double w = 1;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            w *= (x_points[i] - x_points[j]);
        }
        coefficients[i] = y(x_points[i]) / w;
    }

    double approximation = 0;
    std::string out{};

    for (int i = 0; i < n; ++i) {
        double prod = 1;
        std::string x_out{};
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            prod *= (x_precision - x_points[j]);
            x_out += std::format("(x - {:.1f})", x_points[j]);
        }
        approximation += coefficients[i] * prod;

        if (coefficients[i] == 0)
            continue;
        if (coefficients[i] < 0)
            out += " - ";
        else if (i > 0)
            out += " + ";
        out += std::format("{:.5f}{}", std::abs(coefficients[i]), x_out);
    }

    std::cout << "Многочлен Лагранжа:\n" << out << std::endl;
    std::cout << std::format("Значение в x* = {:.1f}: {:.8f}\n", x_precision, approximation);
    std::cout << std::endl;
    return approximation;
}


double Newton_polynomial(const func y, const std::vector<double>& x_points,
                         const double x_precision) {
    const int n = x_points.size();

    // двумерная таблица разделённых разностей
    std::vector<std::vector<double>> table(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        table[i][0] = y(x_points[i]);

    for (int k = 1; k < n; ++k) {
        for (int i = 0; i + k < n; ++i) {
            table[i][k] = (table[i + 1][k - 1] - table[i][k - 1])
                          / (x_points[i + k] - x_points[i]);
        }
    }

    std::cout << "Таблица разделённых разностей:\n";
    std::cout << std::format("{:>8}", "x_i");
    std::cout << std::format("{:>14}", "f(x_i)");
    for (int k = 0; k < n - 1; ++k)
        std::cout << std::format("{:>14}",  std::to_string(k));
    std::cout << "\n";
    for (int i = 0; i < n; ++i) {
        std::cout << std::format("{:8.4f}", x_points[i]);
        for (int k = 0; k < n - i; ++k)
            std::cout << std::format("{:14.6f}", table[i][k]);
        std::cout << "\n";
    }
    std::cout << std::endl;


    double approximation = 0;
    std::string out{};

    for (int i = 0; i < n; ++i) {
        double prod = 1;
        std::string x_out{};
        for (int j = 0; j < i; ++j) {
            prod *= (x_precision - x_points[j]);
            x_out += std::format("(x - {:.1f})", x_points[j]);
        }
        approximation += table[0][i] * prod;

        if (table[0][i] == 0)
            continue;
        if (table[0][i] < 0)
            out += " - ";
        else if (i > 0)
            out += " + ";
        out += std::format("{:.5f}{}", std::abs(table[0][i]), x_out);
    }

    std::cout << "Многочлен Ньютона:\n" << out << std::endl;
    std::cout << std::format("Значение в x* = {:.1f}: {:.8f}\n", x_precision, approximation);
    std::cout << std::endl;
    return approximation;
}