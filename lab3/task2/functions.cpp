#include "main.h"
#include "../../lab1/task2/main.h"

double find_h(std::vector<double> x, int index) {
    return x[index] - x[index - 1];
}

double* cubic_spline_equation(const std::vector<double> &x, const std::vector<double> &f) {
    const int n = x.size() - 1;

    Matrix<double> tridiag(n - 1, n - 1);
    std::vector<double> d_eq(n - 1);

    { // i = 2 методичка
        double h = find_h(x, 2); // h_2 = x_2 - x_1
        double prev_h = find_h(x, 1); // h_1 = x_1 - x_0
        tridiag.coefficients[0][0] = 2 * (prev_h + h); // 2(h_1 + h_2)
        tridiag.coefficients[0][1] = h; // h_2
        d_eq[0] = 3 * ((f[2] - f[1]) / h - (f[1] - f[0]) / prev_h);
    }

    for (int i = 3; i < n; ++i) { // i = 3,...,n-1 методичка
        double h = find_h(x, i);
        double prev_h = find_h(x, i - 1);
        int index = i - 2; // сдвиг т.к. с_1 в solution[0], с_2 в solution[1]
        tridiag.coefficients[index][index - 1] = prev_h; // h_{i-1}
        tridiag.coefficients[index][index] = 2 * (prev_h + h); // 2(h_{i-1}+h_i)
        tridiag.coefficients[index][index + 1] = h; // h_i
        d_eq[index] = 3 * ((f[i] - f[i - 1]) / h - (f[i - 1] - f[i - 2]) / prev_h);
    }

    {  // i = n методичка
        double h = find_h(x, n);
        double prev_h = find_h(x, n - 1);
        int index = n - 2;
        tridiag.coefficients[index][index - 1] = prev_h;
        tridiag.coefficients[index][index] = 2 * (prev_h + h);
        d_eq[index] = 3 * ((f[n] - f[n - 1]) / h - (f[n - 1] - f[n - 2]) / prev_h);
    }

    return thomas_algorithm(tridiag, d_eq.data());
}

// для одного участка i
SplineCoeffs compute_spline_coeffs(const std::vector<double> &x,
                                   const std::vector<double> &f,
                                   double* solution) {
    const int n = x.size() - 1;
    SplineCoeffs s;
    s.a.resize(n); s.b.resize(n); s.c.resize(n); s.d.resize(n);

    s.c[0] = 0;
    for (int i = 0; i < n; ++i) {
        s.a[i] = f[i];
        double h = find_h(x, i + 1);
        s.c[i + 1] = solution[i];
        s.b[i] = (f[i + 1] - f[i]) / h - h / 3 * (s.c[i + 1] + 2 * s.c[i]);
        s.d[i] = (s.c[i + 1] - s.c[i]) / (3 * h);
    }
    return s;
}

double eval_spline(const SplineCoeffs& s, const std::vector<double>& x, int i, double px) {
    double dx = px - x[i];
    return s.a[i] + s.b[i] * dx + s.c[i] * dx * dx + s.d[i] * dx * dx * dx;
}

double cubic_spline(const std::vector<double>& x, const std::vector<double>& f, const double point_x) {
    const int n = x.size() - 1;

    double* solution = cubic_spline_equation(x, f);
    SplineCoeffs s = compute_spline_coeffs(x, f, solution);
    delete[] solution;


    std::cout << "  Исходные узлы и значения функции \n";
    for (int i = 0; i <= n; ++i)
        std::cout << std::format("  x[{}] = {:.6f}, f[{}] = {:.6f}\n", i, x[i], i, f[i]);

    std::cout << "\n  Коэффициенты и интервалы для каждого участка \n";
    for (int i = 0; i < n; ++i) {
        std::cout << std::format(
                "  [{:.4f}, {:.4f}]:  S{}(x) = {:.6f} + {:.6f}(x-{:.4f}) + {:.6f}(x-{:.4f})^2 + {:.6f}(x-{:.4f})^3\n",
                x[i], x[i + 1], i, s.a[i], s.b[i], x[i], s.c[i], x[i], s.d[i], x[i]);
    }

    // поиск интервала для x*
    int idx = -1;
    for (int i = 0; i < n; ++i) {
        if (point_x >= x[i] && point_x <= x[i + 1]) {
            idx = i; break;
        }
    }
    double result = 0.0;
    if (idx < 0) {
        std::cout << "\n x* вне диапазона интерполяции\n";
        idx = (point_x < x[0]) ? 0 : n - 1;
    }
    result = eval_spline(s, x, idx, point_x);


    std::cout << std::format("\n  Значение сплайна в x* \n");
    std::cout << std::format("  x* = {:.6f}\n", point_x);
    std::cout << std::format("  Номер интервала: {}\n", idx + 1);
    std::cout << std::format("  S(x*) = {:.6f}\n", result);


    std::cout << "\n  Проверка прохождения через узлы \n";
    for (int i = 0; i <= n; ++i) {
        double val = (i == n) ? eval_spline(s, x, n - 1, x[n])
                              : eval_spline(s, x, i, x[i]);
        std::cout << std::format("  S({:.4f}) = {:.6f}, f = {:.6f}, |Δ| = {:.2e}\n",
                                 x[i], val, f[i], std::abs(val - f[i]));
    }

    //  вторая производная на концах
    // S"(x) = 2c_i + 6d_i (x - x_i)
    double left2  = 2 * s.c[0];                     // в x[0]
    double right2 = 2 * s.c[n - 1] + 6 * s.d[n - 1] * (x[n] - x[n - 1]); // в x[n]
    std::cout << "\n Вторая производная на концах \n";
    std::cout << std::format("  S''({:.4f}) = {:.6f}\n", x[0], left2);
    std::cout << std::format("  S''({:.4f}) = {:.6f}\n", x[n], right2);
    std::cout << "  (граничные условия: S''(x0)=0, S''(xn)=0 — проверка)\n";

    return result;
}