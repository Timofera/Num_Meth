#include "main.h"

std::vector<double> approximation_polynomial_first_degree(const std::vector<double>& x, const std::vector<double>& y) {
    const int N = x.size();

    LinearEquation<double> equation(2);

    double** coef = equation.matrix.coefficients;

    double xsum = 0, ysum = 0, xxsum = 0, xysum = 0;

    for (int i = 0; i < N; ++i) {
        xsum += x[i];
        xxsum += x[i] * x[i];
        ysum += y[i];
        xysum += x[i] * y[i];
    }

    coef[0][0] = N;
    coef[0][1] = xsum;
    coef[1][0] = xsum;
    coef[1][1] = xxsum;

    equation.b[0] = ysum;
    equation.b[1] = xysum;

    const double* a = solveEquation(equation);


    std::vector<double> result(2);
    result[0] = a[0];
    result[1] = a[1];

    delete[] a;

    double square_error = 0;

    for (int i = 0; i < N; ++i) {
        double Fsum = 0, x_coef = 1;

        for (int j = 0; j < 2; ++j) {
            Fsum += result[j] * x_coef; // a_j * x^i
            x_coef *= x[i]; // x^(j+1)
        }

        double diff = Fsum - y[i]; // F_n(x_i) - y_i

        square_error += diff * diff;
    }

    std::string out;

    for (int i = 0; i < 2; ++i) {

        if (result[i] < 0) {
            out += '-';
            if (i != 0)
                out += ' ';
        } else if (i != 0) {
            out += "+ ";
        }

        out += std::format("{:.4f}", std::abs(result[i]));

        if (i > 0) {
            out += "x";
            if (i > 1) {
                out += "^" + std::to_string(i);
            }
        }

        out += ' ';
    }

    std::cout << out << std::endl;
    std::cout << std::format("Square Error: {:.4f}", square_error) << std::endl;

    return result;
}

std::vector<double> approximation_polynomial_second_degree(const std::vector<double>& x, const std::vector<double>& y) {
    const int N = x.size();

    LinearEquation<double> equation(3);

    double** coef = equation.matrix.coefficients;

    std::vector<double> xsums(4), ysums(3);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 4; ++j) {
            xsums[j] += std::pow(x[i], j + 1); // xsums[j] = sum x^(j+1)
        }
        for (int j = 0; j < 3; ++j) {
            ysums[j] += y[i] * std::pow(x[i], j);
        }
    }

    coef[0][0] = N; // sum x^0
    coef[0][1] = xsums[0]; // sum x^1
    coef[0][2] = xsums[1];

    coef[1][0] = xsums[0]; // sum x^1
    coef[1][1] = xsums[1]; // sum x^2
    coef[1][2] = xsums[2];

    coef[2][0] = xsums[1]; // sum x^2
    coef[2][1] = xsums[2]; // sum x^3
    coef[2][2] = xsums[3];

    equation.b[0] = ysums[0]; // sum y
    equation.b[1] = ysums[1]; // sum y*x
    equation.b[2] = ysums[2]; // sum y*x^2

    const double* a = solveEquation(equation);

    std::vector<double> result;
    result.assign(a, a + 3);

    delete[] a;

    double square_error = 0;

    for (int i = 0; i < N; ++i) { // a_1, a_2, a_3
        double Fsum = 0, x_coef = 1;

        for (int j = 0; j < 3; ++j) {
            Fsum += result[j] * x_coef;
            x_coef *= x[i];
        }

        double diff = Fsum - y[i];

        square_error += diff * diff;
    }

    std::string out;

    for (int i = 0; i < 3; ++i) {

        if (result[i] < 0) {
            out += '-';
            if (i != 0)
                out += ' ';
        } else if (i != 0) {
            out += "+ ";
        }

        out += std::format("{:.4f}", std::abs(result[i]));

        if (i > 0) {
            out += "x";
            if (i > 1) {
                out += "^" + std::to_string(i);
            }
        }

        out += ' ';
    }

    std::cout << out << std::endl;
    std::cout << std::format("Square Error: {:.4f}", square_error) << std::endl;

    return result;


}

std::vector<double> approximation_polynomial_third_degree(const std::vector<double>& x, const std::vector<double>& y) {
    const int N = x.size();

    LinearEquation<double> equation(4);
    double** coef = equation.matrix.coefficients;

    std::vector<double> xsums(6, 0.0), ysums(4, 0.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 6; ++j) {
            xsums[j] += std::pow(x[i], j + 1);
        }
        for (int j = 0; j < 4; ++j) {
            ysums[j] += y[i] * std::pow(x[i], j);
        }
    }

    // Матрица нормальных уравнений 4x4
    coef[0][0] = N;         coef[0][1] = xsums[0]; coef[0][2] = xsums[1]; coef[0][3] = xsums[2];
    coef[1][0] = xsums[0];  coef[1][1] = xsums[1]; coef[1][2] = xsums[2]; coef[1][3] = xsums[3];
    coef[2][0] = xsums[1];  coef[2][1] = xsums[2]; coef[2][2] = xsums[3]; coef[2][3] = xsums[4];
    coef[3][0] = xsums[2];  coef[3][1] = xsums[3]; coef[3][2] = xsums[4]; coef[3][3] = xsums[5];

    equation.b[0] = ysums[0]; // sum y
    equation.b[1] = ysums[1]; // sum y*x
    equation.b[2] = ysums[2]; // sum y*x^2
    equation.b[3] = ysums[3];

    const double* a = solveEquation(equation);
    std::vector<double> result(a, a + 4);
    delete[] a;

    double square_error = 0;
    for (int i = 0; i < N; ++i) { // копим Ф
        double Fsum = 0, x_coef = 1;
        for (int j = 0; j < 4; ++j) {
            Fsum += result[j] * x_coef;
            x_coef *= x[i];
        }
        double diff = Fsum - y[i];
        square_error += diff * diff;
    }
    std::string out;

    for (int i = 0; i < 4; ++i) {

        if (result[i] < 0) {
            out += '-';
            if (i != 0)
                out += ' ';
        } else if (i != 0) {
            out += "+ ";
        }

        out += std::format("{:.4f}", std::abs(result[i]));

        if (i > 0) {
            out += "x";
            if (i > 1) {
                out += "^" + std::to_string(i);
            }
        }

        out += ' ';
    }

    std::cout << out << std::endl;

    std::cout << std::format("Square Error: {:.4f}", square_error) << std::endl;
    return result;
}