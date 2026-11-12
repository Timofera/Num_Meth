#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <iostream>
#include <cmath>
#include <vector>
#include <format>

struct SplineCoeffs {
    std::vector<double> a, b, c, d;
};

double find_h(std::vector<double> x, int index);
SplineCoeffs compute_spline_coeffs(const std::vector<double>& x,
                                   const std::vector<double>& f,
                                   double* solution);
double eval_spline(const SplineCoeffs& s, const std::vector<double>& x, int i, double px);
double cubic_spline(const std::vector<double>& x, const std::vector<double>& f, double point_x);

#endif //NUMERIC_METHODS_MAIN_H