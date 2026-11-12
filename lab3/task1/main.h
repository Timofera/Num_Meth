#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <cmath>
#include <iostream>
#include <vector>
#include <format>
#include <string>

using func = double(*)(double);

double function(double x);

void print_nodes(func y, const std::vector<double>& x_points);

double Lagrange_polynomial(func y, const std::vector<double>& x_points, double x_precision);

double Newton_polynomial(func y, const std::vector<double>& x_points, double x_precision);

#endif //NUMERIC_METHODS_MAIN_H