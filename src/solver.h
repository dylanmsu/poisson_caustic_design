#ifndef SOLVER_H
#define SOLVER_H

#include <thread>
#include <cmath>
#include <string>

void poisson_solver(std::vector<std::vector<double>> &input, std::vector<std::vector<double>> &output, double hx, double hy, int max_iterations, double convergence_threshold, int max_threads);

#endif // SOLVER_H