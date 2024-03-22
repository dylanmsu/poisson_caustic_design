#ifndef SOLVER_H
#define SOLVER_H

#include <thread>
#include <cmath>
#include <string>

void poisson_solver(std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &phi, int width, int height, int max_iterations, double convergence_threshold, int max_threads);

double get_progress();

#endif // SOLVER_H