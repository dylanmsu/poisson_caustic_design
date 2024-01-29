#ifndef SOLVER_H
#define SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

void relax(std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &phi, int width, int height, double omega, double &max_update);
void poisson_solver(std::vector<std::vector<double>> &D, std::vector<std::vector<double>> &phi, int width, int height, int max_iterations, double convergence_threshold);

#ifdef __cplusplus
}
#endif

#endif // SOLVER_H