#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "solver.h"

// perform a relaxation step
void relax(std::vector<std::vector<double>> &output, std::vector<std::vector<double>> &input, int width, int height, double omega, double &max_update) {
    int x, y;
    max_update = 0.0;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            double val = output[y][x];
            double delta;

            // corner 1
            if (x == 0 && y == 0) {
                double val_down = output[(y + 1)][x];
                double val_right = output[(y)][(x + 1)];
                delta = omega / 2 * (val_down + val_right - 2 * val - input[y][x]);
            }

            // corner 2
            else if (x == 0 && y == (height-1)) {
                double val_up = output[(y - 1)][x];
                double val_right = output[(y)][(x + 1)];
                delta = omega / 2 * (val_up + val_right - 2 * val - input[y][x]);
            }

            // corner 3
            else if (x == (width-1) && y == 0) {
                double val_down = output[(y + 1)][x];
                double val_left = output[(y)][(x - 1)];
                delta = omega / 2 * (val_down + val_left - 2 * val - input[y][x]);
            }

            // corner 4
            else if (x == (width-1) && y == (height-1)){
                double val_up = output[(y - 1)][x];
                double val_left = output[(y)][(x - 1)];
                delta = omega / 2 * (val_up + val_left - 2 * val - input[y][x]);
            }
            
            // side 1
            else if (x == 0) {
                double val_up = output[(y - 1)][x];
                double val_down = output[(y + 1)][x];
                double val_right = output[(y)][(x + 1)];
                delta = omega / 3 * (val_up + val_down + val_right - 3 * val - input[y][x]);
            }

            // side 2
            else if (x == (width-1)) {
                double val_up = output[(y - 1)][x];
                double val_down = output[(y + 1)][x];
                double val_left = output[(y)][(x - 1)];
                delta = omega / 3 * (val_up + val_down + val_left - 3 * val - input[y][x]);
            }

            // side 3
            else if (y == 0) {
                double val_down = output[(y + 1)][x];
                double val_left = output[(y)][(x - 1)];
                double val_right = output[(y)][(x + 1)];
                delta = omega / 3 * (val_down + val_left + val_right - 3 * val - input[y][x]);
            }

            // side 4
            else if (y == (height-1)) {
                double val_up = output[(y - 1)][x];
                double val_left = output[(y)][(x - 1)];
                double val_right = output[(y)][(x + 1)];
                delta = omega / 3 * (val_up + val_left + val_right - 3 * val - input[y][x]);
            }

            // all other non-boundary points
            else {
                double val_up = output[(y - 1)][x];
                double val_down = output[(y + 1)][x];
                double val_left = output[(y)][(x - 1)];
                double val_right = output[(y)][(x + 1)];
                delta = omega / 4 * (val_up + val_down + val_left + val_right - 4 * val - input[y][x]);
            }

            // get the minimum max_update as the convergance indicator
            double abs_delta = fabs(delta);
            if (abs_delta > max_update) {
                max_update = abs_delta;
            }

            // increment the value by delta
            output[y][x] += delta;
        }
    }
}

void poisson_solver(std::vector<std::vector<double>> &input, std::vector<std::vector<double>> &output, int width, int height, int max_iterations, double convergence_threshold) {
    double max_update = 0.0;
    double omega = 2.0 / (1.0 + 3.14159265 / width);
    
    // set the initial guess for the solution to all zero's
    /*for (int i=0; i<width; i++) {
        for (int j=0; j<width; j++) {
            output[j][i] = 0.0f;
        }
    }*/

    for (int i = 0; i < max_iterations; i++) {
        relax(output, input, width, height, omega, max_update);

        if (i % 100 == 0) {
            printf("\33[2K\r");
            printf("Max update: %f\n\r", max_update);
        }

        if (max_update < convergence_threshold) {
            printf("Convergence reached at iteration %d with max_update of %lf\n", i, max_update);
            break;
        }
    }
}