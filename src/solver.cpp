#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "solver.h"

// perform a relaxation step
double patial_relax(std::vector<std::vector<double>> &output, std::vector<std::vector<double>> &input, int width, int height, double omega, int start_x, int start_y, int end_x, int end_y) {
    int x, y;
    double max_update = 0.0;

    for (y = start_y; y < end_y; y++) {
        for (x = start_x; x < end_x; x++) {
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

    return max_update;
}

void poisson_solver(std::vector<std::vector<double>> &input, std::vector<std::vector<double>> &output, int width, int height, int max_iterations, double convergence_threshold) {
    double omega = 2.0 / (1.0 + 3.14159265 / width);

    int num_segments_x = 4;
    int num_segments_y = 4;

    std::vector<std::thread> threads(num_segments_x * num_segments_y);
    
    // set the initial guess for the solution to all zero's
    /*for (int i=0; i<width; i++) {
        for (int j=0; j<width; j++) {
            output[j][i] = 0.0f;
        }
    }*/

    for (int i = 0; i < max_iterations; i++) {
        double max_update = 0.0;

        //max_update = patial_relax(output, input, width, height, omega, 0, 0, width, height);

        // Function to process a portion of the grid
        auto process_grid_part = [&](int start_x, int start_y, int end_x, int end_y) {
            double local_max_update = 0.0;
            local_max_update = patial_relax(output, input, width, height, omega, start_x, start_y, end_x, end_y);
            if (max_update < local_max_update) {
                max_update = local_max_update;
            }
        };

        for (int y = 0; y < num_segments_x; y++) {
            for (int x = 0; x < num_segments_x; x++) {
                threads[x + num_segments_x * y] = std::thread(process_grid_part, 
                    x * (width/num_segments_x), 
                    y * (height/num_segments_y), 
                    (x + 1) * (width/num_segments_x), 
                    (y + 1) * (height/num_segments_y));
            }
        }

        //threads[0] = std::thread(process_grid_part, 0,       0,          width/2,    height/2);
        //threads[1] = std::thread(process_grid_part, width/2, 0,          width,      height/2);
        //threads[2] = std::thread(process_grid_part, 0,       height/2,   width/2,    height);
        //threads[3] = std::thread(process_grid_part, width/2, height/2,   width,      height);

        // Join threads
        for (int i = 0; i < num_segments_x * num_segments_y; ++i) {
            threads[i].join();
        }//*/

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