#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <mutex>

#include "solver.h"

// perform a relaxation step
double patial_relax(std::vector<std::vector<double>> &output, std::vector<std::vector<double>> &input, int width, int height, double omega, int start_x, int start_y, int end_x, int end_y) {
    int x, y;
    double max_update = 0.0;

    for (y = start_y; y < end_y; y++) {
        for (x = start_x; x < end_x; x++) {
            double val = output[y][x];
            double delta;

            double neighbor_sum = 0.0f;
            double neighbor_cnt = 0.0f;

            // boundary conditions
            if (x != 0 && !std::isnan(input[y][x - 1]))             {neighbor_cnt += 1.0f; neighbor_sum += output[y][x - 1];}
            if (y != 0 && !std::isnan(input[y - 1][x]))             {neighbor_cnt += 1.0f; neighbor_sum += output[y - 1][x];}
            if (x != (width-1) && !std::isnan(input[y][x + 1]))     {neighbor_cnt += 1.0f; neighbor_sum += output[y][x + 1];}
            if (y != (height-1) && !std::isnan(input[y + 1][x]))    {neighbor_cnt += 1.0f; neighbor_sum += output[y + 1][x];}

            // calculate delta
            delta = omega / neighbor_cnt * (neighbor_sum - neighbor_cnt * val - input[y][x]);

            // get max update
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

void poisson_solver(std::vector<std::vector<double>> &input, std::vector<std::vector<double>> &output, int width, int height, int max_iterations, double convergence_threshold, int max_threads) {
    double omega = 2.0 / (1.0 + 3.14159265 / width);

    int num_threads = std::thread::hardware_concurrency();

    num_threads = fmin(max_threads, num_threads);

    // process grid into all equal subsets
    int num_segments_x = floor(sqrt(num_threads));
    int num_segments_y = num_segments_x;

    std::vector<std::thread> threads(num_segments_x * num_segments_y);
    
    // set the initial guess for the solution to all zero's -> removed
    /*for (int i=0; i<width; i++) {
        for (int j=0; j<width; j++) {
            output[j][i] = 0.0f;
        }
    }*/

    for (int i = 0; i < max_iterations; i++) {
        double max_update = 0.0;

        std::mutex mtx;

        // Function to process a portion of the grid
        auto process_grid_part = [&](int start_x, int start_y, int end_x, int end_y) {
            double local_max_update = 0.0;
            local_max_update = patial_relax(output, input, width, height, omega, start_x, start_y, end_x, end_y);
            
            mtx.lock();
            if (max_update < local_max_update) {
                max_update = local_max_update;
            }
            mtx.unlock();
        };

        // Create tasks to process the grid in chunks
        for (int y = 0; y < num_segments_y; y++) {
            for (int x = 0; x < num_segments_x; x++) {
                threads[num_segments_y * x + y] = std::thread(process_grid_part, 
                    x * (width/num_segments_x),
                    y * (height/num_segments_y),
                    (x + 1) * (width/num_segments_x),
                    (y + 1) * (height/num_segments_y));
            }
        }

        // Join threads
        for (int i = 0; i < num_segments_x * num_segments_y; ++i) {
            threads[i].join();
        }//*/

        if (i % 100 == 0) {
            printf("\33[2K\r");
            printf("Solver max_update: %.5e, convergence at %.2e\r", max_update, convergence_threshold);
        }

        if (max_update < convergence_threshold) {
            printf("Convergence reached at iteration %d with max_update of %.5e\n", i, max_update);
            break;
        }
    }
}