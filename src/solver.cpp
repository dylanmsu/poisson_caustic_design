#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <mutex>

#include "solver.h"

double solver_progress = 0.0f;

double get_progress() {
    return solver_progress;
}

// perform a relaxation step
double patial_relax(std::vector<std::vector<double>> &output, std::vector<std::vector<double>> &input, int width, int height, double omega, int start_x, int start_y, int end_x, int end_y) {
    int x, y;
    double max_update = 0.0;

    // Precompute width and height - 1
    int width_minus_1 = width - 1;
    int height_minus_1 = height - 1;

    for (y = start_y; y < end_y; y++) {
        for (x = start_x; x < end_x; x++) {
            double val = output[y][x];
            double delta;

            double neighbor_sum = 0.0;
            double neighbor_cnt = 0.0;

            // Neumann boundary conditions for domain edges and NAN values
            if (x != 0 && !std::isnan(input[y][x - 1])) {
                neighbor_cnt += 1.0;
                neighbor_sum += output[y][x - 1];
            }
            if (y != 0 && !std::isnan(input[y - 1][x])) {
                neighbor_cnt += 1.0;
                neighbor_sum += output[y - 1][x];
            }
            if (x != width_minus_1 && !std::isnan(input[y][x + 1])) {
                neighbor_cnt += 1.0;
                neighbor_sum += output[y][x + 1];
            }
            if (y != height_minus_1 && !std::isnan(input[y + 1][x])) {
                neighbor_cnt += 1.0;
                neighbor_sum += output[y + 1][x];
            }

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

void calculate_progress(int value, int minValue, int maxValue) {
    // Calculate the percentage completion
    solver_progress = static_cast<double>((value - minValue) / (maxValue - minValue));
}

void poisson_solver(std::vector<std::vector<double>> &input, std::vector<std::vector<double>> &output, int width, int height, int max_iterations, double convergence_threshold, int max_threads) {
    double omega = 2.0 / (1.0 + 3.14159265 / width);

    int num_threads = std::thread::hardware_concurrency();

    double initial_max_update = 0.0f;

    num_threads = fmin(max_threads, num_threads);

    // process grid into all equal subsets
    int num_segments_x = floor(sqrt(num_threads));
    int num_segments_y = num_segments_x;

    std::vector<std::thread> threads(num_segments_x * num_segments_y);

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

        // Create tasks to process the grid in equal chunks
        for (int y = 0; y < num_segments_y; y++) {
            for (int x = 0; x < num_segments_x; x++) {
                //threads[num_segments_y * x + y] = std::thread(process_grid_part,0,0,width, height);
                threads[num_segments_y * x + y] = std::thread(process_grid_part, 
                    floor(((double)(x)) * ((double)(width/num_segments_x))),
                    floor(((double)(y)) * ((double)(height/num_segments_y))),
                    ceil(((double)(x + 1)) * ((double)(width/num_segments_x))),
                    ceil(((double)(y + 1)) * ((double)(height/num_segments_y))));
            }
        }

        // Join threads
        for (int i = 0; i < num_segments_x * num_segments_y; ++i) {
            threads[i].join();
        }//*/

        if (i == 0) {
            initial_max_update = max_update;
        }

        // print progress to terminal
        if (i % 100 == 0) {
            printf("\33[2K\r");
            printf("Solver max_update: %f, convergence at %.2e\r", log(1.0f / max_update), convergence_threshold);

            calculate_progress(log(1.0f / max_update), log(1.0f / initial_max_update), log(1.0f / convergence_threshold));
        }

        // check for convergance
        if (max_update < convergence_threshold) {
            printf("\r\n");
            solver_progress = -1.0f;
            break;
        }
    }
}