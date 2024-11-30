#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <mutex>

#include "solver.h"

double solver_progress = 0.0f;

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

void calculate_progress(double value, double minValue, double maxValue) {
    // Calculate the percentage completion
    double solver_progress = static_cast<double>(value - minValue) / (maxValue - minValue);

    // Define the width of the progress bar
    int barWidth = 50;

    // Calculate how much of the progress bar should be filled
    int position = barWidth * solver_progress;

    // Render the progress bar
    printf("[");
    for (int i = 0; i < barWidth; ++i) {
        if (i < position) printf("=");      // Completed part of the bar
        else if (i == position) printf(">");// Current progress marker
        else printf(" ");                   // Remaining part of the bar
    }

    // Print the percentage at the end
    printf("] %.2f%%\r", solver_progress * 100.0);
    fflush(stdout);  // Ensure the output is displayed immediately
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

        // Create tasks to process the grid in equal chunks
        for (int y = 0; y < num_segments_y; y++) {
            for (int x = 0; x < num_segments_x; x++) {
                int start_x = floor(((double)(x)) * ((double)(width)/(double)(num_segments_x)));
                int start_y = floor(((double)(y)) * ((double)(height)/(double)(num_segments_y)));
                int end_x = ceil(((double)(x + 1)) * ((double)(width)/(double)(num_segments_x)));
                int end_y = ceil(((double)(y + 1)) * ((double)(height)/(double)(num_segments_y)));

                // Ensure end coordinates don't go beyond grid boundaries
                end_x = std::min(end_x, width);
                end_y = std::min(end_y, height);
                
                //threads[num_segments_y * x + y] = std::thread(process_grid_part,0,0,width, height);
                threads[num_segments_y * x + y] = std::thread(process_grid_part, start_x, start_y, end_x, end_y);
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
        if (i % 10 == 0) {
            //printf("\33[2K\r");
            calculate_progress(log(1.0f / max_update), log(1.0f / initial_max_update), log(1.0f / convergence_threshold));
            //printf("Solver progress: %f%\r\n", solver_progress * 100.0f);
        }

        // check for convergance
        if (max_update < convergence_threshold) {
            break;
            printf("\r\n");
        }
    }
}