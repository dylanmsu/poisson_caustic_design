#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "src/caustic_design.h"

#include <thread>

#define cimg_use_png
//#define cimg_use_jpeg

#include "cimg/CImg.h"
#include "lsqcpp/lsqcpp.hpp"

Caustic_design caustic_design;

void image_to_grid(const cimg_library::CImg<unsigned char>& image, std::vector<std::vector<double>>& image_grid) {
    for (int i = 0; i < image.height(); ++i) {
        std::vector<double> row;
        for (int j = 0; j < image.width(); ++j) {
            double r = image(j, i, 0) / 255.0; // Normalize R component
            double g = image(j, i, 1) / 255.0; // Normalize G component
            double b = image(j, i, 2) / 255.0; // Normalize B component
            double gray = (0.299 * r) + (0.587 * g) + (0.114 * b); // Calculate grayscale value using luminosity method
            row.push_back(gray);
        }
        image_grid.push_back(row);
    }
}

void save_grid_as_image(const std::vector<std::vector<double>>& img, int resolution_x, int resolution_y, const std::string& filename) {
    // Create an empty CImg object with the specified resolution
    cimg_library::CImg<unsigned char> image(resolution_x, resolution_y);

    // Copy the grid data to the image
    for (int i = 0; i < resolution_y; ++i) {
        for (int j = 0; j < resolution_x; ++j) {
            image(j, i) = static_cast<unsigned char>(img[i][j] * 255); // Scale to [0, 255] and cast to unsigned char
        }
    }

    // Save the image as a PNG
    image.save(filename.c_str());
}

std::unordered_map<std::string, std::string> parse_arguments(int argc, char const *argv[]) {
    // Define a map to store the parsed arguments
    std::unordered_map<std::string, std::string> args;

    // Iterate through command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string key, value;

        // Check if argument starts with '--'
        if (arg.substr(0, 2) == "--") {
            // Split argument by '=' to separate key and value
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                key = arg.substr(2, pos - 2);
                value = arg.substr(pos + 1);
            }
            else {
                key = arg.substr(2);
                value = ""; // No value provided
            }
        }
        // Check if argument starts with '-'
        else if (arg[0] == '-') {
            // The next argument is the value
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                key = arg.substr(1);
                value = argv[++i];
            }
            else {
                key = arg.substr(1);
                value = ""; // No value provided
            }
        }
        // Invalid argument format
        else {
            //std::cerr << "Invalid argument format: " << arg << std::endl;
            //return 1;
        }

        // Store key-value pair in the map
        args[key] = value;
    }

    return args;
}

struct ParabolicError
{
    static constexpr bool ComputesJacobian = false; // Set to true to indicate that the functor computes the Jacobian.

    template<typename Scalar, int Inputs, int Outputs>
    void operator()(const Eigen::Matrix<Scalar, Inputs, 1> &xval, Eigen::Matrix<Scalar, Outputs, 1> &fval) const {
        fval.resize(xval.size());

        // set the heights of the mesh based on the input vector
        for(lsqcpp::Index i = 0; i < xval.size(); ++i) {
            caustic_design.mesh->source_points[i][2] = xval(i);
        }

        Scalar E_int = 0.0f;
        for (lsqcpp::Index i = 0; i < xval.size(); ++i) {
            // Calculate the current normal
            std::vector<double> normal = caustic_design.calculate_vertex_normal(caustic_design.mesh->source_points, i);

            // Get the target normal
            std::vector<double> normal_trg = { caustic_design.normals[0][i], caustic_design.normals[1][i], caustic_design.normals[2][i]};

            normal_trg = normalize(normal_trg);
            normal = normalize(normal);

            // Calculate the difference
            std::vector<double> diff = vector_subtract(normal, normal_trg);

            fval(i) = diff[0];

            // Sum up the squared components of the difference
            //double energy = diff[0] * diff[0];// + diff[1] * diff[1] + diff[2] * diff[2];

            //E_int += energy;
        }

        //fval(0) = E_int;
    }
};







int main(int argc, char const *argv[])
{
    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    // Load image
    cimg_library::CImg<unsigned char> image(args["intput_png"].c_str());
    double aspect_ratio = (double)image.width() / (double)image.height();

    // Print parsed arguments
    /*for (const auto& pair : args) {
        //std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
        printf("Key: %s, Value: %s\r\n", pair.first.c_str(), pair.second.c_str());
    }*/

    int mesh_resolution_x = atoi(args["res_w"].c_str());
    double mesh_width = std::stod(args["width"]);

    caustic_design.set_mesh_resolution(mesh_resolution_x, mesh_resolution_x / aspect_ratio);
    caustic_design.set_domain_resolution(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio);

    double mesh_height = floor((mesh_resolution_x) / aspect_ratio) * (mesh_width / (mesh_resolution_x));

    caustic_design.set_mesh_size(mesh_width, mesh_height);

    caustic_design.set_lens_focal_length(std::stod(args["focal_l"]));
    caustic_design.set_lens_thickness(std::stod(args["thickness"]));
    caustic_design.set_solver_max_threads(atoi(args["max_threads"].c_str()));

    image = image.resize(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, -100, -100, 3); // Resize using linear interpolation

    // Convert image to grid
    std::vector<std::vector<double>> pixels;
    
    image_to_grid(image, pixels);

    caustic_design.initialize_solvers(pixels);

    printf("a\r\n");

    //caustic_design.export_paramererization_to_svg("../parameterization_0.svg", 0.5f);

    printf("b\r\n");
    
    for (int itr=0; itr<30; itr++) {
        printf("starting iteration %i\r\n", itr);
        
        double step_size = caustic_design.perform_transport_iteration();

        //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(caustic_design.vertex_gradient[0], 0.0f, 1.0f), "../x_grad.svg");

        //save_grid_as_image(scale_matrix_proportional(caustic_design.gradient[0], 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../grad_x_" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.gradient[1], 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../grad_y_" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.raster, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../raster_" + std::to_string(itr) + ".png");

        caustic_design.export_paramererization_to_svg("../parameterization_" + std::to_string(itr + 1) + ".svg", 1.0f);
        caustic_design.export_inverted_transport_map("../inverted.svg", 1.0f);

        printf("step_size = %f\r\n", step_size);

        // convergence treshold for the parameterization
        if (step_size < 0.005) break;
    }

    printf("\033[0;32mTransport map solver done! Starting height solver.\033[0m\r\n");

    for (int itr=0; itr<3; itr++) {
        //caustic_design.perform_height_map_iteration(itr);
        //save_grid_as_image(scale_matrix_proportional(caustic_design.h, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../h" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.divergence, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../div" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.norm_x, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "norm_x" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.norm_y, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "norm_y" + std::to_string(itr) + ".png");
    }

    caustic_design.normals.clear();
    caustic_design.normals = caustic_design.mesh->calculate_refractive_normals_uniform(caustic_design.focal_l * 16, 1.49);

    // Create GaussNewton optimizer object with ParabolicError functor as objective.
    // There are GradientDescent, GaussNewton and LevenbergMarquardt available.
    //
    // You can specify a StepSize functor as template parameter.
    // There are ConstantStepSize, BarzilaiBorwein, ArmijoBacktracking
    // WolfeBacktracking available. (Default for GaussNewton is ArmijoBacktracking)
    //
    // You can additionally specify a Callback functor as template parameter.
    //
    // You can additionally specify a FiniteDifferences functor as template
    // parameter. There are Forward-, Backward- and CentralDifferences
    // available. (Default is CentralDifferences)
    //
    // For GaussNewton and LevenbergMarquardt you can additionally specify a
    // linear equation system solver.
    // There are DenseSVDSolver and DenseCholeskySolver available.
    lsqcpp::GaussNewtonX<double, ParabolicError, lsqcpp::ArmijoBacktracking> optimizer;
    //lsqcpp::LevenbergMarquardtX<double, ParabolicError> optimizer;

    // Set number of iterations as stop criterion.
    // Set it to 0 or negative for infinite iterations (default is 0).
    optimizer.setMaximumIterations(100);

    // Set the minimum length of the gradient.
    // The optimizer stops minimizing if the gradient length falls below this
    // value.
    // Set it to 0 or negative to disable this stop criterion (default is 1e-9).
    optimizer.setMinimumGradientLength(1e-6);

    // Set the minimum length of the step.
    // The optimizer stops minimizing if the step length falls below this
    // value.
    // Set it to 0 or negative to disable this stop criterion (default is 1e-9).
    optimizer.setMinimumStepLength(1e-6);

    // Set the minimum least squares error.
    // The optimizer stops minimizing if the error falls below this
    // value.
    // Set it to 0 or negative to disable this stop criterion (default is 0).
    optimizer.setMinimumError(0);

    // Set the parameters of the step refiner (Armijo Backtracking).
    //optimizer.setMethodParameters({1.0, 2.0, 0.5, 100});

    // Turn verbosity on, so the optimizer prints status updates after each
    // iteration.
    optimizer.setVerbosity(2);

    std::vector<double> guess(caustic_design.mesh->source_points.size(), 0.0f);

    // Set initial guess.
    Eigen::VectorXd initialGuess = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(guess.data(), guess.size());

    // Start the optimization.
    auto result = optimizer.minimize(initialGuess);

    //printf("Done! Converged: %s Iterations: %d\n", result.converged ? "true" : "false", result.iterations);

    // do something with final function value
    //printf("Final fval: %s\n", result.fval.transpose());

    /*bool miss = false;
    std::vector<std::vector<double>> norm_x = caustic_design.mesh->interpolate_raster_source(x_normals, 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, miss);
    std::vector<std::vector<double>> norm_y = caustic_design.mesh->interpolate_raster_source(y_normals, 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, miss);
    save_grid_as_image(scale_matrix_proportional(norm_x, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../x_normals.png");
    save_grid_as_image(scale_matrix_proportional(norm_y, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../y_normals.png");
    //*/

    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(E_int, 0.0f, 1.0f), "../integration_energy.svg");
    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(y_normals_trg, 0.0f, 1.0f), "../y_normals_trg.svg");
    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(x_normals, 0.0f, 1.0f), "../x_normals.svg");
    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(y_normals, 0.0f, 1.0f), "../y_normals.svg");

    printf("Height solver done! Exporting as solidified obj\r\n");

    //caustic_design.save_solid_obj_target("../output.obj");
    caustic_design.save_solid_obj_source("../output.obj");

    return 0;
}
