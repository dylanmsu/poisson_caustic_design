#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "src/caustic_design.h"
#include "src/normal_integration.h"

#include <thread>

#define cimg_use_png
//#define cimg_use_jpeg

#include "cimg/CImg.h"

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


int main(int argc, char const *argv[])
{
    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    // Load image
    std::string image_path = args["input_png"];
    cimg_library::CImg<unsigned char> image(image_path.c_str());
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
        if (step_size < 0.003) break;
    }

    printf("\033[0;32mTransport map solver done! Starting height solver.\033[0m\r\n");

    normal_integration normal_int;

    normal_int.initialize_data(*caustic_design.mesh);

    std::vector<std::vector<double>> desired_normals;

    for (int i=0; i<20; i++)
    {
        caustic_design.normals.clear();
        caustic_design.normals = caustic_design.mesh->calculate_refractive_normals_uniform(caustic_design.focal_l, 1.49);

        desired_normals.clear();

        // make a copy of the original positions of the vertices
        for (int i = 0; i < caustic_design.mesh->source_points.size(); i++) {
            std::vector<double> trg_normal = {caustic_design.normals[0][i], caustic_design.normals[1][i], caustic_design.normals[2][i]};
            desired_normals.push_back(trg_normal);
        }

        normal_int.perform_normal_integration(*caustic_design.mesh, desired_normals);
    }

    printf("\033[0;32mHeight solver done! Exporting as solidified obj\033[0m\r\n");

    //caustic_design.save_solid_obj_target("../output.obj");
    caustic_design.save_solid_obj_source("../output.obj");

    return 0;
}
