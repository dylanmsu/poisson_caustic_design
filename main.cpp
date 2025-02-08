#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <png.h>

#include "src/caustic_design.h"

#include <thread>

void image_to_grid(const std::string& filename, std::vector<std::vector<double>>& image_grid) {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cout << "FILENAME : " << filename << std::endl;
        throw std::runtime_error("Failed to open PNG file.");
    }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fclose(fp);
        throw std::runtime_error("Failed to create PNG read struct.");
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_read_struct(&png, NULL, NULL);
        fclose(fp);
        throw std::runtime_error("Failed to create PNG info struct.");
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        throw std::runtime_error("Error during PNG read initialization.");
    }

    png_init_io(png, fp);
    png_read_info(png, info);

    int width = png_get_image_width(png, info);
    int height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth = png_get_bit_depth(png, info);

    if (bit_depth == 16) {
        png_set_strip_16(png);
    }

    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png);
    }

    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(png);
    }

    if (png_get_valid(png, info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(png);
    }

    png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

    png_read_update_info(png, info);

    std::vector<png_bytep> row_pointers(height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = (png_bytep)malloc(png_get_rowbytes(png, info));
    }

    png_read_image(png, row_pointers.data());

    fclose(fp);

    for (int i = 0; i < height; ++i) {
        std::vector<double> row;
        for (int j = 0; j < width; ++j) {
            png_bytep px = &row_pointers[i][j * 4];
            double r = px[0] / 255.0;
            double g = px[1] / 255.0;
            double b = px[2] / 255.0;
            double gray = (0.299 * r) + (0.587 * g) + (0.114 * b);
            row.push_back(gray);
        }
        image_grid.push_back(row);
    }

    for (int y = 0; y < height; y++) {
        free(row_pointers[y]);
    }

    png_destroy_read_struct(&png, &info, NULL);
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

        // Store key-value pair in the map
        args[key] = value;
    }

    return args;
}

int main(int argc, char const *argv[]) {
    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    // Load image to grid
    std::vector<std::vector<double>> pixels;
    image_to_grid(args["input_png"], pixels);
    double aspect_ratio = (double)pixels[0].size() / (double)pixels.size();

    Caustic_design caustic_design;

    int mesh_resolution_x = atoi(args["res_w"].c_str());
    double mesh_width = std::stod(args["width"]);

    caustic_design.set_mesh_resolution(mesh_resolution_x, mesh_resolution_x / aspect_ratio);
    caustic_design.set_domain_resolution(4 * mesh_resolution_x, 4 * mesh_resolution_x / aspect_ratio);

    double mesh_height = floor((mesh_resolution_x) / aspect_ratio) * (mesh_width / (mesh_resolution_x));

    caustic_design.set_mesh_size(mesh_width, mesh_height);

    caustic_design.set_lens_focal_length(std::stod(args["focal_l"]));
    caustic_design.set_lens_thickness(std::stod(args["thickness"]));
    caustic_design.set_solver_max_threads(atoi(args["max_threads"].c_str()));

    caustic_design.initialize_solvers(pixels);

    caustic_design.export_paramererization_to_svg("../parameterization_0.svg", 0.5f);

    for (int itr = 0; itr < 30; itr++) {
        printf("starting iteration %i\r\n", itr);

        double step_size = caustic_design.perform_transport_iteration();

        caustic_design.export_paramererization_to_svg("../parameterization_" + std::to_string(itr + 1) + ".svg", 1.0f);
        caustic_design.export_inverted_transport_map("../inverted.svg", 1.0f);

        printf("step_size = %f\r\n", step_size);

        if (step_size < 0.01) break;
    }

    printf("\033[0;32mTransport map solver done! Starting height solver.\033[0m\r\n");

    for (int itr = 0; itr < 3; itr++) {
        caustic_design.perform_height_map_iteration(itr);
    }

    printf("Height solver done! Exporting as solidified obj\r\n");

    caustic_design.save_solid_obj_source("../output.obj");

    return 0;
}
