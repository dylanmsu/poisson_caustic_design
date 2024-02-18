#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "src/solver.h"
#include "src/mesh.h"
#include "src/polygon_utils.h"
#include "src/utils.h"

#include <thread>

#define cimg_use_png
//#define cimg_use_jpeg

#include<X11/Xlib.h>
#include "/home/dylan/caustic_engineering/CImg.h"

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

int main(int argc, char const *argv[])
{
    printf("hello\r\n"); fflush(stdout);

    int mesh_res_x = 200;
    int mesh_res_y = 200;

    int resolution_x = 3*mesh_res_x;
    int resolution_y = 3*mesh_res_y;

    double width = 0.5;
    double height = 0.5;

    double focal_l = 1.5f;

    double thickness = 0.1;

    std::vector<std::vector<double>> phi;
    std::vector<double> errors;

    phi.clear();
    for (int i = 0; i < resolution_x; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_y; ++j) {
            row.push_back(0.0f);
        }
        phi.push_back(row);
    }

    printf("loading image\r\n");

    // Load image
    cimg_library::CImg<unsigned char> image("../img/lena.png");

    image = image.resize(resolution_x, resolution_y, -100, -100, 3); // Resize using linear interpolation

    // Convert image to grid
    std::vector<std::vector<double>> pixels;
    image_to_grid(image, pixels);

    printf("converted image to grid\r\n");

    //std::string image_path = cv::samples::findFile("C:/Users/dylan/Documents/caustic_engineering/img/face.png");
    //cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    /*if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }*/

    /*cv::Mat scaledImg;
    cv::resize(img, scaledImg, cv::Size(resolution_x, resolution_y), cv::INTER_LINEAR);*/

    //cv::bitwise_not(scaledImg, scaledImg);

    // convert image to grayscale values
    //std::vector<std::vector<double>> pixels;

    //printf("hello2\r\n"); fflush(stdout);
    //image_to_grid(scaledImg, pixels);

    pixels = scale_matrix_proportional(pixels, 0, 1.0f);

    Mesh mesh(width, height, mesh_res_x, mesh_res_y);

    printf("generated mesh\r\n");

    //std::cout << "built mesh" << std::endl;

    std::vector<std::vector<std::vector<double>>> target_cells;
    std::vector<std::vector<std::vector<double>>> source_cells;
    mesh.build_target_dual_cells(target_cells);
    mesh.build_source_dual_cells(source_cells);

    printf("generated dual cells\r\n");

    std::vector<double> target_areas = get_target_areas(pixels, target_cells, resolution_x, resolution_y, width, height);
    
    for (int itr=0; itr<500; itr++) {
        
        // build median dual mesh of the updated parameterization
        target_cells.clear();
        mesh.build_target_dual_cells(target_cells);

        // calculate 
        std::vector<double> source_areas = get_source_areas(target_cells);

        calculate_errors(source_areas, target_areas, target_cells, errors);

        export_cells_as_svg(target_cells, scale_array_proportional(errors, 0.0f, 1.0f), "../cells.svg");

        mesh.build_bvh(1, 30);

        std::vector<std::vector<double>> raster = mesh.interpolate_raster(errors, resolution_x, resolution_y);

        printf("interpolated\r\n");

        subtractAverage(raster);
        poisson_solver(raster, phi, resolution_x, resolution_y, 100000, 0.0000001);

        printf("generated dual cells\r\n");
        std::vector<std::vector<std::vector<double>>> grad = calculate_gradient(phi);

        std::vector<std::vector<double>> gradient = integrate_cell_gradients(grad, target_cells, resolution_x, resolution_y, width, height);

        printf("integrated\r\n");

        double min_step = mesh.step_grid(gradient[0], gradient[1], 0.05);
        printf("stepped grid\r\n");

        // iteration done, exporting visualizations:

        mesh.export_to_svg("../mesh.svg", 0.5);

        std::string svg_filename = "../parameterization.svg";
        mesh.export_paramererization_to_svg(svg_filename, 1.0f);

        mesh.build_bvh(1, 30);
        mesh.calculate_inverted_transport_map("../inverted.svg", 1.0f);

        //std::string png_filename = "../param/interpolated_" + std::to_string(itr) + "_" + std::to_string(min_step) + ".png";
        std::string png_filename = "../interpolated.png";
        raster = scale_matrix_proportional(raster, 0, 1.0f);
        save_grid_as_image(raster, resolution_x, resolution_y, png_filename);//*/

        printf("min_step = %f\r\n", min_step*(resolution_x/width));

        // convergence treshold for the parameterization
        if (min_step*(resolution_x/width) < 0.005) break;
    }

    for (int itr=0; itr<3; itr++) {
        std::vector<std::vector<double>> normals = mesh.calculate_refractive_normals(resolution_x / width * focal_l, 1.49);

        mesh.build_bvh(1, 30);
        std::vector<std::vector<double>> norm_x = mesh.interpolate_raster(normals[0], resolution_x, resolution_y);
        std::vector<std::vector<double>> norm_y = mesh.interpolate_raster(normals[1], resolution_x, resolution_y);

        std::vector<std::vector<double>> divergance = calculate_divergence(norm_x, norm_y, resolution_x, resolution_y);
    
        std::vector<std::vector<double>> h;
        for (int i = 0; i < resolution_x; ++i) {
            std::vector<double> row;
            for (int j = 0; j < resolution_y; ++j) {
                row.push_back(0.0f);
            }
            h.push_back(row);
        }

        subtractAverage(divergance);
        poisson_solver(divergance, h, resolution_x, resolution_y, 100000, 0.0000001);

        std::vector<std::vector<std::vector<double>>> source_cells;
        mesh.build_source_dual_cells(source_cells);

        std::vector<double> interpolated_h = integrate_grid_into_cells(h, source_cells, resolution_x, resolution_y, width, height);

        mesh.set_source_heights(interpolated_h);
    }

    mesh.save_solid_obj(thickness, "../output.obj");

    return 0;
}
