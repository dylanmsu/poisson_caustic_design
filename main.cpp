#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#define OPENCV_DISABLE_THREAD_SUPPORT ON

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "src/solver.h"
#include "src/mesh.h"
#include "src/polygon_utils.h"
#include "src/utils.h"

void image_to_grid(cv::Mat image, std::vector<std::vector<double>>& image_grid) {
    for (int i = 0; i < image.rows; ++i) {
        std::vector<double> row;
        for (int j = 0; j < image.cols; ++j) {
            cv::Vec3b intensity = image.at<cv::Vec3b>(i, j); // Retrieve BGR values of the pixel
            double b = intensity[0] / 255.0; // Normalize B component
            double g = intensity[1] / 255.0; // Normalize G component
            double r = intensity[2] / 255.0; // Normalize R component
            double gray = (0.299 * r) + (0.587 * g) + (0.114 * b); // Calculate grayscale value using luminosity method
            row.push_back(gray);
        }
        image_grid.push_back(row);
    }
}

void show_grid(std::vector<std::vector<double>> image_grid, int width, int height) {
    cv::Mat grayImg(height, width, CV_8U);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            grayImg.at<uchar>(y, x) = static_cast<uchar>(image_grid[y][x] * 255);
        }
    }

    cv::Mat grayImg8bit;
    cv::convertScaleAbs(grayImg, grayImg8bit);
    cv::imshow("Grayscale Image", grayImg8bit);
    int k = cv::waitKey(0);
}

void save_grid_as_image(std::vector<std::vector<double>>& img, int resolution_x, int resolution_y, std::string filename) {
    // Convert the 'raster' data to a suitable format for OpenCV
    cv::Mat image(resolution_y, resolution_x, CV_64FC1);

    for (int i = 0; i < resolution_y; ++i) {
        for (int j = 0; j < resolution_x; ++j) {
            image.at<double>(i, j) = img[i][j];
        }
    }

    // Normalize the values to the range [0, 255]
    cv::normalize(image, image, 0, 255, cv::NORM_MINMAX);

    // Convert to 8-bit unsigned integer format
    image.convertTo(image, CV_8UC1);

    // Save the image as a PNG
    cv::imwrite(filename, image);
}

int main(int argc, char const *argv[])
{
    int resolution_x = 750;
    int resolution_y = 750;

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

    std::string image_path = cv::samples::findFile("../img/einstein.png");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }

    cv::Mat scaledImg;
    cv::resize(img, scaledImg, cv::Size(resolution_x, resolution_y), cv::INTER_LINEAR);

    //cv::bitwise_not(scaledImg, scaledImg);

    // convert image to grayscale values
    std::vector<std::vector<double>> pixels;
    image_to_grid(scaledImg, pixels);

    pixels = scale_matrix_proportional(pixels, 0, 1.0f);

    Mesh mesh(1.0f, 1.0f, 256, 256);

    std::cout << "built mesh" << std::endl;

    std::vector<std::vector<std::vector<double>>> cells;
    mesh.build_target_dual_cells(cells);

    std::cout << "generated dual cells" << std::endl;

    std::vector<double> target_areas = get_target_areas(pixels, cells, resolution_x, resolution_y, 1.0f, 1.0f);
    
    for (int itr=0; itr<200; itr++) {
        cells.clear();
        mesh.build_target_dual_cells(cells);

        std::vector<double> source_areas = get_source_areas(cells);

        calculate_errors(source_areas, target_areas, cells, errors);

        export_cells_as_svg(cells, scale_array_proportional(errors, 0.0f, 1.0f), "../cells.svg");

        mesh.build_bvh(1, 30);

        std::vector<std::vector<double>> raster = mesh.interpolate_raster(errors, resolution_x, resolution_y);

        printf("interpolated\r\n");

        //mesh.export_to_svg("../output.svg", 1.0f);
    
        //show_grid(scale_matrix_proportional(raster, 0.0f, 1.0f), resolution_x, resolution_y);

        subtractAverage(raster);
        poisson_solver(raster, phi, resolution_x, resolution_y, 100000, 0.0000001);

        //save_grid_as_image(phi, resolution_x, resolution_y, "../phi.png");

        //show_grid(scale_matrix_proportional(phi, 0.0f, 1.0f), resolution_x, resolution_y);

        std::vector<std::vector<std::vector<double>>> grad = calculate_gradient(phi);

        //save_grid_as_image(grad[0], resolution_x, resolution_y, "../grad_x.png");

        //show_grid(scale_matrix_proportional(grad[0], 0.0f, 1.0f), resolution_x, resolution_y);
        //show_grid(scale_matrix_proportional(grad[1], 0.0f, 1.0f), resolution_x, resolution_y);
        
        std::vector<std::vector<double>> gradient = integrate_cell_gradients(grad, cells, resolution_x, resolution_y, 1.0f, 1.0f);

        //export_cells_as_svg(cells, scale_array_proportional(gradient[0], 0.0f, 1.0f), "../grad_x.svg");
        //export_cells_as_svg(cells, scale_array_proportional(gradient[1], 0.0f, 1.0f), "../grad_y.svg");

        printf("integrated\r\n");

        double min_step = mesh.step_grid(gradient[0], gradient[1], 0.1);

        std::string svg_filename = "../parameterization.svg";
        mesh.export_paramererization_to_svg(svg_filename, 1.0f);

        mesh.build_bvh(1, 30);
        mesh.calculate_inverted_transport_map("../inverted.svg", 1.0f);

        std::string png_filename = "../interpolated_" + std::to_string(itr) + "_" + std::to_string(min_step) + ".png";
        //std::string png_filename = "../interpolated.png";
        save_grid_as_image(raster, resolution_x, resolution_y, png_filename);


        if (min_step < 1e-6) break;
    }

    for (int itr=0; itr<3; itr++) {
        std::vector<std::vector<double>> normals = mesh.calculate_refractive_normals(resolution_x / 1.0f * 2.0f, 1.49);

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

        std::vector<double> interpolated_h = integrate_grid_into_cells(h, source_cells, resolution_x, resolution_y, 1.0f, 1.0f);

        mesh.set_source_heights(interpolated_h);
    }

    mesh.save_solid_obj(0.2f, "../output.obj");

    return 0;
}
