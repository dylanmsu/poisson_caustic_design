#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "src/solver.h"
#include "src/mesh.h"
#include "src/polygon_utils.h"

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

void export_cells_as_svg(std::vector<polygon_t> cells, std::vector<double> intensities, std::string filename) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)1 / (double)1) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    for (int i=0; i<cells.size(); i++) {
        polygon_t cell = cells[i];
        std::string path_str = "M";
        for (std::size_t j = 0; j < cell.size(); ++j) {
            const auto& point = cell[j];
            path_str += std::to_string((point[0] / (double)1) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)1) * 1000.0f * ((double)1 / (double)1));

            if (j < cell.size() - 1)
                path_str += "L";
        }
        path_str += "Z";
        svg_file << "<path d=\"" << path_str << "\" fill=\"" << "rgb(" << intensities[i]*255 << ", " << intensities[i]*255 << ", " << intensities[i]*255 << ")\" stroke=\"red\" stroke-width=\"" << 1.0 << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();

}

std::vector<double> scale_array_proportional(const std::vector<double>& arr, double min_value, double max_value) {
    // Find the min and max values in the array
    auto arr_minmax = std::minmax_element(arr.begin(), arr.end());
    double arr_min = *arr_minmax.first;
    double arr_max = *arr_minmax.second;

    // Perform proportional scaling
    std::vector<double> scaled_array(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        scaled_array[i] = min_value + (max_value - min_value) * (arr[i] - arr_min) / (arr_max - arr_min);
    }

    return scaled_array;
}

std::vector<std::vector<double>> scale_matrix_proportional(const std::vector<std::vector<double>>& matrix, double min_value, double max_value) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    // Find the min and max values in the matrix
    double matrix_min = matrix[0][0];
    double matrix_max = matrix[0][0];

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix_min = std::min(matrix_min, matrix[i][j]);
            matrix_max = std::max(matrix_max, matrix[i][j]);
        }
    }

    // Perform proportional scaling
    std::vector<std::vector<double>> scaled_matrix(rows, std::vector<double>(cols));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            scaled_matrix[i][j] = min_value + (max_value - min_value) * (matrix[i][j] - matrix_min) / (matrix_max - matrix_min);
        }
    }

    return scaled_matrix;
}

std::vector<std::vector<std::vector<double>>> calculate_gradient(const std::vector<std::vector<double>>& grid) {
    int width = static_cast<int>(grid[0].size());
    int height = static_cast<int>(grid.size());

    std::vector<std::vector<std::vector<double>>> gradient3D(2, std::vector<std::vector<double>>(height, std::vector<double>(width, 0.0)));

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Calculate x gradient
            gradient3D[0][y][x] = (grid[y][std::min(x + 1, width - 1)] - grid[y][std::max(x - 1, 0)]) / 2.0;

            // Calculate y gradient
            gradient3D[1][y][x] = (grid[std::min(y + 1, height - 1)][x] - grid[std::max(y - 1, 0)][x]) / 2.0;
        }
    }

    return gradient3D;
}

void subtractAverage(std::vector<std::vector<double>>& raster) {
    // Calculate the average of the raster
    double sum = 0.0;
    int count = 0;

    for (const auto& row : raster) {
        for (double value : row) {
            sum += value;
            count++;
        }
    }

    double average = sum / count;

    // Subtract the average from each element of the raster
    for (auto& row : raster) {
        std::transform(row.begin(), row.end(), row.begin(), [average](double value) {
            return value - average;
        });
    }
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

std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny) {
    std::vector<std::vector<double>> divergence(ny, std::vector<double>(nx, 0.0));

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            // Forward and backward differences for divergence with boundary checks
            double div_x = (x < nx - 1) ? (Nx[y][x + 1] - Nx[y][x]) : 0.0;
            double div_y = (y < ny - 1) ? (Ny[y + 1][x] - Ny[y][x]) : 0.0;

            // Accumulate the divergence
            if (x == 0 || x == nx - 1 || y == 0 || y == ny - 1) {
                divergence[y][x] = 0.0;
            } else {
                divergence[y][x] = div_x + div_y;
            }
        }
    }

    return divergence;
}

int main(int argc, char const *argv[])
{
    int resolution_x = 750;
    int resolution_y = 750;

    std::vector<std::vector<double>> phi;
    std::vector<double> errors;

    std::string image_path = cv::samples::findFile("../img/nerdland.png");
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

        errors.clear();
        for (int i=0; i<target_areas.size(); i++) {
            errors.push_back(target_areas[i] - source_areas[i]);
        }

        double error_sum = 0;
        for (int i=0; i<target_areas.size(); i++) {
            errors[i] = errors[i] / calculate_polygon_area_vec(cells[i]);
            error_sum += errors[i];
        }
        double average = error_sum / target_areas.size();

        for (int i=0; i<target_areas.size(); i++) {
            errors[i] -= average;
        }

        //export_cells_as_svg(cells, scale_array_proportional(errors, 0.0f, 1.0f), "../cells.svg");
        
        printf("hello world!\r\n");

        mesh.build_bvh(1, 30);

        std::vector<std::vector<double>> raster = mesh.interpolate_raster(errors, resolution_x, resolution_y);

        printf("interpolated\r\n");

        //mesh.export_to_svg("../output.svg", 1.0f);
        
        std::string png_filename = "../interpolated_" + std::to_string(itr) + ".png";
        //std::string png_filename = "../interpolated.png";
        save_grid_as_image(raster, resolution_x, resolution_y, png_filename);

        //show_grid(scale_matrix_proportional(raster, 0.0f, 1.0f), resolution_x, resolution_y);

        phi.clear();
        for (int i = 0; i < resolution_x; ++i) {
            std::vector<double> row;
            for (int j = 0; j < resolution_y; ++j) {
                row.push_back(0.0f);
            }
            phi.push_back(row);
        }

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

        //mesh.export_to_svg("../output.svg", 1.0f);
        std::string svg_filename = "../parameterization.svg";
        mesh.export_paramererization_to_svg(svg_filename, 1.0f);

        mesh.build_bvh(1, 30);
        mesh.calculate_inverted_transport_map("../inverted.svg", 1.0f);

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
