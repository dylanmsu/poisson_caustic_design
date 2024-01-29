#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "solver.h"
#include "mesh.h"

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

int main(int argc, char const *argv[])
{
    /*int resolution_x = 500;
    int resolution_y = 500;

    std::string image_path = cv::samples::findFile("../img/lena.png");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }

    cv::Mat scaledImg;
    cv::resize(img, scaledImg, cv::Size(resolution_x, resolution_y), cv::INTER_LINEAR);

    // convert image to grayscale values
    std::vector<std::vector<double>> pixels;
    image_to_grid(scaledImg, pixels);

    show_grid(pixels, resolution_x, resolution_y);*/

    //std::vector<std::vector<double>> image_grid;
    //for (int i = 0; i < resolution_x; ++i) {
    //    std::vector<double> row;
    //    for (int j = 0; j < resolution_y; ++j) {
    //        row.push_back(0.0f);
    //    }
    //    image_grid.push_back(row);
    //}

    //poisson_solver(pixels, image_grid, resolution_x, resolution_y, 1000, 0.0000001);

    //show_grid(image_grid, resolution_x, resolution_y);

    Mesh mesh(1.0f, 1.0f, 500, 500);

    std::cout << "built mesh" << std::endl;

    mesh.export_to_svg("../output.svg", 1.0f);

    std::vector<polygon_t> cells = mesh.build_dual_cells();

    std::ofstream svg_file("../cell.svg", std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << "../cell.svg" << std::endl;
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
                        std::to_string((1 - point[1] / (double)1) * 1000.0f * ((double)1 / (double)1));

            if (j < cell.size() - 1)
                path_str += "L";
        }
        path_str += "Z";
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << 1.0 << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();

    printf("hello world!\r\n");
    return 0;
}
