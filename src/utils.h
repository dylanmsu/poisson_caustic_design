#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "polygon_utils.h"

void subtractAverage(std::vector<std::vector<double>>& raster);
std::vector<std::vector<double>> scale_matrix_proportional(const std::vector<std::vector<double>>& matrix, double min_value, double max_value);
std::vector<double> scale_array_proportional(const std::vector<double>& arr, double min_value, double max_value);

void calculate_errors(std::vector<double> &source_areas, std::vector<double> &target_areas, std::vector<std::vector<std::vector<double>>> cells, std::vector<double> &errors);

// Vector math
std::vector<std::vector<std::vector<double>>> calculate_gradient(const std::vector<std::vector<double>>& grid);
std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny);

// SVG export
void export_cells_as_svg(std::vector<std::vector<std::vector<double>>> cells, std::vector<double> intensities, std::string filename);
void export_grid_to_svg(std::vector<std::vector<double>> &points, double width, double height, int res_x, int res_y, std::string filename, double stroke_width);
void export_triangles_to_svg(std::vector<std::vector<double>> &points, std::vector<std::vector<int>> &triangles, double width, double height, int res_x, int res_y, std::string filename, double stroke_width);

// 3d export
void save_solid_obj(std::vector<std::vector<double>> &front_points, std::vector<std::vector<double>> &back_points, std::vector<std::vector<int>> &triangles, double thickness, double width, double height, int res_x, int res_y, const std::string& filename);

// linear algebra
std::vector<double> vector_subtract(const std::vector<double>& p1, const std::vector<double>& p2);
std::vector<double> cross_product(const std::vector<double>& p1, const std::vector<double>& p2);
double dot_product(const std::vector<double>& p1, const std::vector<double>& p2);

std::vector<double> normalize(std::vector<double> p1);

void calculate_angle_and_normal_from_triangle(std::vector<double> &p1, std::vector<double> &p2,std::vector<double> &p3, std::vector<double> &normal_out, double &angle_out);

#endif