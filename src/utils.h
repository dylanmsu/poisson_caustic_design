#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "polygon_utils.h"

std::vector<std::vector<std::vector<double>>> calculate_gradient(const std::vector<std::vector<double>>& grid);
void subtractAverage(std::vector<std::vector<double>>& raster);
std::vector<std::vector<double>> scale_matrix_proportional(const std::vector<std::vector<double>>& matrix, double min_value, double max_value);
std::vector<double> scale_array_proportional(const std::vector<double>& arr, double min_value, double max_value);
void export_cells_as_svg(std::vector<std::vector<std::vector<double>>> cells, std::vector<double> intensities, std::string filename);
std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny);
void calculate_errors(std::vector<double> &source_areas, std::vector<double> &target_areas, std::vector<std::vector<std::vector<double>>> cells, std::vector<double> &errors);

#endif