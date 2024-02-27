#include "utils.h"

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
            if (!std::isnan(value)) {
                sum += value;
                count++;
            }
        }
    }

    double average = sum / count;

    // Subtract the average from each element of the raster
    for (auto& row : raster) {
        std::transform(row.begin(), row.end(), row.begin(), [average](double value) {
            if (std::isnan(value)) {
                return value;
            } else {
                return value - average;
            }
        });
    }
}

std::vector<std::vector<double>> scale_matrix_proportional(const std::vector<std::vector<double>>& matrix, double min_value, double max_value) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    // Find the min and max values in the matrix
    double matrix_min = matrix[0][0];
    double matrix_max = matrix[0][0];

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (!std::isnan(matrix[i][j])) {
                matrix_min = std::min(matrix_min, matrix[i][j]);
                matrix_max = std::max(matrix_max, matrix[i][j]);
            }
        }
    }

    // Perform proportional scaling
    std::vector<std::vector<double>> scaled_matrix(rows, std::vector<double>(cols));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (!std::isnan(matrix[i][j])) {
                scaled_matrix[i][j] = min_value + (max_value - min_value) * (matrix[i][j] - matrix_min) / (matrix_max - matrix_min);
            }
        }
    }

    return scaled_matrix;
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

void export_cells_as_svg(std::vector<std::vector<std::vector<double>>> cells, std::vector<double> intensities, std::string filename) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)1 / (double)1) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    for (int i=0; i<cells.size(); i++) {
        std::vector<std::vector<double>> cell = cells[i];
        std::string path_str = "M";
        for (std::size_t j = 0; j < cell.size(); ++j) {
            const auto& point = cell[j];
            path_str += std::to_string((point[0] / (double)1) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)1) * 1000.0f * ((double)1 / (double)1));

            if (j < cell.size() - 1)
                path_str += "L";
        }
        path_str += "Z";
        svg_file << "<path d=\"" << path_str << "\" fill=\"" << "rgb(" << intensities[i]*255 << ", " << intensities[i]*255 << ", " << intensities[i]*255 << ")\" stroke=\"none\" stroke-width=\"" << 1.0 << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny) {
    std::vector<std::vector<double>> divergence(ny, std::vector<double>(nx, 0.0));

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            // Forward and backward differences for divergence with boundary checks
            double div_x = 0.0;
            double div_y = 0.0;

            // Check if neighboring cells are valid (not NaN) before calculating differences
            if (x < nx - 1 && !std::isnan(Nx[y][x]) && !std::isnan(Nx[y][x + 1])) {
                div_x = Nx[y][x + 1] - Nx[y][x];
            }

            if (y < ny - 1 && !std::isnan(Ny[y][x]) && !std::isnan(Ny[y + 1][x])) {
                div_y = Ny[y + 1][x] - Ny[y][x];
            }

            // Accumulate the divergence
            if (std::isnan(Nx[y][x]) || std::isnan(Ny[y][x]) || std::isnan(div_x) || std::isnan(div_y)) {
                divergence[y][x] = std::numeric_limits<double>::quiet_NaN(); // Set divergence to NaN if any neighbor is NaN
            } else if (x == 0 || x == nx - 1 || y == 0 || y == ny - 1) {
                divergence[y][x] = 0.0; // Set divergence to 0 for boundary cells
            } else {
                divergence[y][x] = div_x + div_y; // Calculate and set divergence
            }
        }
    }

    return divergence;
}

void calculate_errors(std::vector<double> &source_areas, std::vector<double> &target_areas, std::vector<std::vector<std::vector<double>>> cells, std::vector<double> &errors) {
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
}