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

    if (rows == 0) {
        throw std::invalid_argument("Input matrix is empty.");
    }

    size_t cols = matrix[0].size();

    // Check if all inner vectors have the same size
    for (size_t i = 1; i < rows; ++i) {
        if (matrix[i].size() != cols) {
            throw std::invalid_argument("Input matrix has inconsistent row sizes.");
        }
    }

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

/*std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny) {
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
}*/

std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny) {
    std::vector<std::vector<double>> divergence(ny, std::vector<double>(nx, 0.0));

    std::vector<std::vector<std::vector<double>>> grad_x = calculate_gradient(Nx);
    std::vector<std::vector<std::vector<double>>> grad_y = calculate_gradient(Ny);

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            if (x == 0 || x == nx - 1 || y == 0 || y == ny - 1) {
                divergence[y][x] = 0.0; // Set divergence to 0 for boundary cells
            } else {
                divergence[y][x] = grad_x[0][y][x] + grad_y[1][y][x]; // Calculate and set divergence
            }
        }
    }

    return divergence;
}//*/

/*std::vector<std::vector<double>> calculate_divergence(const std::vector<std::vector<double>>& Nx, const std::vector<std::vector<double>>& Ny, int nx, int ny) {
    int height = ny;
    int width = nx;

    std::vector<std::vector<double>> divergence(height, std::vector<double>(width, 0.0));

    // Iterate over each point in the grid
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Calculate divergence at each point using gradient components Nx and Ny
            double div_x = (x < width - 1) ? Nx[y][x + 1] - Nx[y][x] : 0.0;
            double div_y = (y < height - 1) ? Ny[y + 1][x] - Ny[y][x] : 0.0;
            
            // Calculate the total divergence
            divergence[y][x] = div_x + div_y;
        }
    }

    return divergence;
}//*/

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

void find_perimeter_vertices(int nx, int ny, std::vector<int> &perimeter_vertices) {
    // Top row
    for (int i = 0; i < nx; ++i) {
        perimeter_vertices.push_back(i);
    }

    // Right column
    for (int i = nx - 1; i < nx * ny; i += nx) {
        perimeter_vertices.push_back(i);
    }

    // Bottom row
    for (int i = nx * (ny - 1) + nx - 1; i > nx * (ny - 1) - 1; --i) {
        perimeter_vertices.push_back(i);
    }

    // Left column
    for (int i = nx * (ny - 1) - nx; i > nx - 1; i -= nx) {
        perimeter_vertices.push_back(i);
    }
}

void save_solid_obj(std::vector<std::vector<double>> &front_points, std::vector<std::vector<double>> &back_points, std::vector<std::vector<int>> &triangles, double thickness, double width, double height, int res_x, int res_y, const std::string& filename) {
    int num_points = front_points.size();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Find maximum height
    double min_h = 0;
    for (int i=0; i<num_points; i++) {
        double h = front_points[i][2];

        if (min_h > h) {
            min_h = h;
        }
    }

    file << "# Generated by the software algorithm written by Dylan Missuwe" << "\n";
    file << "# The algorithm used to create the lens is based on the paper " 
        << "\"Poisson-Based Continuous Surface Generation for Goal-Based Caustics\" " 
        << "by Yue et al (2014)" << "\n";

    // Curved mesh verts on the bottom
    for (const auto& point : front_points) {
        file << "v " << width - point[0] << " " << height - point[1] << " " << -(point[2]) << "\n";
    }

    // Flat mesh verts on the bottom
    for (const auto& point : back_points) {
        file << "v " << width - point[0] << " " << height - point[1] << " " << -thickness - min_h << "\n";
    }

    // Curved mesh triangles on the top
    for (const auto& triangle : triangles) {
        file << "f " << triangle[0] + 1 << " " << triangle[1] + 1 << " " << triangle[2] + 1 << "\n";
    }

    // Flat mesh triangles on the bottom
    for (const auto& triangle : triangles) {
        file << "f " << triangle[0] + num_points + 1 << " " << triangle[2] + num_points + 1 << " " << triangle[1] + num_points + 1 << "\n";
    }

    // Generate triangles connecting top and bottom mesh
    std::vector<int> perimeter_verts;
    find_perimeter_vertices(res_x, res_y, perimeter_verts);

    for (size_t i = 0; i < perimeter_verts.size(); ++i) {
        int top_idx = perimeter_verts[i];
        int bottom_idx = perimeter_verts[i] + num_points;
        int next_top_idx = perimeter_verts[(i + 1) % perimeter_verts.size()];
        int next_bottom_idx = perimeter_verts[(i + 1) % perimeter_verts.size()] + num_points;


        file << "f " << top_idx + 1 << " " << bottom_idx + 1 << " " << next_bottom_idx + 1 << "\n";
        file << "f " << top_idx + 1 << " " << next_bottom_idx + 1 << " " << next_top_idx + 1 << "\n";
        //file << "f " << top_idx + 1 << " " << next_bottom_idx + 1 << " " << bottom_idx + 1 << "\n";
        //file << "f " << top_idx + 1 << " " << next_top_idx + 1 << " " << next_bottom_idx + 1 << "\n";
    }

    //std::cout << "Exported model \"" << filename << "\"" << std::endl;
}

void export_grid_to_svg(std::vector<std::vector<double>> &points, double width, double height, int res_x, int res_y, std::string filename, double stroke_width) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * (height / width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    for (int j = 0; j < res_y; j++) {
        std::string path_str = "M";
        for (int i = 0; i < res_x; i++) {
            int idx = i + j * res_x;

            const auto& point = points[idx];
            path_str += std::to_string((point[0] / width) * 1000.0f) + "," +
                        std::to_string((point[1] / height) * 1000.0f * (height / width));

            if (i < res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    for (int j = 0; j < res_x; j++) {
        std::string path_str = "M";
        for (int i = 0; i < res_y; i++) {
            int idx = j + i * res_x;

            const auto& point = points[idx];
            path_str += std::to_string((point[0] / width) * 1000.0f) + "," +
                        std::to_string((point[1] / height) * 1000.0f * (height / width));

            if (i < res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

void export_triangles_to_svg(std::vector<std::vector<double>> &points, std::vector<std::vector<int>> &triangles, double width, double height, int res_x, int res_y, std::string filename, double stroke_width) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * (height / width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
    
    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Draw polygons
    for (const auto& polygon : triangles) {
        std::vector<std::vector<double>> poly_points;
        for (const auto& point_idx : polygon) {
            poly_points.push_back(points[point_idx]);
        }

        std::string path_str = "M";
        for (std::size_t j = 0; j < poly_points.size(); ++j) {
            const auto& point = poly_points[j];
            path_str += std::to_string((point[0] / width) * 1000.0f) + "," +
                        std::to_string((point[1] / height) * 1000.0f * (height / width));

            if (j < poly_points.size() - 1)
                path_str += "L";
        }
        path_str += "Z";
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

std::vector<double> vector_subtract(const std::vector<double>& p1, const std::vector<double>& p2) {
    std::vector<double> difference(p1.size());
    
    for (size_t i = 0; i < p1.size(); ++i) {
        difference[i] = p1[i] - p2[i];
    }

    return difference;
}

std::vector<double> cross_product(const std::vector<double>& p1, const std::vector<double>& p2) {
    std::vector<double> cross(3);
    
    cross[0] = p1[1]*p2[2] - p1[2]*p2[1];
    cross[1] = p1[2]*p2[0] - p1[0]*p2[2];
    cross[2] = p1[0]*p2[1] - p1[1]*p2[0];

    return cross;
}

double dot_product(const std::vector<double>& p1, const std::vector<double>& p2) {
    return p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
}

std::vector<double> normalize(std::vector<double> p1) {
    std::vector<double> vec(3);
    double squared_len = 0;
    for (int i=0; i<p1.size(); i++) {
        squared_len += p1[i] * p1[i];
    }

    double len = std::sqrt(squared_len);

    for (int i=0; i<p1.size(); i++) {
        vec[i] = p1[i] / len;
    }

    return vec;
}

void calculate_angle_and_normal_from_triangle(std::vector<double> &p1, std::vector<double> &p2,std::vector<double> &p3, std::vector<double> &normal_out, double &angle_out) {
    std::vector<double> u = vector_subtract(p2, p1);
    std::vector<double> v = vector_subtract(p3, p1);
    
    normal_out = normalize(cross_product(u, v));

    std::vector<double> nu = normalize(u);
    std::vector<double> nv = normalize(v);

    double res = dot_product(nu, nv);

    // Ensure res is within valid range [-1, 1]
    if (res <= -1.0) {
        angle_out = M_PI; // angle is 180 degrees
    } else if (res >= 1.0) {
        angle_out = 0.0; // angle is 0 degrees
    } else {
        angle_out = acos(res);
    }
}