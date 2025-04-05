#include "caustic_design.h"

Caustic_design::Caustic_design(/* args */)
{
    this->mesh_res_x = 0;
    this->mesh_res_y = 0;
    this->resolution_x = 0;
    this->resolution_y = 0;
    this->width = 0.0f;
    this->height = 0.0f;
    this->focal_l = 0.0f;
    this->thickness = 0.0f;
    this->nthreads = 0;
}

Caustic_design::~Caustic_design()
{
}

void Caustic_design::export_inverted_transport_map(std::string filename, double stroke_width) {
    mesh->calculate_and_export_inverted_transport_map(filename, stroke_width);
}

void Caustic_design::export_paramererization_to_svg(const std::string& filename, double line_width) {
    mesh->export_paramererization_to_svg(filename, line_width);
}

void Caustic_design::set_mesh_resolution(int width, int height) {
    this->mesh_res_x = width;
    this->mesh_res_y = height;
}

void Caustic_design::set_domain_resolution(int width, int height) {
    this->resolution_x = width;
    this->resolution_y = height;
}

void Caustic_design::set_mesh_size(double width, double height) {
    this->width = width;
    this->height = height;
}

void Caustic_design::set_lens_focal_length(double focal_length) {
    this->focal_l = focal_length;
}

void Caustic_design::set_lens_thickness(double thickness) {
    this->thickness = thickness;
}

void Caustic_design::set_solver_max_threads(int n_threads) {
    this->nthreads = n_threads;
}

void Caustic_design::save_solid_obj_target(const std::string& filename) {
    this->mesh->save_solid_obj_target(thickness, filename);
}

void Caustic_design::save_solid_obj_source(const std::string& filename) {
    this->mesh->save_solid_obj_source(thickness, filename);
}

void clamp(int &value, int min, int max) {
    value = std::max(std::min(value, max), min);
}

// Bilinear interpolation function
double bilinearInterpolation(const std::vector<std::vector<double>>& image, double x, double y) {
    int x0 = floor(x);
    int y0 = floor(y);
    int x1 = ceil(x);
    int y1 = ceil(y);

    clamp(x0, 0, image[0].size() - 1);
    clamp(x1, 0, image[0].size() - 1);
    clamp(y0, 0, image.size() - 1);
    clamp(y1, 0, image.size() - 1);

    // Check if the point is outside the image bounds
    if (x0 < 0 || y0 < 0 || x1 >= image[0].size() || y1 >= image.size()) {
        printf("interpolation out of range: x: %f, y: %f\r\n", x, y);

        printf("x0: %i, y0: %i, x1: %i, y1: %i\r\n", x0, y0, x1, y1);
        // Handle out-of-bounds condition
        return 0.0;  // Default value
    }

    // Interpolate along x-axis
    double fx1 = x - x0;
    double fx0 = 1.0 - fx1;

    // Interpolate along y-axis
    double fy1 = y - y0;
    double fy0 = 1.0 - fy1;

    // Perform bilinear interpolation
    double top = fx0 * image[y0][x0] + fx1 * image[y0][x1];
    double bottom = fx0 * image[y1][x0] + fx1 * image[y1][x1];
    return fy0 * top + fy1 * bottom;
}

std::vector<std::vector<std::vector<double>>> calculate_second_derivatives(const std::vector<std::vector<double>>& grid) {
    int width = static_cast<int>(grid[0].size());
    int height = static_cast<int>(grid.size());

    std::vector<std::vector<double>> u_xx(height, std::vector<double>(width, 0.0));
    std::vector<std::vector<double>> u_yy(height, std::vector<double>(width, 0.0));
    std::vector<std::vector<double>> u_xy(height, std::vector<double>(width, 0.0));
    std::vector<std::vector<double>> u_x(height, std::vector<double>(width, 0.0));
    std::vector<std::vector<double>> u_y(height, std::vector<double>(width, 0.0));

    // Correct grid spacing calculation
    double h_x = 2.0 / (width);
    double h_y = 2.0 / (height);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Central differences with boundary clamping
            int x_prev = std::max(x - 1, 0);
            int x_next = std::min(x + 1, width - 1);
            int y_prev = std::max(y - 1, 0);
            int y_next = std::min(y + 1, height - 1);

            u_xx[y][x] = -(grid[y][x_next] + grid[y][x_prev] - 2 * grid[y][x]) / (h_x * h_x);
            u_yy[y][x] = -(grid[y_next][x] + grid[y_prev][x] - 2 * grid[y][x]) / (h_y * h_y);
            u_xy[y][x] = -(grid[y_next][x_next] + grid[y_prev][x_prev] - grid[y_next][x_prev] - grid[y_prev][x_next]) / (4 * h_x * h_y);
            u_x[y][x] =  -(grid[y][x_next] - grid[y][x_prev] ) / (2 * h_x);
            u_y[y][x] =  -(grid[y_next][x] - grid[y_prev][x] ) / (2 * h_y);
        }
    }

    return {u_xx, u_yy, u_xy, u_x, u_y};
}

void _normalize(std::vector<std::vector<double>>& F) {
    // Validate input matrix
    if (F.empty() || F[0].empty()) {
        throw std::invalid_argument("Input matrix F cannot be empty");
    }
    
    const int height = F.size();
    const int width = F[0].size();
    
    // Check matrix consistency
    for (const auto& row : F) {
        if (row.size() != static_cast<size_t>(width)) {
            throw std::invalid_argument("All rows in F must have the same width");
        }
    }

    // Calculate scaling factors
    const double hx = 2.0 / width;
    const double hy = 2.0 / height;

    // Calculate weighted sum
    double sum = 0.0;
    for (int j = 0; j < height; ++j) {        // Rows (y-axis)
        for (int i = 0; i < width; ++i) {      // Columns (x-axis)
            sum += F[j][i] * hx * hy;
        }
    }

    // Avoid division by zero
    if (sum == 0.0) {
        throw std::runtime_error("Cannot normalize: sum of elements is zero");
    }

    // Apply normalization factor
    const double factor = 4.0 / sum;
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            F[j][i] *= factor;
        }
    }
}

#include <png.h>

void grid_to_image(const std::vector<std::vector<double>>& image_grid, const std::string& filename) {
    if (image_grid.empty()) {
        throw std::runtime_error("Image grid is empty.");
    }
    size_t height = image_grid.size();
    size_t width = image_grid[0].size();
    if (width == 0 || height == 0) {
        throw std::runtime_error("Image grid has invalid dimensions (zero width or height).");
    }
    for (const auto& row : image_grid) {
        if (row.size() != width) {
            throw std::runtime_error("Image grid rows have inconsistent lengths.");
        }
    }

    FILE* fp = fopen(filename.c_str(), "wb");
    if (!fp) {
        throw std::runtime_error("Failed to open file for writing.");
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png) {
        fclose(fp);
        throw std::runtime_error("Failed to create PNG write struct.");
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        png_destroy_write_struct(&png, nullptr);
        fclose(fp);
        throw std::runtime_error("Failed to create PNG info struct.");
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        throw std::runtime_error("Error during PNG creation.");
    }

    png_init_io(png, fp);

    png_set_IHDR(
        png,
        info,
        width,
        height,
        8,
        PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );

    png_write_info(png, info);

    std::vector<png_bytep> row_pointers(height);
    for (size_t y = 0; y < height; ++y) {
        row_pointers[y] = static_cast<png_bytep>(png_malloc(png, png_get_rowbytes(png, info)));
        for (size_t x = 0; x < width; ++x) {
            double gray = image_grid[y][x];
            // Clamp the grayscale value to [0.0, 1.0]
            gray = std::max(0.0, std::min(gray, 1.0));
            png_byte value = static_cast<png_byte>(gray * 255.0);
            png_bytep pixel = &(row_pointers[y][x * 4]);
            pixel[0] = value; // Red
            pixel[1] = value; // Green
            pixel[2] = value; // Blue
            pixel[3] = 255;   // Alpha (fully opaque)
        }
    }

    png_write_image(png, row_pointers.data());
    png_write_end(png, nullptr);

    // Cleanup
    for (size_t y = 0; y < height; ++y) {
        png_free(png, row_pointers[y]);
    }
    png_destroy_write_struct(&png, &info);
    fclose(fp);
}

std::vector<std::vector<double>> Caustic_design::calculate_updated_distribution(std::vector<std::vector<double>> image, std::vector<std::vector<double>> u_x, std::vector<std::vector<double>> u_y) {
    // bilinear interpolating of the gradients into the vertices
    std::vector<double> vertex_gradient_x;
    std::vector<double> vertex_gradient_y;
    for (int i=0; i<mesh->target_points.size(); i++) {
        vertex_gradient_x.push_back(bilinearInterpolation(u_x, 
            (mesh->target_points[i][0] / mesh->width) * (resolution_x) - 0.5, 
            (mesh->target_points[i][1] / mesh->height) * (resolution_y) - 0.5
        ));

        vertex_gradient_y.push_back(bilinearInterpolation(u_y, 
            (mesh->target_points[i][0] / mesh->width) * (resolution_x) - 0.5, 
            (mesh->target_points[i][1] / mesh->height) * (resolution_y) - 0.5
        ));
    }
    
    // compute barycentric cells based on refrence triangular mesh
    std::vector<std::vector<polygon_t>> start_cells;
    for (int i=0; i<mesh->target_points.size(); i++) {
        std::vector<polygon_t> cell = mesh->get_partitioned_barycentric_dual_cell(i, mesh->target_points);
        start_cells.push_back(cell);
    }

    // get areas of the start cells
    std::vector<double> start_areas = get_partitioned_source_areas(start_cells);

    // integrate the pixels into the dual cells to get vertex intensity values
    std::vector<double> vertex_values;
    for (int i=0; i<start_cells.size(); i++) {
		double total_cell_area = 0.0f;
		for (int j=0; j<start_cells[i].size(); j++) {
			total_cell_area += integrate_cell_intensities(image, start_cells[i][j], resolution_x, resolution_y, width);
		}
		vertex_values.push_back(total_cell_area / start_areas[i]);
	}

    /*std::vector<std::vector<double>> velocities;
    for (size_t i = 0; i < mesh->target_points.size(); ++i) {
        std::vector<double> velocity;
        velocity.push_back(vertex_gradient_x[i]);
        velocity.push_back(vertex_gradient_y[i]);

        velocities.push_back(velocity);
    }

    double min_t = mesh->find_min_delta_t(velocities);

    std::cout << min_t << std::endl;*/

    // move the mesh along the gradient of the potential
    std::vector<point_t> moved_points;
    for (int i = 0; i < mesh->target_points.size(); i++) {
        int y = i / mesh->res_x;
        int x = i % mesh->res_x;
        double dx = vertex_gradient_x[i];
        double dy = vertex_gradient_y[i];

        // Apply boundary conditions to lock points
        if ((x == 0 && y == 0) || 
            (x == 0 && y == mesh->res_y - 1) || 
            (x == mesh->res_x - 1 && y == 0) || 
            (x == mesh->res_x - 1 && y == mesh->res_y - 1)) {
            // Corners: no movement
            dx = 0;
            dy = 0;
        } else if (x == 0 || x == mesh->res_x - 1) {
            // Left/Right edges: only move in y-direction
            dx = 0;
        } else if (y == 0 || y == mesh->res_y - 1) {
            // Top/Bottom edges: only move in x-direction
            dy = 0;
        }

        moved_points.push_back({
            mesh->target_points[i][0] + dx*0.3,
            mesh->target_points[i][1] + dy*0.3
        });
    }

    export_grid_to_svg(moved_points, mesh->width, mesh->height, mesh->res_x, mesh->res_y, "moved_grid.svg", 1.0);

    // compute barycentric cells of the moved cells
    std::vector<std::vector<polygon_t>> moved_cells;
    for (int i=0; i<moved_points.size(); i++) {
        moved_cells.push_back(mesh->get_partitioned_barycentric_dual_cell(i, moved_points));
    }

    // get areas of the moved cells
    std::vector<double> moved_areas = get_partitioned_source_areas(moved_cells);

    // multiply the vertex values by the jacobian determinant
    //for (int i = 0; i < moved_cells.size(); i++)
    //{
    //    vertex_values[i] *= (start_areas[i] / moved_areas[i])*0.5;
    //}

    export_cells_as_svg(moved_cells, scale_array_proportional(vertex_values, 0.0f, 1.0f), "../cells.svg");

    // rasterize the moved mesh into a grid image 
    bool triangle_miss = false;
    std::vector<std::vector<double>> interpolation = mesh->interpolate_raster(vertex_values, moved_points, mesh->triangles, resolution_x, resolution_y, mesh->width, mesh->height, triangle_miss);//*/

    return interpolation;
}

#include <algorithm>

double Caustic_design::perform_transport_iteration() {
    std::vector<std::vector<double>> potential;

    for (int i = 0; i < resolution_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_x; ++j) {
            row.push_back(0.0f);
        }
        potential.push_back(row);
    }

    std::vector<std::vector<std::vector<double>>> derivatives = calculate_second_derivatives(kantorovich_potential);

    std::vector<std::vector<double>> u_xx = derivatives[0];
    std::vector<std::vector<double>> u_yy = derivatives[1];
    std::vector<std::vector<double>> u_xy = derivatives[2];
    std::vector<std::vector<double>> u_x = derivatives[3];
    std::vector<std::vector<double>> u_y = derivatives[4];

    std::vector<std::vector<double>> updated_source = calculate_updated_distribution(pixels_trg, u_x, u_y);

    _normalize(updated_source);

    std::vector<std::vector<double>> rhs;
    for (int i = 0; i < resolution_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_x; ++j) {
            if (updated_source[i][j] < 1e-10) {
                updated_source[i][j] = 1e-10;
            }

            double dxx = u_xx[i][j] + 1.0;
            double dyy = u_yy[i][j] + 1.0;
            double dxy = u_xy[i][j];
            double f = (pixels[i][j]) / (updated_source[i][j]);

            double d = pow(dxx, 2) + pow(dyy, 2) + 2.0 * pow(dxy, 2) + 2.0 * f;
            d = sqrt(d) - 2.0;

            row.push_back(d);
        }
        rhs.push_back(row);
    }

    grid_to_image(scale_matrix_proportional(updated_source, 0.0f, 1.0f), "interpolation.png");
    grid_to_image(scale_matrix_proportional(pixels, 0.0f, 1.0f), "pixels_trg.png");
    grid_to_image(scale_matrix_proportional(rhs, 0.0f, 1.0f), "rhs.png");

    // Calculate scaling factors
    const double hx = 2.0 / resolution_x;
    const double hy = 2.0 / resolution_y;

    subtractAverage(rhs);
    poisson_solver(rhs, potential, hx, hy, 1000000, 1.0E-7, nthreads);

    for (int i = 0; i < resolution_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_x; ++j) {
            kantorovich_potential[i][j] = potential[i][j];
            //kantorovich_potential[i][j] *= 0.5;
        }
    }
    //*/

    return 0.0f;

    //std::vector<std::vector<double>> vertex_gradient;
    /*double min_step = 0.0f;

    // build median dual mesh of the updated parameterization
    target_cells.clear();
    mesh->build_target_partitioned_dual_cells(target_cells);
    std::vector<double> source_areas = get_partitioned_source_areas(target_cells);

    errors.clear();
    
    // calculate difference D (equation 2 in the paper)
    for (int i=0; i<target_areas.size(); i++) {
        errors.push_back(target_areas[i] - source_areas[i]);
    }

    // scale errors by the inverse cell area. Keeps error magnitude consistent
    for (int i=0; i<target_areas.size(); i++) {
        errors[i] = errors[i] / calculate_partitioned_cell_area(target_cells[i]);
    }

    // rasterize the mesh into a uniform rectangular matrix
    bool triangle_miss = false;
    raster = mesh->interpolate_raster_target(errors, resolution_x, resolution_y, triangle_miss);
    
    if (triangle_miss) {
        mesh->laplacian_smoothing(mesh->target_points, 0.1f);
        return NAN;
    }

    // solve the poisson equation 3 in the paper
    subtractAverage(raster);
    poisson_solver(raster, phi, resolution_x, resolution_y, 100000, 0.0000001, nthreads);

    // calculate the gradient given by equation 4
    gradient = calculate_gradient(phi);

    // bilinear interpolating the gradients
    std::vector<double> vertex_gradient_x;
    std::vector<double> vertex_gradient_y;
    for (int i=0; i<mesh->target_points.size(); i++) {
        vertex_gradient_x.push_back(bilinearInterpolation(gradient[0], 
            (mesh->target_points[i][0] / mesh->width) * (resolution_x) - 0.5, 
            (mesh->target_points[i][1] / mesh->height) * (resolution_y) - 0.5
        ));

        vertex_gradient_y.push_back(bilinearInterpolation(gradient[1], 
            (mesh->target_points[i][0] / mesh->width) * (resolution_x) - 0.5, 
            (mesh->target_points[i][1] / mesh->height) * (resolution_y) - 0.5
        ));
    }

    vertex_gradient.clear();
    vertex_gradient.push_back(vertex_gradient_x);
    vertex_gradient.push_back(vertex_gradient_y);
    
    std::vector<std::vector<double>> old_points;
    std::copy(mesh->target_points.begin(), mesh->target_points.end(), back_inserter(old_points));

    // step the mesh vertices in the direction of their gradient vector
    mesh->step_grid(vertex_gradient[0], vertex_gradient[1], 0.05f);

    // calculate the mesh movement size for convergence status
    for (int i=0; i<old_points.size(); i++) {
        double dx = (old_points[i][0] - mesh->target_points[i][0]);
        double dy = (old_points[i][1] - mesh->target_points[i][1]);
        double dz = (old_points[i][2] - mesh->target_points[i][2]);

        double dist = sqrt(dx*dx + dy*dy + dz*dz);

        if (min_step < dist) {
            min_step = dist;
        }
    }

    return min_step / width;*/
}

// uses uniform grid as caustic surface
void Caustic_design::perform_height_map_iteration(int itr) {
    // calculate the target normals
    normals = mesh->calculate_refractive_normals_uniform(resolution_x / width * focal_l, 1.49);

    // interpolates the vertex normals into a large uniform grid
    mesh->build_source_bvh(5, 30);
    bool triangle_miss = false;
    norm_x = mesh->interpolate_raster_source(normals[0], resolution_x, resolution_y, triangle_miss);
    norm_y = mesh->interpolate_raster_source(normals[1], resolution_x, resolution_y, triangle_miss);

    if (triangle_miss) {
        return;
    }

    /*std::vector<std::vector<double>> curl = calculate_curl(norm_x, norm_y, resolution_x, resolution_y);
    std::vector<double> vertex_curl;
    for (int i=0; i<mesh->source_points.size(); i++) {
        vertex_curl.push_back(bilinearInterpolation(curl,
            (mesh->source_points[i][0] / mesh->width) * (resolution_x) - 0.5, 
            (mesh->source_points[i][1] / mesh->height) * (resolution_y) - 0.5
        ));
    }

    std::vector<polygon_t> cells;
    std::vector<double> colors;

    for (int i = 0; i < first_target_cells.size(); i++)
    {
        for (int j = 0; j < first_target_cells[i].size(); j++)
        {
            cells.push_back(first_target_cells[i][j]);
            colors.push_back(vertex_curl[i]);
        }
    }

    export_cells_as_svg(cells, scale_array_proportional(colors, 0.0f, 1.0f), "../curl.svg");*/

    // calculates the divergance of the interpolated normals
    divergence = calculate_divergence(norm_x, norm_y, resolution_x, resolution_y);
    subtractAverage(divergence);

    // solve the poisson equation for the divergance
    poisson_solver(divergence, h, resolution_x, resolution_y, 100000, 0.00000001, nthreads);

    /*std::vector<double> interpolated_h;
    for (int i=0; i<mesh->target_points.size(); i++) {
        interpolated_h.push_back(bilinearInterpolation(h, mesh->target_points[i][0] * ((resolution_x) / mesh->width), mesh->target_points[i][1] * ((resolution_y) / mesh->height)));
    }
    double max_update = mesh->set_target_heights(interpolated_h);
    printf("height max update %.5e\r\n", max_update);*/

    double epsilon = 1e-8;

    // get the heights on the vertex positions
    std::vector<double> interpolated_h;
    for (int i=0; i<mesh->source_points.size(); i++) {
        interpolated_h.push_back(bilinearInterpolation(h, 
            (mesh->source_points[i][0] / mesh->width) * (resolution_x) - 0.5, 
            (mesh->source_points[i][1] / mesh->height) * (resolution_y) - 0.5
        ));
    }
    double max_update = mesh->set_source_heights(interpolated_h);
    printf("height max update %.5e\r\n", max_update);
}

void Caustic_design::initialize_solvers(std::vector<std::vector<double>> image, std::vector<std::vector<double>> target) {
    pixels = scale_matrix_proportional(image, 0, 1.0f);
    pixels_trg = scale_matrix_proportional(target, 0, 1.0f);

    _normalize(pixels);
    _normalize(pixels_trg);

    //printf("scaled\r\n");

    mesh = new Mesh(width, height, mesh_res_x, mesh_res_y);

    std::cout << "built mesh" << std::endl;

    target_cells.clear();
    mesh->build_target_partitioned_dual_cells(target_cells);

    first_target_cells = target_cells;

    target_areas = get_target_partitioned_areas(pixels, target_cells, resolution_x, resolution_y, width, height);

    std::cout << target_areas.size() << std::endl;

    //export_cells_as_svg(target_cells, scale_array_proportional(target_areas, 0.0f, 1.0f), "../cells.svg");

    phi.clear();
    h.clear();
    kantorovich_potential.clear();
    for (int i = 0; i < resolution_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_x; ++j) {
            row.push_back(0.0f);
        }
        phi.push_back(row);
        h.push_back(row);
        kantorovich_potential.push_back(row);
    }
}