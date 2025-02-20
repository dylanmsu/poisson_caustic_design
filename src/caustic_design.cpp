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

// Function to calculate the approximate vertex normal
std::vector<double> Caustic_design::calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index) {
    std::vector<double> avg_normal = {0.0, 0.0, 0.0}; // Initialize normal to zero vector
    
    int left_vtx = 0;
    int right_vtx = 0;
    int top_vtx = 0;
    int bot_vtx = 0;

    //printf("aa\r\n");
    
    mesh->get_vertex_neighbor_ids(vertex_index, left_vtx, right_vtx, top_vtx, bot_vtx);

    //printf("ab\r\n");

    if (left_vtx != -1 && top_vtx != -1) {
        std::vector<double> normal;
        double angle_out;

        calculate_angle_and_normal_from_triangle(points[vertex_index], points[left_vtx], points[top_vtx], normal, angle_out);

        //printf("ac1\r\n"); fflush(stdout);
        //printf("normal = {%f, %f, %f}, angle = %f\r\n", normal[0], normal[1], normal[2], angle_out); fflush(stdout);

        avg_normal[0] += normal[0] * angle_out;
        avg_normal[1] += normal[1] * angle_out;
        avg_normal[2] += normal[2] * angle_out;
    }

    if (left_vtx != -1 && bot_vtx != -1) {
        std::vector<double> normal;
        double angle_out;

        calculate_angle_and_normal_from_triangle(points[vertex_index], points[bot_vtx], points[left_vtx], normal, angle_out);

        //printf("ac2\r\n"); fflush(stdout);
        //printf("normal = {%f, %f, %f}, angle = %f\r\n", normal[0], normal[1], normal[2], angle_out); fflush(stdout);

        avg_normal[0] += normal[0] * angle_out;
        avg_normal[1] += normal[1] * angle_out;
        avg_normal[2] += normal[2] * angle_out;
    }

    if (right_vtx != -1 && bot_vtx != -1) {
        std::vector<double> normal;
        double angle_out;

        calculate_angle_and_normal_from_triangle(points[vertex_index], points[right_vtx], points[bot_vtx], normal, angle_out);

        //printf("ac3\r\n"); fflush(stdout);
        //printf("normal = {%f, %f, %f}, angle = %f\r\n", normal[0], normal[1], normal[2], angle_out); fflush(stdout);
        
        avg_normal[0] += normal[0] * angle_out;
        avg_normal[1] += normal[1] * angle_out;
        avg_normal[2] += normal[2] * angle_out;
    }

    if (right_vtx != -1 && top_vtx != -1) {
        std::vector<double> normal;
        double angle_out;

        calculate_angle_and_normal_from_triangle(points[vertex_index], points[top_vtx], points[right_vtx], normal, angle_out);

        //printf("ac4\r\n"); fflush(stdout);
        //printf("normal = {%f, %f, %f}, angle = %f\r\n", normal[0], normal[1], normal[2], angle_out); fflush(stdout);
        
        avg_normal[0] += normal[0] * angle_out;
        avg_normal[1] += normal[1] * angle_out;
        avg_normal[2] += normal[2] * angle_out;
    }

    //printf("ad\r\n"); fflush(stdout);

    // Calculate magnitude
    double magnitude = sqrt(avg_normal[0] * avg_normal[0] + avg_normal[1] * avg_normal[1] + avg_normal[2] * avg_normal[2]);

    // Avoid division by zero
    if (magnitude > 1e-12) {
        avg_normal[0] /= -magnitude;
        avg_normal[1] /= -magnitude;
        avg_normal[2] /= magnitude;
    }

    //printf("ae\r\n"); fflush(stdout);

    return avg_normal;
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

double Caustic_design::perform_transport_iteration() {
    //std::vector<std::vector<double>> vertex_gradient;
    double min_step = 0.0f;

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

    return min_step / width;
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

void Caustic_design::initialize_solvers(std::vector<std::vector<double>> image) {
    pixels = scale_matrix_proportional(image, 0, 1.0f);

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
    for (int i = 0; i < resolution_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_x; ++j) {
            row.push_back(0.0f);
        }
        phi.push_back(row);
        h.push_back(row);
    }
}