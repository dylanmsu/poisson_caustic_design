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

std::vector<double> calculate_normal_from_points(std::vector<double> &p0, std::vector<double> &p1,std::vector<double> &p2) {
    std::vector<double> u = vector_subtract(p1, p0);
    std::vector<double> v = vector_subtract(p2, p0);

    return cross_product(u, v);
}

// Function to calculate the approximate vertex normal
std::vector<double> Caustic_design::calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index) {
    std::vector<double> avg_normal = {0.0, 0.0, 0.0}; // Initialize normal to zero vector
    
    int left_vtx = NAN;
    int right_vtx = NAN;
    int top_vtx = NAN;
    int bot_vtx = NAN;
    
    mesh->get_vertex_neighbor_ids(vertex_index, left_vtx, right_vtx, top_vtx, bot_vtx);
    
    if (!std::isnan(left_vtx) && !std::isnan(top_vtx)) {
        std::vector<double> normal = calculate_normal_from_points(points[vertex_index], points[left_vtx], points[top_vtx]);

        avg_normal[0] += normal[0];
        avg_normal[1] += normal[1];
        avg_normal[2] += normal[2];
    }

    if (!std::isnan(left_vtx) && !std::isnan(bot_vtx)) {
        std::vector<double> normal = calculate_normal_from_points(points[vertex_index], points[bot_vtx], points[left_vtx]);

        avg_normal[0] += normal[0];
        avg_normal[1] += normal[1];
        avg_normal[2] += normal[2];
    }

    if (!std::isnan(right_vtx) && !std::isnan(bot_vtx)) {
        std::vector<double> normal = calculate_normal_from_points(points[vertex_index], points[right_vtx], points[bot_vtx]);

        avg_normal[0] += normal[0];
        avg_normal[1] += normal[1];
        avg_normal[2] += normal[2];
    }

    if (!std::isnan(right_vtx) && !std::isnan(top_vtx)) {
        std::vector<double> normal = calculate_normal_from_points(points[vertex_index], points[top_vtx], points[right_vtx]);

        avg_normal[0] += normal[0];
        avg_normal[1] += normal[1];
        avg_normal[2] += normal[2];
    }

    // Calculate magnitude
    double magnitude = sqrt(avg_normal[0] * avg_normal[0] + avg_normal[1] * avg_normal[1] + avg_normal[2] * avg_normal[2]);

    // Avoid division by zero
    if (magnitude > 1e-12) {
        avg_normal[0] /= magnitude;
        avg_normal[1] /= -magnitude;
        avg_normal[2] /= magnitude;
    }

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
    std::vector<std::vector<double>> vertex_gradient;
    double min_step;

    // build median dual mesh of the updated parameterization
    target_cells.clear();
    mesh->build_target_dual_cells(target_cells);

    // calculate difference D (interpretation of equation 2)
    std::vector<double> source_areas = get_source_areas(target_cells);
    calculate_errors(source_areas, target_areas, target_cells, errors);

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

    // calculate the gradient vectors corresponding to each vertex in the mesh

    // bilinear interpolating the gradients (negligibly faster, but gives lower contrast results)
    /*std::vector<std::vector<double>> gradient(2);
    for (int i=0; i<mesh.target_points.size(); i++) {
        gradient[0].push_back(bilinearInterpolation(grad[0], mesh.target_points[i][0] * ((resolution_x - 2) / mesh.width), mesh.target_points[i][1] * ((resolution_y - 2) / mesh.height)));
        gradient[1].push_back(bilinearInterpolation(grad[1], mesh.target_points[i][0] * ((resolution_x - 2) / mesh.width), mesh.target_points[i][1] * ((resolution_y - 2) / mesh.height)));
    }//*/

    // integrate the gradient grid into the dual cells of the vertices (slower but better contrast)
    vertex_gradient = integrate_cell_gradients(gradient, target_cells, resolution_x, resolution_y, width, height);

    std::vector<std::vector<double>> old_points;

    std::copy(mesh->target_points.begin(), mesh->target_points.end(), back_inserter(old_points));

    // step the mesh vertices in the direction of their gradient vector
    mesh->step_grid(vertex_gradient[0], vertex_gradient[1], 0.95f);

    min_step = 0.0f;

    for (int i=0; i<old_points.size(); i++) {
        double dx = (old_points[i][0] - mesh->target_points[i][0]);
        double dy = (old_points[i][1] - mesh->target_points[i][1]);
        double dz = (old_points[i][2] - mesh->target_points[i][2]);

        double dist = sqrt(dx*dx + dy*dy + dz*dz);

        if (min_step < dist) {
            min_step = dist;
        }
    }

    //mesh->laplacian_smoothing(mesh->target_points, min_step*(resolution_x/width));
    mesh->laplacian_smoothing(mesh->target_points, 0.5f);

    return min_step / width;

    //return min_step*(resolution_x/width);
}

// normalize a vector such that its length equals 1
/*void normalize(std::vector<double> &vec) {
    double squared_len = 0;
    for (int i=0; i<vec.size(); i++) {
        squared_len += vec[i] * vec[i];
    }

    double len = std::sqrt(squared_len);

    for (int i=0; i<vec.size(); i++) {
        vec[i] /= len;
    }
}*/

// uses uniform grid as caustic surface
void Caustic_design::perform_height_map_iteration(int itr) {
    // calculate the target normals
    std::vector<std::vector<double>> normals = mesh->calculate_refractive_normals_uniform(resolution_x / width * focal_l, 1.49);

    // interpolates the vertex normals into a large uniform grid
    mesh->build_source_bvh(5, 30);
    bool triangle_miss = false;
    norm_x = mesh->interpolate_raster_source(normals[0], resolution_x, resolution_y, triangle_miss);
    norm_y = mesh->interpolate_raster_source(normals[1], resolution_x, resolution_y, triangle_miss);

    if (triangle_miss) {
        return;
    }

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

    double epsilon = 1e-6;

    // get the heights on the vertex positions
    std::vector<double> interpolated_h;
    for (int i=0; i<mesh->source_points.size(); i++) {
        interpolated_h.push_back(bilinearInterpolation(h, 
            mesh->source_points[i][0] * ((resolution_x - epsilon) / mesh->width + 0.5 * epsilon), 
            mesh->source_points[i][1] * ((resolution_y - epsilon) / mesh->height + 0.5 * epsilon)));
    }
    double max_update = mesh->set_source_heights(interpolated_h);
    printf("height max update %.5e\r\n", max_update);
}

void Caustic_design::initialize_solvers(std::vector<std::vector<double>> image) {
    pixels = scale_matrix_proportional(image, 0, 1.0f);

    printf("scaled\r\n");

    mesh = new Mesh(width, height, mesh_res_x, mesh_res_y);

    printf("generated mesh\r\n");

    mesh->export_to_svg("../mesh.svg", 1);

    //std::cout << "built mesh" << std::endl;

    //std::vector<std::vector<std::vector<double>>> circ_target_cells;
    mesh->build_target_dual_cells(target_cells);
    mesh->build_source_dual_cells(source_cells);
    //mesh.build_circular_target_dual_cells(circ_target_cells);

    //std::vector<double> target_areas = get_target_areas(pixels, circ_target_cells, resolution_x, resolution_y, width, height);
    target_areas = get_target_areas(pixels, target_cells, resolution_x, resolution_y, width, height);

    export_cells_as_svg(target_cells, scale_array_proportional(target_areas, 0.0f, 1.0f), "../cells.svg");

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