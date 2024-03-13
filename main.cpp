#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "src/solver.h"
#include "src/mesh.h"
#include "src/polygon_utils.h"
#include "src/utils.h"

#include <thread>

#define cimg_use_png
//#define cimg_use_jpeg

#if defined(_WIN32) || defined(_WIN64)
    #include "CImg.h"
#else
    #include<X11/Xlib.h>
    #include "/home/dylan/caustic_engineering/CImg.h"
#endif

void image_to_grid(const cimg_library::CImg<unsigned char>& image, std::vector<std::vector<double>>& image_grid) {
    for (int i = 0; i < image.height(); ++i) {
        std::vector<double> row;
        for (int j = 0; j < image.width(); ++j) {
            double r = image(j, i, 0) / 255.0; // Normalize R component
            double g = image(j, i, 1) / 255.0; // Normalize G component
            double b = image(j, i, 2) / 255.0; // Normalize B component
            double gray = (0.299 * r) + (0.587 * g) + (0.114 * b); // Calculate grayscale value using luminosity method
            row.push_back(gray);
        }
        image_grid.push_back(row);
    }
}

void save_grid_as_image(const std::vector<std::vector<double>>& img, int resolution_x, int resolution_y, const std::string& filename) {
    // Create an empty CImg object with the specified resolution
    cimg_library::CImg<unsigned char> image(resolution_x, resolution_y);

    // Copy the grid data to the image
    for (int i = 0; i < resolution_y; ++i) {
        for (int j = 0; j < resolution_x; ++j) {
            image(j, i) = static_cast<unsigned char>(img[i][j] * 255); // Scale to [0, 255] and cast to unsigned char
        }
    }

    // Save the image as a PNG
    image.save(filename.c_str());
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

std::unordered_map<std::string, std::string> parse_arguments(int argc, char const *argv[]) {
    // Define a map to store the parsed arguments
    std::unordered_map<std::string, std::string> args;

    // Iterate through command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string key, value;

        // Check if argument starts with '--'
        if (arg.substr(0, 2) == "--") {
            // Split argument by '=' to separate key and value
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                key = arg.substr(2, pos - 2);
                value = arg.substr(pos + 1);
            }
            else {
                key = arg.substr(2);
                value = ""; // No value provided
            }
        }
        // Check if argument starts with '-'
        else if (arg[0] == '-') {
            // The next argument is the value
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                key = arg.substr(1);
                value = argv[++i];
            }
            else {
                key = arg.substr(1);
                value = ""; // No value provided
            }
        }
        // Invalid argument format
        else {
            //std::cerr << "Invalid argument format: " << arg << std::endl;
            //return 1;
        }

        // Store key-value pair in the map
        args[key] = value;
    }

    return args;
}

Mesh *mesh;
std::vector<std::vector<double>> phi;
std::vector<double> errors;
std::vector<std::vector<std::vector<double>>> target_cells;
std::vector<std::vector<std::vector<double>>> source_cells;
std::vector<double> target_areas;
std::vector<std::vector<double>> pixels;
std::vector<std::vector<double>> raster;
std::vector<std::vector<std::vector<double>>> gradient;

int mesh_res_x;
int mesh_res_y;

int resolution_x;
int resolution_y;

double width;
double height;

double focal_l;
double thickness;
int nthreads;

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
std::vector<double> calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index) {
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

double perform_transport_iteration() {
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
    raster = mesh->interpolate_raster(errors, resolution_x, resolution_y, triangle_miss);

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

    // step the mesh vertices in the direction of their gradient vector
    min_step = mesh->step_grid(vertex_gradient[0], vertex_gradient[1], 0.95f);

    mesh->laplacian_smoothing(mesh->target_points, min_step*(resolution_x/width) / 2);
    //mesh->laplacian_smoothing(mesh->target_points, 0.1f);

    return min_step*(resolution_x/width);
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
void perform_height_map_iteration(int itr) {
    // calculate the target normals
    std::vector<std::vector<double>> normals = mesh->calculate_refractive_normals(resolution_x / width * focal_l, 1.49);

    // interpolates the vertex normals into a large uniform grid
    mesh->build_bvh(1, 30);
    bool triangle_miss = false;
    std::vector<std::vector<double>> norm_x = mesh->interpolate_raster(normals[0], resolution_x, resolution_y, triangle_miss);
    std::vector<std::vector<double>> norm_y = mesh->interpolate_raster(normals[1], resolution_x, resolution_y, triangle_miss);

    if (triangle_miss) {
        return;
    }

    // calculates the divergance of the interpolated normals
    std::vector<std::vector<double>> divergence = calculate_divergence(norm_x, norm_y, resolution_x, resolution_y);
    subtractAverage(divergence);

    // solve the poisson equation for the divergance
    std::vector<std::vector<double>> h(resolution_y, std::vector<double> (resolution_x, 0.0));
    poisson_solver(divergence, h, resolution_x, resolution_y, 100000, 0.0000001, nthreads);

    std::vector<double> interpolated_h;
    for (int i=0; i<mesh->target_points.size(); i++) {
        interpolated_h.push_back(bilinearInterpolation(h, mesh->target_points[i][0] * ((resolution_x) / mesh->width), mesh->target_points[i][1] * ((resolution_y) / mesh->height)));
    }
    double max_update = mesh->set_target_heights(interpolated_h);
    printf("height max update %.5e\r\n", max_update);

    // get the heights on the vertex positions
    //std::vector<double> interpolated_h;
    //for (int i=0; i<mesh->source_points.size(); i++) {
    //    interpolated_h.push_back(bilinearInterpolation(h, mesh->source_points[i][0] * ((resolution_x) / mesh->width), mesh->source_points[i][1] * ((resolution_y) / mesh->height)));
    //}
    //mesh->set_source_heights(interpolated_h);
}

void initialize_solvers() {
    pixels = scale_matrix_proportional(pixels, 0, 1.0f);

    printf("scaled\r\n");

    mesh = new Mesh(width, height, mesh_res_x, mesh_res_y);

    printf("generated mesh\r\n");

    //std::cout << "built mesh" << std::endl;

    //std::vector<std::vector<std::vector<double>>> circ_target_cells;
    mesh->build_target_dual_cells(target_cells);
    mesh->build_source_dual_cells(source_cells);
    //mesh.build_circular_target_dual_cells(circ_target_cells);

    //std::vector<double> target_areas = get_target_areas(pixels, circ_target_cells, resolution_x, resolution_y, width, height);
    target_areas = get_target_areas(pixels, target_cells, resolution_x, resolution_y, width, height);
}

int main(int argc, char const *argv[])
{
    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    // Load image
    cimg_library::CImg<unsigned char> image(args["intput_png"].c_str());
    double aspect_ratio = (double)image.width() / (double)image.height();

    // Print parsed arguments
    /*for (const auto& pair : args) {
        //std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
        printf("Key: %s, Value: %s\r\n", pair.first.c_str(), pair.second.c_str());
    }*/

    mesh_res_x = atoi(args["res_w"].c_str());
    mesh_res_y = atoi(args["res_w"].c_str()) / aspect_ratio;

    resolution_x = 4*mesh_res_x;
    resolution_y = 4*mesh_res_y;

    width = std::stod(args["width"]);
    height = width / aspect_ratio;

    focal_l = std::stod(args["focal_l"]);

    thickness = std::stod(args["thickness"]);

    nthreads = atoi(args["max_threads"].c_str());

    phi.clear();
    for (int i = 0; i < resolution_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < resolution_x; ++j) {
            row.push_back(0.0f);
        }
        phi.push_back(row);
    }

    image = image.resize(resolution_x, resolution_y, -100, -100, 3); // Resize using linear interpolation

    // Convert image to grid
    image_to_grid(image, pixels);

    initialize_solvers();
    
    for (int itr=0; itr<100; itr++) {
        printf("starting iteration %i\r\n", itr);
        
        double step_size = perform_transport_iteration();

        // export dual cells as svg
        export_cells_as_svg(target_cells, scale_array_proportional(errors, 0.0f, 1.0f), "../cells.svg");
        //save_grid_as_image(scale_matrix_proportional(raster, 0, 1.0f), resolution_x, resolution_y, "../raster.png");
        //save_grid_as_image(scale_matrix_proportional(phi, 0, 1.0f), resolution_x, resolution_y, "../phi.png");
        //save_grid_as_image(scale_matrix_proportional(gradient[0], 0, 1.0f), resolution_x, resolution_y, "../gradient_x.png");
        //save_grid_as_image(scale_matrix_proportional(gradient[1], 0, 1.0f), resolution_x, resolution_y, "../gradient_y.png");

        // export parameterization
        mesh->export_paramererization_to_svg("../parameterization.svg", 1.0f);

        // export inverted transport map (can be used for dithering)
        mesh->calculate_inverted_transport_map("../inverted.svg", 1.0f);

        printf("step_size = %f\r\n", step_size);

        // convergence treshold for the parameterization
        if (step_size < 0.005) break;
    }

    printf("\033[0;32mTransport map solver done! Starting height solver.\033[0m\r\n");

    //mesh.target_points = mesh.circular_transform(mesh.target_points);
    //mesh.source_points = mesh.circular_transform(mesh.source_points);

    for (int itr=0; itr<10; itr++) {
        perform_height_map_iteration(itr);
    }

    printf("Height solver done! Exporting as solidified obj\r\n");

    //mesh->save_solid_obj_source(thickness, "../output.obj");
    mesh->save_solid_obj_target(thickness, "../output.obj");
    
    //*/

    return 0;
}
