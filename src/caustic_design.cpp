#include "caustic_design.h"

Caustic_design::Caustic_design(int width, int height, int mesh_res_x, int mesh_res_y) : mesh(width, height, mesh_res_x, mesh_res_y)
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
    mesh.calculate_and_export_inverted_transport_map(filename, stroke_width);
}

void Caustic_design::export_paramererization_to_svg(const std::string& filename, double line_width) {
    mesh.export_paramererization_to_svg(filename, line_width);
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
    this->mesh.save_solid_obj_target(thickness, filename);
}

void Caustic_design::save_solid_obj_source(const std::string& filename) {
    this->mesh.save_solid_obj_source(thickness, filename);
}

// Function to calculate the approximate vertex normal
std::vector<double> Caustic_design::calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index) {
    std::vector<double> avg_normal = {0.0, 0.0, 0.0}; // Initialize normal to zero vector
    
    int left_vtx = 0;
    int right_vtx = 0;
    int top_vtx = 0;
    int bot_vtx = 0;

    //printf("aa\r\n");
    
    mesh.get_vertex_neighbor_ids(vertex_index, left_vtx, right_vtx, top_vtx, bot_vtx);

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

std::vector<double> cross(std::vector<double> v1, std::vector<double> v2){
  std::vector<double> result(3);
  result[0] = v1[1]*v2[2] - v1[2]*v2[1];
  result[1] = v1[2]*v2[0] - v1[0]*v2[2];
  result[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return result;
}

double dot(std::vector<double> a, std::vector<double> b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

std::vector<double> mult(double a, std::vector<double> b) {
  return {a*b[0], a*b[1], a*b[2]};
}

std::vector<double> add(std::vector<double> a, std::vector<double> b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::vector<double> sub(std::vector<double> a, std::vector<double> b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double magnitude(std::vector<double> a) {
  return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

double cot(const std::vector<double>& a, const std::vector<double>& b) {
    auto cross_product = cross(a, b);
    double cross_magnitude = magnitude(cross_product);

    if (cross_magnitude < 1e-12) {
        //throw std::invalid_argument("Vectors are parallel or one is a zero vector, cotangent undefined.");
        cross_magnitude = 1e-12;
    }

    return dot(a, b) / cross_magnitude;
}

double angle(std::vector<double> v1, std::vector<double> v2)
{
    return std::acos(dot(v1, v2));
}

std::vector<double> compute_laplacian(Mesh &mesh, std::vector<int> adjacent_triangles, std::vector<int> neighboring_vertices, int i) {
    std::vector<double> laplacian(neighboring_vertices.size(), 0.0f);

    for (int j_index = 0; j_index < neighboring_vertices.size(); ++j_index) {
        int j = neighboring_vertices[j_index];

        // Find triangles shared between `i` and `j`
        std::vector<int> shared_triangles;
        for (int triangle : adjacent_triangles) {
          // Check if `j` is one of the vertices in this triangle
          const auto& vertices = mesh.triangles[triangle];
          if (std::find(vertices.begin(), vertices.end(), j) != vertices.end()) {
              shared_triangles.push_back(triangle);
          }
        }

        // Handle cases based on the number of shared triangles
        if (shared_triangles.size() == 2) {
            // Interior case: Two triangles are connected
            std::vector<int> k_vertices;
            for (int triangle : shared_triangles) {
                for (int vertex : mesh.triangles[triangle]) {
                    if (vertex != i && vertex != j) {
                        k_vertices.push_back(vertex);
                        break; // Only one `k` per triangle
                    }
                }
            }

            // Ensure we found two `k` vertices
            if (k_vertices.size() != 2) {
                throw std::runtime_error("Error identifying k vertices in triangles.");
            }

            int k1 = k_vertices[0];
            int k2 = k_vertices[1];

            std::vector<double> edge1;
            std::vector<double> edge2;

            edge1 = sub(mesh.target_points[k1], mesh.target_points[j]);
            edge2 = sub(mesh.target_points[k1], mesh.target_points[i]);
            double cot_k1 = cot(edge1, edge2);

            edge1 = sub(mesh.target_points[k2], mesh.target_points[j]);
            edge2 = sub(mesh.target_points[k2], mesh.target_points[i]);
            double cot_k2 = cot(edge1, edge2);

            laplacian[j_index] = 0.5 * (cot_k1 + cot_k2);

            //std::cout << "k1=" << k1 << ", k2=" << k2 << std::endl;

        } else if (shared_triangles.size() == 1) {
            // Boundary case: Only one triangle is connected
            int triangle = shared_triangles[0];
            int k = -1;

            // Find the single `k` vertex
            for (int vertex : mesh.triangles[triangle]) {
                if (vertex != i && vertex != j) {
                    k = vertex;
                    break;
                }
            }

            if (k == -1) {
                throw std::runtime_error("Error identifying k vertex in boundary triangle.");
            }

            std::vector<double> edge1;
            std::vector<double> edge2;

            edge1 = sub(mesh.target_points[k], mesh.target_points[j]);
            edge2 = sub(mesh.target_points[k], mesh.target_points[i]);
            double cot_k = cot(edge1, edge2);

            laplacian[j_index] = 1;

            //std::cout << "k=" << k << std::endl;

        } else {
            throw std::runtime_error("No shared triangles between i and j; invalid mesh or disconnected vertex.");
        }
    }

    return laplacian;
}

#include <thread>
#include <vector>
#include <mutex>
#include <atomic>
#include <cmath>
#include <iostream>

void poisson_solver_unstructured(
    Mesh &mesh, std::vector<double> &input, std::vector<double> &solution, int num_threads = 1) {
    
    double omega = 1.0;
    const int max_iterations = 1000;
    const double tolerance = 0.01;

    std::mutex max_update_mutex;
    std::atomic<double> max_update;

    std::vector<std::vector<double>> laplacians;

    for (int i = 0; i < mesh.target_points.size(); i++)
    {
        laplacians.push_back(compute_laplacian(mesh, mesh.vertex_adjecent_triangles[i], mesh.vertex_adjecent_vertices[i], i));
    }
    

    auto worker = [&](int start, int end, double &local_max_update) {
    for (int i = start; i < end; ++i) {
        double delta;
        double neighbor_sum = 0.0;
        double neighbor_cnt = 0.0;

        auto neighboring_vertices = mesh.vertex_adjecent_vertices[i];

        std::vector<double> laplacian = laplacians[i];

        for (int j = 0; j < neighboring_vertices.size(); j++)
        {
            double weight = laplacian[j] / magnitude(sub(mesh.target_points[i], mesh.target_points[j]));
            neighbor_cnt += weight;
            neighbor_sum += weight * solution[neighboring_vertices[j]]; // Access via neighbor indices
        }

        // calculate delta
        delta = omega / neighbor_cnt * (neighbor_sum - neighbor_cnt * solution[i] - input[i]);

        //printf("neighbor_cnt = %f, neighbor_sum = %f\r\n", neighbor_cnt, neighbor_sum);

        // get max update
        double abs_delta = fabs(delta);
        if (abs_delta > max_update) {
            max_update = abs_delta;
        }

        // increment the value by delta
        solution[i] += delta;
    }
    };

    for (int itr = 0; itr < max_iterations; ++itr) {
        max_update = 0.0;

        // Divide the workload among threads
        std::vector<std::thread> threads;
        std::vector<double> local_max_updates(num_threads, 0.0);

        int chunk_size = (mesh.target_points.size() + num_threads - 1) / num_threads;

        for (int t = 0; t < num_threads; ++t) {
            int start = t * chunk_size;
            int end = std::min(start + chunk_size, (int)mesh.target_points.size());
            if (start >= end) break; // No work for this thread
            threads.emplace_back(worker, start, end, std::ref(local_max_updates[t]));
        }

        // Join all threads
        for (auto &thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        // Calculate global max_update
        for (const auto &local_update : local_max_updates) {
            if (local_update > max_update) {
                max_update = local_update;
            }
        }

        //std::cout << "Iteration " << itr << ", max_update = " << max_update << "\n";

        // Check for convergence
        if (max_update < tolerance) {
            break;
        }
    }
}


std::vector<std::vector<double>> calculate_mesh_gradient(Mesh &mesh, std::vector<double> &input) {
    std::vector<std::vector<double>> triangle_gradients(mesh.triangles.size());
    std::vector<double> triangle_areas(mesh.triangles.size());
    std::vector<std::vector<double>> vertex_gradients(mesh.target_points.size());

    for (int i = 0; i < mesh.triangles.size(); i++) {
        std::vector<int> triangle = mesh.triangles[i];
        std::vector<double> vertex_i = mesh.target_points[triangle[0]];
        std::vector<double> vertex_j = mesh.target_points[triangle[1]];
        std::vector<double> vertex_k = mesh.target_points[triangle[2]];

        // get triangle area
        std::vector<std::vector<double>> polygon = {vertex_i, vertex_j, vertex_k};
        double area = calculate_polygon_area_vec(polygon);

        // get edges
        std::vector<double> edge_ik = sub(vertex_i, vertex_k);
        std::vector<double> edge_ji = sub(vertex_j, vertex_i);

        // rotate by 90°
        std::vector<double> edge_ik_t = {-edge_ik[1], edge_ik[0], edge_ik[2]};
        std::vector<double> edge_ji_t = {-edge_ji[1], edge_ji[0], edge_ji[2]};

        // function gradient on edges
        double f_ji = input[triangle[1]] - input[triangle[0]];
        double f_ki = input[triangle[2]] - input[triangle[0]];

        // gradient
        std::vector<double> a = mult(f_ji / (2.0*area), edge_ik_t);
        std::vector<double> b = mult(f_ki / (2.0*area), edge_ji_t);
        triangle_gradients[i] = add(a, b);
        triangle_areas[i] = area;
    }

    for (int i = 0; i < mesh.target_points.size(); i++) {
        auto [adjacent_edges, adjacent_triangles, neighboring_vertices] = mesh.find_adjacent_elements(i);

        std::vector<double> vertex_gradient(3, 0.0);
        double total_weight = 0.0f;

        for (int j = 0; j < adjacent_triangles.size(); j++) {
            std::vector<double> edge_ij = sub(mesh.target_points[mesh.triangles[adjacent_triangles[j]][0]], mesh.target_points[mesh.triangles[adjacent_triangles[j]][1]]);
            std::vector<double> edge_ik = sub(mesh.target_points[mesh.triangles[adjacent_triangles[j]][0]], mesh.target_points[mesh.triangles[adjacent_triangles[j]][2]]);
            //double weight = 1.0;// / angle(edge_ij, edge_ik);
            double weight = triangle_areas[j];
            vertex_gradient = add(vertex_gradient, mult(weight, triangle_gradients[adjacent_triangles[j]]));
            total_weight += weight;
        }

        vertex_gradients[i] = mult(1.0 / total_weight, vertex_gradient);
    }
    
    return vertex_gradients;
}

std::vector<double> compute_divergence(Mesh &mesh, const std::vector<std::vector<double>> &vector_field) {
    std::vector<double> divergence(mesh.target_points.size(), 0.0);

    for (int i = 0; i < mesh.triangles.size(); i++) {
        std::vector<int> triangle = mesh.triangles[i];
        std::vector<double> vertex_i = mesh.target_points[triangle[0]];
        std::vector<double> vertex_j = mesh.target_points[triangle[1]];
        std::vector<double> vertex_k = mesh.target_points[triangle[2]];

        // Compute area
        std::vector<std::vector<double>> polygon = {vertex_i, vertex_j, vertex_k};
        double area = calculate_polygon_area_vec(polygon);

        // Compute shape function gradients (already defined in your gradient function)
        std::vector<double> edge_ik = sub(vertex_i, vertex_k);
        std::vector<double> edge_ji = sub(vertex_j, vertex_i);

        // Rotated by 90 degrees
        std::vector<double> edge_ik_t = {-edge_ik[1], edge_ik[0], edge_ik[2]};
        std::vector<double> edge_ji_t = {-edge_ji[1], edge_ji[0], edge_ji[2]};

        // Compute per-triangle divergence contribution
        for (int j = 0; j < 3; j++) {
            int v_id = triangle[j];

            // Shape function gradient for vertex v_id
            std::vector<double> grad_phi;
            if (j == 0) {
                grad_phi = add(mult(1.0 / (2.0 * area), edge_ik_t), mult(1.0 / (2.0 * area), edge_ji_t));
            } else if (j == 1) {
                grad_phi = mult(-1.0 / (2.0 * area), edge_ik_t);
            } else {
                grad_phi = mult(-1.0 / (2.0 * area), edge_ji_t);
            }

            // Compute divergence at v_id
            divergence[v_id] += dot_product(grad_phi, vector_field[i]) * area;
        }
    }

    return divergence;
}

void subtractAverageVec(std::vector<double>& values) {
    // Calculate the average
    double sum = 0.0;
    for (size_t i = 0; i < values.size(); ++i) {
        sum += values[i];
    }
    double average = sum / values.size();

    // Subtract the average from each element
    for (size_t i = 0; i < values.size(); ++i) {
        values[i] -= average;
    }
}

double Caustic_design::perform_transport_iteration() {
    //std::vector<std::vector<double>> vertex_gradient;
    double min_step;

    // build median dual mesh of the updated parameterization
    target_cells.clear();
        mesh.build_target_dual_cells(target_cells);

    // calculate difference D (interpretation of equation 2)
    std::vector<double> source_areas = get_source_areas(target_cells);
    calculate_errors(source_areas, target_areas, target_cells, errors);

    // rasterize the mesh into a uniform rectangular matrix
    bool triangle_miss = false;
    raster = mesh.interpolate_raster_target(errors, resolution_x, resolution_y, triangle_miss);
    
    if (triangle_miss) {
        mesh.laplacian_smoothing(mesh.target_points, 0.1f);
        return NAN;
    }

    // solve the poisson equation 3 in the paper
    subtractAverage(raster);
    poisson_solver(raster, phi, resolution_x, resolution_y, 100000, 0.0000001, nthreads);

    // calculate the gradient given by equation 4
    gradient = calculate_gradient(phi);

    // calculate the gradient vectors corresponding to each vertex in the mesh

    // bilinear interpolating the gradients (negligibly faster, but gives lower contrast results)
    double epsilon = 1e-8;
    

    std::vector<double> vertex_gradient_x;
    std::vector<double> vertex_gradient_y;
    for (int i=0; i<mesh.target_points.size(); i++) {
        vertex_gradient_x.push_back(bilinearInterpolation(gradient[0], 
            (mesh.target_points[i][0] / mesh.width) * (resolution_x) - 0.5, 
            (mesh.target_points[i][1] / mesh.height) * (resolution_y) - 0.5
        ));

        vertex_gradient_y.push_back(bilinearInterpolation(gradient[1], 
            (mesh.target_points[i][0] / mesh.width) * (resolution_x) - 0.5, 
            (mesh.target_points[i][1] / mesh.height) * (resolution_y) - 0.5
        ));
    }//*/

    export_cells_as_svg(target_cells, scale_array_proportional(vertex_gradient_x, 0.0f, 1.0f), "../vertex_gradient_x.svg");
    export_cells_as_svg(target_cells, scale_array_proportional(vertex_gradient_y, 0.0f, 1.0f), "../vertex_gradient_y.svg");

    vertex_gradient.clear();
    vertex_gradient.push_back(vertex_gradient_x);
    vertex_gradient.push_back(vertex_gradient_y);
    
    //*/

    // integrate the gradient grid into the dual cells of the vertices (slower but better contrast)
    //vertex_gradient = integrate_cell_gradients(gradient, target_cells, resolution_x, resolution_y, width, height);

    std::vector<std::vector<double>> old_points;

    std::copy(mesh.target_points.begin(), mesh.target_points.end(), back_inserter(old_points));

    // step the mesh vertices in the direction of their gradient vector
    //mesh.step_grid(vertex_gradient[0], vertex_gradient[1], 0.0005f);
    mesh.step_grid(vertex_gradient[0], vertex_gradient[1], 0.05f);

    //mesh.laplacian_smoothing(mesh.target_points, 0.5f);

    min_step = 0.0f;

    for (int i=0; i<old_points.size(); i++) {
        double dx = (old_points[i][0] - mesh.target_points[i][0]);
        double dy = (old_points[i][1] - mesh.target_points[i][1]);
        double dz = (old_points[i][2] - mesh.target_points[i][2]);

        double dist = sqrt(dx*dx + dy*dy + dz*dz);

        if (min_step < dist) {
            min_step = dist;
        }
    }

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
    normals = mesh.calculate_refractive_normals_uniform(resolution_x / width * focal_l, 1.49);

    std::vector<std::vector<double>> vector_field;

    for (int i = 0; i < normals[0].size(); i++)
    {
        vector_field.push_back({
            normals[0][i],
            normals[1][i],
            1.0f
        });

        vector_field[i] = normalize(vector_field[i]);
    }
    
    std::vector<double> vertex_divergance = compute_divergence(mesh, vector_field);
    //export_cells_as_svg(target_cells, scale_array_proportional(vertex_divergance, 0.0f, 1.0f), "../vertex_divergance.svg");

    // interpolates the vertex normals into a large uniform grid
    mesh.build_source_bvh(5, 30);
    bool triangle_miss = false;
    norm_x = mesh.interpolate_raster_source(normals[0], resolution_x, resolution_y, triangle_miss);
    norm_y = mesh.interpolate_raster_source(normals[1], resolution_x, resolution_y, triangle_miss);

    if (triangle_miss) {
        return;
    }

    // calculates the divergance of the interpolated normals
    divergence = calculate_divergence(norm_x, norm_y, resolution_x, resolution_y);
    subtractAverage(divergence);

    // solve the poisson equation for the divergance
    poisson_solver(divergence, h, resolution_x, resolution_y, 100000, 0.00000001, nthreads);

    /*std::vector<double> interpolated_h;
    for (int i=0; i<mesh.target_points.size(); i++) {
        interpolated_h.push_back(bilinearInterpolation(h, mesh.target_points[i][0] * ((resolution_x) / mesh.width), mesh.target_points[i][1] * ((resolution_y) / mesh.height)));
    }
    double max_update = mesh.set_target_heights(interpolated_h);
    printf("height max update %.5e\r\n", max_update);*/

    double epsilon = 1e-8;

    // get the heights on the vertex positions
    std::vector<double> interpolated_h;
    for (int i=0; i<mesh.source_points.size(); i++) {
        interpolated_h.push_back(bilinearInterpolation(h, 
            (mesh.source_points[i][0] / mesh.width) * (resolution_x) - 0.5, 
            (mesh.source_points[i][1] / mesh.height) * (resolution_y) - 0.5
        ));
    }
    double max_update = mesh.set_source_heights(interpolated_h);
    printf("height max update %.5e\r\n", max_update);
}

void Caustic_design::initialize_solvers(std::vector<std::vector<double>> image) {
    pixels = scale_matrix_proportional(image, 0, 1.0f);

    printf("scaled\r\n");

    //Mesh mesh(width, height, mesh_res_x, mesh_res_y);

    // TODO put mesh into local private object

    printf("generated mesh\r\n");

    mesh.export_to_svg("../mesh.svg", 1);

    //std::cout << "built mesh" << std::endl;

    //std::vector<std::vector<std::vector<double>>> circ_target_cells;
    mesh.build_target_dual_cells(target_cells);
    mesh.build_source_dual_cells(source_cells);
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
