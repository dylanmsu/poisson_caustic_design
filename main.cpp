#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "src/caustic_design.h"

#include <thread>

#define cimg_use_png
//#define cimg_use_jpeg

#include "cimg/CImg.h"
#include "lsqcpp/lsqcpp.hpp"

#include "ceres/ceres.h"
#include "costFunctor.h"

#include "glm/glm.hpp"

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

Caustic_design caustic_design;

vector<glm::vec3> x_sources;

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

template <typename T>
std::vector<T> vector_subtract_ceres(const std::vector<T>& v1, const std::vector<T>& v2) {
    return { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
}

template<typename T> 
std::vector<T> cross_ceres(std::vector<T> v1, std::vector<T> v2) {
    std::vector<T> result(3);
    result[0] = T(v1[1]*v2[2] - v1[2]*v2[1]);
    result[1] = T(v1[2]*v2[0] - v1[0]*v2[2]);
    result[2] = T(v1[0]*v2[1] - v1[1]*v2[0]);
    return result;
}

template <typename T>
std::vector<T> normalize_ceres(const std::vector<T>& v) {
    T norm = ceres::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return { v[0] / norm, v[1] / norm, v[2] / norm };
}


// Compute the angle between two vectors
template <typename T>
T compute_angle(const std::vector<T> &v1, const std::vector<T> &v2) {
    T dot_product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    T length_v1 = ceres::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    T length_v2 = ceres::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    return ceres::acos(dot_product / (length_v1 * length_v2));
}

// Find the normal vector for the given vertex
template <typename T>
std::vector<T> calculate_vertex_normal_ceres(Caustic_design* cd, std::vector<std::vector<T>> &points, int vertex_index) {
    std::vector<T> final_normal = {T(0.0), T(0.0), T(0.0)};
    T total_weight = T(0.0);

    // Find triangles connected to the vertex
    auto triangles_containing_vertex = cd->mesh->vertex_to_triangles.find(vertex_index);
    if (triangles_containing_vertex != cd->mesh->vertex_to_triangles.end()) {
        // Iterate over each triangle containing the vertex
        for (int triangle_index : triangles_containing_vertex->second) {
            std::vector<int> triangle = cd->mesh->triangles[triangle_index];

            // Get the three vertices of the triangle
            std::vector<T> p0 = points[triangle[0]];
            std::vector<T> p1 = points[triangle[1]];
            std::vector<T> p2 = points[triangle[2]];

            // Determine which vertex in the triangle corresponds to the vertex_index
            std::vector<T> vertex, v1, v2;
            if (triangle[0] == vertex_index) {
                vertex = p0; v1 = p1; v2 = p2;
            } else if (triangle[1] == vertex_index) {
                vertex = p1; v1 = p2; v2 = p0;
            } else {
                vertex = p2; v1 = p0; v2 = p1;
            }

            // Compute two edges of the triangle
            std::vector<T> edge1 = {v1[0] - vertex[0], v1[1] - vertex[1], v1[2] - vertex[2]};
            std::vector<T> edge2 = {v2[0] - vertex[0], v2[1] - vertex[1], v2[2] - vertex[2]};

            // Compute the normal of the triangle (cross product of two edges)
            std::vector<T> normal = cross_ceres(edge1, edge2);
            normal = normalize_ceres(normal);

            // Compute the angle at the given vertex (between edge1 and edge2)
            T angle = compute_angle(edge1, edge2);

            // Add the weighted normal to the final normal
            final_normal[0] += normal[0] * angle;
            final_normal[1] += normal[1] * angle;
            final_normal[2] += normal[2] * angle;

            // Accumulate total weight (angle)
            total_weight += angle;
        }

        // Normalize the final normal (weighted by angles)
        final_normal[0] /= total_weight;
        final_normal[1] /= total_weight;
        final_normal[2] /= total_weight;

        // Return the normalized final normal
        return normalize_ceres(final_normal);
    }

    // If no triangles are connected, return a zero vector
    return {T(0.0), T(0.0), T(0.0)};
}

struct ParabolicError {
    int vertex_index;
    Caustic_design* cd;  // Assuming caustic_design.mesh is a pointer to the mesh.

    // Constructor to initialize the vertex index and mesh pointer.
    ParabolicError(int index, Caustic_design* cd) : vertex_index(index), cd(cd) { }

    template<typename T>
    bool operator()(const T* const zval, T* residual) const {
        T E_int = T(0.0);

        std::vector<std::vector<T>> temp_source_points(cd->mesh->source_points.size());

        // Copy the mesh source points from double to T
        for (size_t i = 0; i < cd->mesh->source_points.size(); ++i) {
            temp_source_points[i].resize(3);  // Assume 3D points
            temp_source_points[i][0] = T(cd->mesh->source_points[i][0]);  // Copy x-coordinate
            temp_source_points[i][1] = T(cd->mesh->source_points[i][1]);  // Copy y-coordinate
            temp_source_points[i][2] = T(cd->mesh->source_points[i][2]);  // Copy z-coordinate
        }

        // Set the z-value of the source point for the vertex
        temp_source_points[vertex_index][2] = zval[0];

            // Calculate the target normal at the neighboring vertex
            std::vector<T> normal_trg = {
                T(-cd->normals[0][vertex_index]),
                T(-cd->normals[1][vertex_index]),
                T(cd->normals[2][vertex_index])
            };

            // Calculate the current normal at the neighboring vertex
            std::vector<T> normal = calculate_vertex_normal_ceres(cd, temp_source_points, vertex_index);

            // Normalize both normal vectors
            normal_trg = normalize_ceres(normal_trg);
            normal = normalize_ceres(normal);

            // Compute the difference between the target and current normals
            std::vector<T> diff = vector_subtract_ceres<T>(normal, normal_trg);

            // Calculate energy contribution for this neighbor (squared difference)
            T energy = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
            E_int += energy;

        // Set the final parabolic error value (residual)
        residual[0] = E_int;
        return true;
    }
};

/**
 * @brief TargetOptimization::gatherVertexInformation Gathers information about neighbors of current vertex.
 * @param vertex The current vertex
 * @param vertexIndex The index of the current vertex
 * @param neighborList [out] A list of indices to the neighbors of the current vertex. Each neighbor is only included once
 * @param neighborMap [out] A list of indices to indices to neighbors. The value is an index to the neighborList. Is used to handle several references to one neighbor (e.g. for neighboring faces)
 * @param gridNeighbors [out] A list indices to of horizontal and vertical neighbors
 */
void gatherVertexInformation(Mesh &mesh, uint vertexIndex, vector<int> &neighborList, vector<int> &neighborMap, vector<int> & gridNeighbors)
{
    mesh.find_vertex_connectivity(vertexIndex, neighborList, neighborMap);

    int left_vtx = 0;
    int right_vtx = 0;
    int top_vtx = 0;
    int bot_vtx = 0;
    
    mesh.get_vertex_neighbor_ids(vertexIndex, left_vtx, right_vtx, top_vtx, bot_vtx);

    if (left_vtx != -1) {
        gridNeighbors.push_back(left_vtx);
    }

    if (right_vtx != -1) {
        gridNeighbors.push_back(right_vtx);
    }

    if (top_vtx != -1) {
        gridNeighbors.push_back(top_vtx);
    }

    if (bot_vtx != -1) {
        gridNeighbors.push_back(bot_vtx);
    }
}

/**
 * @brief TargetOptimization::addResidualBlocks Adds residuals blocks for the given vertex to the given problem.
 * @param problem The problem to be solved
 * @param vertexIndex The vertex index (in the edgeVertices-vector of the mesh)
 * @param neighbors The neighboring-vertices of the current vertex. These include only the ones that share a face with the current vertex
 * @param neighborMap The index of the neighbors. Two neighbors make one face (with the vertex itself). So vertex + neighbors[neighborMap[0]] + neighbors[neighborMap[1]] forms one face
 * @param gridNeighbors The horizontal and vertical neighbors of the current vertex
 * @param vertices The reference to the vertices. Each vertex uses three doubles (x,y,z).
 */
void addResidualBlocks(Problem *problem, uint vertexIndex, vector<int> &neighbors, vector<int> &neighborMap, vector<int> &gridNeighbors, double *vertices, glm::vec3 &trg_normal)
{
    float weightMult = 1.0;
    if(caustic_design.mesh->is_border(vertexIndex)) // we have an edge. Set weight for edir extremely high
        weightMult = 10000;

    // EDir depends on the original position
    CostFunction* cost_function_edir =
            new AutoDiffCostFunction<CostFunctorEdir2, 3, 3>(new CostFunctorEdir2(&x_sources[vertexIndex], weightMult));

    problem->AddResidualBlock(
                cost_function_edir,
                NULL, // no loss function
                &vertices[vertexIndex*3]
                );//*/


    if(true){//caustic_design.mesh->is_border(vertexIndex)){ //not an edge we optimize the normals
        // For eint we have several functors, each for a different amount of neighbors
        switch(neighbors.size()){
        case 2:
        {

            /*CostFunction* cost_function_eint2 =
                new AutoDiffCostFunction<CostFunctorEint2Neighbors, 3, 3, 3, 3>(new CostFunctorEint2Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint2, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3]);//*/


            break;
            }
        case 3:
        {

            /*CostFunction* cost_function_eint3 =
                new AutoDiffCostFunction<CostFunctorEint3Neighbors, 3, 3, 3, 3, 3>(new CostFunctorEint3Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint3, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3],
                                       &vertices[neighbors[2]*3]);//*/


            break;
            }

        case 4:
        {

            /*CostFunction* cost_function_eint4 =
                new AutoDiffCostFunction<CostFunctorEint4Neighbors, 3, 3, 3, 3, 3, 3>(new CostFunctorEint4Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint4, NULL,
                                       &vertices[vertexIndex*3], // vertex
                                       &vertices[neighbors[0]*3], // and the neighbors..
                                       &vertices[neighbors[1]*3],
                                       &vertices[neighbors[2]*3],
                                       &vertices[neighbors[3]*3]);//*/


            break;
            }

        case 5:
        {
            /*CostFunction* cost_function_eint5 =
                new AutoDiffCostFunction<CostFunctorEint5Neighbors, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint5Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock(cost_function_eint5, NULL,
                           &vertices[vertexIndex*3], // vertex
                           &vertices[neighbors[0]*3], // and the neighbors..
                           &vertices[neighbors[1]*3],
                           &vertices[neighbors[2]*3],
                           &vertices[neighbors[3]*3],
                           &vertices[neighbors[4]*3]);//*/



            break;
        }

        case 6:
        {
            CostFunction* cost_function_eint6 =
                new AutoDiffCostFunction<CostFunctorEint6Neighbors, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint6Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock( cost_function_eint6, NULL,
                       &vertices[vertexIndex*3], // vertex
                       &vertices[neighbors[0]*3], // and the neighbors..
                       &vertices[neighbors[1]*3],
                       &vertices[neighbors[2]*3],
                       &vertices[neighbors[3]*3],
                       &vertices[neighbors[4]*3],
                       &vertices[neighbors[5]*3]);//*/



            break;
        }

        case 7:
        {
            /*CostFunction* cost_function_eint7 =
                new AutoDiffCostFunction<CostFunctorEint7Neighbors, 3, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint7Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock( cost_function_eint7, NULL,
                       &vertices[vertexIndex*3], // vertex
                       &vertices[neighbors[0]*3], // and the neighbors..
                       &vertices[neighbors[1]*3],
                       &vertices[neighbors[2]*3],
                       &vertices[neighbors[3]*3],
                       &vertices[neighbors[4]*3],
                       &vertices[neighbors[5]*3],
                       &vertices[neighbors[6]*3]);//*/


            break;
        }

        case 8:
        {
            /*CostFunction* cost_function_eint8 =
                new AutoDiffCostFunction<CostFunctorEint8Neighbors, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3>(new CostFunctorEint8Neighbors(trg_normal, neighborMap));

            problem->AddResidualBlock( cost_function_eint8, NULL,
                           &vertices[vertexIndex*3], // vertex
                           &vertices[neighbors[0]*3], // and the neighbors..
                           &vertices[neighbors[1]*3],
                           &vertices[neighbors[2]*3],
                           &vertices[neighbors[3]*3],
                           &vertices[neighbors[4]*3],
                           &vertices[neighbors[5]*3],
                           &vertices[neighbors[6]*3],
                           &vertices[neighbors[7]*3]);//*/

            break;

        }

        } // switch end
    }//*/

    // regularization term with 4-neighborhood
    CostFunction* ereg;
    switch(gridNeighbors.size())
    {
    case 2:

        ereg = new AutoDiffCostFunction<CostFunctorEreg2Neighbors, 3, 3, 3, 3>(new CostFunctorEreg2Neighbors(gridNeighbors));
        problem->AddResidualBlock(
                    ereg,
                    NULL,
                    &vertices[vertexIndex*3],
                    &vertices[gridNeighbors[0]*3],
                    &vertices[gridNeighbors[1]*3]
                    );
        break;
    case 3:

        ereg = new AutoDiffCostFunction<CostFunctorEreg3Neighbors, 3, 3, 3, 3, 3>(new CostFunctorEreg3Neighbors(gridNeighbors));
        problem->AddResidualBlock(
                    ereg,
                    NULL,
                    &vertices[vertexIndex*3],
                    &vertices[gridNeighbors[0]*3],
                    &vertices[gridNeighbors[1]*3],
                    &vertices[gridNeighbors[2]*3]
                    );
        break;
    case 4:

        ereg = new AutoDiffCostFunction<CostFunctorEreg4Neighbors, 3, 3, 3, 3, 3, 3>(new CostFunctorEreg4Neighbors(gridNeighbors));
        problem->AddResidualBlock(
                    ereg,
                    NULL,
                    &vertices[vertexIndex*3],
                    &vertices[gridNeighbors[0]*3],
                    &vertices[gridNeighbors[1]*3],
                    &vertices[gridNeighbors[2]*3],
                    &vertices[gridNeighbors[3]*3]
                    );
        break;
    }//*/
}


int main(int argc, char const *argv[])
{
    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    // Load image
    std::string image_path = args["input_png"];
    cimg_library::CImg<unsigned char> image(image_path.c_str());
    double aspect_ratio = (double)image.width() / (double)image.height();

    // Print parsed arguments
    /*for (const auto& pair : args) {
        //std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
        printf("Key: %s, Value: %s\r\n", pair.first.c_str(), pair.second.c_str());
    }*/

    int mesh_resolution_x = atoi(args["res_w"].c_str());
    double mesh_width = std::stod(args["width"]);

    caustic_design.set_mesh_resolution(mesh_resolution_x, mesh_resolution_x / aspect_ratio);
    caustic_design.set_domain_resolution(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio);

    double mesh_height = floor((mesh_resolution_x) / aspect_ratio) * (mesh_width / (mesh_resolution_x));

    caustic_design.set_mesh_size(mesh_width, mesh_height);

    caustic_design.set_lens_focal_length(std::stod(args["focal_l"]));
    caustic_design.set_lens_thickness(std::stod(args["thickness"]));
    caustic_design.set_solver_max_threads(atoi(args["max_threads"].c_str()));

    image = image.resize(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, -100, -100, 3); // Resize using linear interpolation

    // Convert image to grid
    std::vector<std::vector<double>> pixels;
    
    image_to_grid(image, pixels);

    caustic_design.initialize_solvers(pixels);

    printf("a\r\n");

    //caustic_design.export_paramererization_to_svg("../parameterization_0.svg", 0.5f);

    printf("b\r\n");
    
    for (int itr=0; itr<6; itr++) {
        printf("starting iteration %i\r\n", itr);
        
        double step_size = caustic_design.perform_transport_iteration();

        //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(caustic_design.vertex_gradient[0], 0.0f, 1.0f), "../x_grad.svg");

        //save_grid_as_image(scale_matrix_proportional(caustic_design.gradient[0], 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../grad_x_" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.gradient[1], 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../grad_y_" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.raster, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../raster_" + std::to_string(itr) + ".png");

        caustic_design.export_paramererization_to_svg("../parameterization_" + std::to_string(itr + 1) + ".svg", 1.0f);
        caustic_design.export_inverted_transport_map("../inverted.svg", 1.0f);

        printf("step_size = %f\r\n", step_size);

        // convergence treshold for the parameterization
        if (step_size < 0.01) break;
    }

    printf("\033[0;32mTransport map solver done! Starting height solver.\033[0m\r\n");

    for (int itr=0; itr<3; itr++) {
        //caustic_design.perform_height_map_iteration(itr);
        //save_grid_as_image(scale_matrix_proportional(caustic_design.h, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../h" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.divergence, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../div" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.norm_x, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "norm_x" + std::to_string(itr) + ".png");
        //save_grid_as_image(scale_matrix_proportional(caustic_design.norm_y, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "norm_y" + std::to_string(itr) + ".png");
    }

    caustic_design.normals.clear();
    caustic_design.normals = caustic_design.mesh->calculate_refractive_normals_uniform(caustic_design.focal_l, 1.49);

    std::vector<glm::vec3> desired_normals;

    // make a copy of the original positions of the vertices
    for (int i = 0; i < caustic_design.mesh->source_points.size(); i++) {
        glm::vec3 v;
        v.x = caustic_design.mesh->source_points[i][0];
        v.y = caustic_design.mesh->source_points[i][1];
        v.z = caustic_design.mesh->source_points[i][2];
        x_sources.push_back(v);

        glm::vec3 trg_normal = {caustic_design.normals[0][i], caustic_design.normals[1][i], caustic_design.normals[2][i]};

        desired_normals.push_back(trg_normal);
    }

    vector<vector<int> > neighborsPerVertex;
    neighborsPerVertex.resize(caustic_design.mesh->source_points.size());

    vector<vector<int> > neighborMapPerVertex;
    neighborMapPerVertex.resize(caustic_design.mesh->source_points.size());

    vector<vector<int> > eightNeighborsPerVertex;
    eightNeighborsPerVertex.resize(caustic_design.mesh->source_points.size());

    // gather information for each vertex to optimize
    for(uint i = 0; i < caustic_design.mesh->source_points.size(); i++)
    {
        vector<int> neighbors;
        vector<int> neighborMap;
        vector<int> eightNeighbors;

        gatherVertexInformation(*caustic_design.mesh, i, neighbors, neighborMap, eightNeighbors);

        neighborsPerVertex[i] = neighbors;
        neighborMapPerVertex[i] = neighborMap;
        eightNeighborsPerVertex[i] = eightNeighbors;
    }

    // put all positions in one big list that we access later
    double* vertices = new double[3*caustic_design.mesh->source_points.size()];
    for(uint i=0; i<caustic_design.mesh->source_points.size(); i++)
    {
        vertices[3*i + 0] = caustic_design.mesh->source_points[i][0];
        vertices[3*i + 1] = caustic_design.mesh->source_points[i][1];
        vertices[3*i + 2] = caustic_design.mesh->source_points[i][2];
    }

    Problem prob;

    // iterate over all vertices and add the corresponding residual blocks
    for(uint i=0; i<caustic_design.mesh->source_points.size(); i++)
    {
        /*std::cout << "block " << i << std::endl;
        std::cout << "neighborsPerVertex " << neighborsPerVertex[i].size() <<  std::endl;
        std::cout << "neighborMapPerVertex " << neighborMapPerVertex[i].size() <<  std::endl;
        std::cout << "eightNeighborsPerVertex " << eightNeighborsPerVertex[i].size() <<  std::endl;*/
        addResidualBlocks(&prob, i, neighborsPerVertex[i], neighborMapPerVertex[i], eightNeighborsPerVertex[i], vertices, desired_normals[i]);
    }

    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::ITERATIVE_SCHUR;
    options.max_num_iterations = 200;
    options.dense_linear_algebra_library_type = ceres::LAPACK;
    options.num_threads = 16;

    string error;
    if(!options.IsValid(&error))
    {
        std::cout << "Options not valid: " << error << std::endl;
    }

    Solver::Summary summary;
    Solve(options, &prob, &summary);

    std::cout << summary.FullReport() << std::endl;

    glm::vec3 * pos;
    for(uint i=0; i<caustic_design.mesh->source_points.size(); i++)
    {
        caustic_design.mesh->source_points[i][0] = vertices[3*i + 0];
        caustic_design.mesh->source_points[i][1] = vertices[3*i + 1];
        caustic_design.mesh->source_points[i][2] = vertices[3*i + 2];
    }

    /*// Create a Ceres Problem
    ceres::Problem problem;

    // Loop through all vertices of the mesh and add residual blocks
    for (int i = 0; i < caustic_design.mesh->source_points.size(); ++i) {

        // Create the cost function for the vertex `i`
        ceres::CostFunction* cost_function =
            new ceres::AutoDiffCostFunction<ParabolicError, 1, 1>(new ParabolicError(i, &caustic_design));

        // Add a residual block for this vertex to the problem
        problem.AddResidualBlock(cost_function, nullptr, &caustic_design.mesh->source_points[i][2]);
    }

    // Set up solver options
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.linear_solver_type = ceres::ITERATIVE_SCHUR;
    options.max_num_iterations = 200;
    options.dense_linear_algebra_library_type = ceres::LAPACK;
    options.num_threads = 16;

    // Run the solver
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // Output solver summary
    std::cout << summary.FullReport() << "\n";*/

    //printf("Done! Converged: %s Iterations: %d\n", result.converged ? "true" : "false", result.iterations);

    // do something with final function value
    //printf("Final fval: %s\n", result.fval.transpose());

    /*bool miss = false;
    std::vector<std::vector<double>> norm_x = caustic_design.mesh->interpolate_raster_source(x_normals, 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, miss);
    std::vector<std::vector<double>> norm_y = caustic_design.mesh->interpolate_raster_source(y_normals, 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, miss);
    save_grid_as_image(scale_matrix_proportional(norm_x, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../x_normals.png");
    save_grid_as_image(scale_matrix_proportional(norm_y, 0.0f, 1.0f), 4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, "../y_normals.png");
    //*/

    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(E_int, 0.0f, 1.0f), "../integration_energy.svg");
    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(y_normals_trg, 0.0f, 1.0f), "../y_normals_trg.svg");
    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(x_normals, 0.0f, 1.0f), "../x_normals.svg");
    //export_cells_as_svg(caustic_design.source_cells, scale_array_proportional(y_normals, 0.0f, 1.0f), "../y_normals.svg");

    printf("Height solver done! Exporting as solidified obj\r\n");

    //caustic_design.save_solid_obj_target("../output.obj");
    caustic_design.save_solid_obj_source("../output.obj");

    return 0;
}
