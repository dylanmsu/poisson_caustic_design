#include "mesh.h"

Mesh::Mesh(double width, double height, int res_x, int res_y)
{
    // set physical size of mesh
    this->width = width;
    this->height = height;

    // set poisson domain resolution
    this->res_x = res_x;
    this->res_y = res_y;

    // Build the parameterization mesh
    generate_structured_mesh(res_x, res_y, width, height, this->triangles, this->target_points);
    build_vertex_to_triangles();

    //circular_transform(this->target_points);

    // Duplicate mesh points
    for (int i=0; i<this->target_points.size(); i++) {
        this->source_points.push_back(this->target_points[i]);
    }

    // Create instance of the bvh class used for interpolation
    target_bvh = new Bvh(triangles, target_points);
    source_bvh = new Bvh(triangles, source_points);
}

Mesh::~Mesh()
{
    delete(target_bvh);
    delete(source_bvh);
}

// Build the BVH tree for the target mesh
void Mesh::build_target_bvh(int targetCellSize, int maxDepth) {
    target_bvh->build(targetCellSize, maxDepth);
}

void Mesh::build_source_bvh(int targetCellSize, int maxDepth) {
    source_bvh->build(targetCellSize, maxDepth);
}

// generates a structured triangulation used for the parameterization
void Mesh::generate_structured_mesh(int nx, int ny, double width, double height, std::vector<std::vector<int>> &triangles, std::vector<point_t> &points) {
    printf("%i, %i, %f, %f\r\n", nx, ny, width, height);
    // Generate points
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            double x = static_cast<double>(j) * width / (nx - 1);
            double y = static_cast<double>(i) * height / (ny - 1);
            points.push_back({x, y, 0.0});
        }
    }

    // Generate triangles
    for (int i = 0; i < ny - 1; ++i) {
        for (int j = 0; j < nx - 1; ++j) {
            int idx = i * nx + j;
            triangles.push_back({idx, idx + 1, idx + nx});
            triangles.push_back({idx + nx, idx + 1, idx + nx + 1});
        }
    }
}

// transforms a square grid into a circular grid -> to support circular lenses in the future
/*void Mesh::circular_transform(std::vector<point_t> &input_points) {
    for (int i = 0; i < input_points.size(); i++) {
        double x = input_points[i][0] - this->width/2.0f;
        double y = input_points[i][1] - this->height/2.0f;

        input_points[i][0] = x * sqrt(1.0 - 2.0*(y * y));
        input_points[i][1] = y * sqrt(1.0 - 2.0*(x * x));

        input_points[i][0] += this->width/2.0f;
        input_points[i][1] += this->height/2.0f;
    }
}*/

std::vector<point_t> Mesh::circular_transform(std::vector<point_t> &input_points) {
    std::vector<point_t> transformed_points;
    for (int i = 0; i < input_points.size(); i++) {
        point_t transformed_point(3);

        double x = input_points[i][0] - this->width/2.0f;
        double y = input_points[i][1] - this->height/2.0f;

        transformed_point[0] = x * sqrt(1.0 - 2.0*(y * y));
        transformed_point[1] = y * sqrt(1.0 - 2.0*(x * x));

        transformed_point[0] += this->width/2.0f;
        transformed_point[1] += this->height/2.0f;
        transformed_point[2] = input_points[i][2];

        transformed_points.push_back(transformed_point);
    }
    return transformed_points;
}

// export triangular mesh (target) as svg
void Mesh::export_to_svg(std::string filename, double stroke_width) {
    export_triangles_to_svg(this->target_points, this->triangles, this->width, this->height, this->res_x, this->res_y, filename, stroke_width);
}

// build mapping from vertices to adjecent triangles -> used for creating dual cells
void Mesh::build_vertex_to_triangles() {
    for (int i = 0; i < this->triangles.size(); ++i) {
        const std::vector<int>& triangle = this->triangles[i];
        
        for (int vertex : triangle) {
            if (this->vertex_to_triangles.find(vertex) == this->vertex_to_triangles.end()) {
                this->vertex_to_triangles[vertex] = std::vector<int>();
            }

            this->vertex_to_triangles[vertex].push_back(i);
        }
    }
}

// find triangles and edges connected to a specific vertex by index
std::pair<std::vector<std::pair<int, int>>, std::vector<int>> Mesh::find_adjacent_elements(int vertex_index) {
    std::unordered_set<std::pair<int, int>, HashPair> adjacent_edges;
    std::unordered_set<int> adjacent_triangles;

    // Find triangles containing the vertex
    auto triangles_containing_vertex = vertex_to_triangles.find(vertex_index);
    if (triangles_containing_vertex != vertex_to_triangles.end()) {
        for (int triangle_index : triangles_containing_vertex->second) {
            adjacent_triangles.insert(triangle_index);
            const std::vector<int>& triangle = triangles[triangle_index];

            // Find edges directly connected to the vertex
            for (int j = 0; j < 3; ++j) {
                std::pair<int, int> edge = std::make_pair(triangle[j], triangle[(j + 1) % 3]);
                if (vertex_index == edge.first || vertex_index == edge.second) {
                    adjacent_edges.insert(std::make_pair(std::min(edge.first, edge.second), std::max(edge.first, edge.second)));
                }
            }
        }
    }

    // Convert sets to vectors
    std::vector<std::pair<int, int>> adjacent_edges_vector(adjacent_edges.begin(), adjacent_edges.end());
    std::vector<int> adjacent_triangles_vector(adjacent_triangles.begin(), adjacent_triangles.end());

    return std::make_pair(adjacent_edges_vector, adjacent_triangles_vector);
}

// Function to calculate the signed area of a triangle
double calculate_signed_area(const point_t& v1, const point_t& v2, const point_t& v3) {
    return 0.5 * (v1[0] * (v2[1] - v3[1]) + v2[0] * (v3[1] - v1[1]) + v3[0] * (v1[1] - v2[1]));
}

// Find neighboring vertices by vertex index
void Mesh::find_vertex_connectivity(int vertex_index, std::vector<int> &neighborList, std::vector<int> &neighborMap) {
    std::unordered_set<int> neighboring_vertices;

    // Find triangles containing the vertex
    auto triangles_containing_vertex = vertex_to_triangles.find(vertex_index);
    if (triangles_containing_vertex != vertex_to_triangles.end()) {
        // Iterate over each triangle that contains the vertex
        for (int triangle_index : triangles_containing_vertex->second) {
            const std::vector<int>& triangle = triangles[triangle_index];

            // Find the other vertices in the same triangle
            for (int j = 0; j < 3; ++j) {
                if (triangle[j] != vertex_index) {
                    neighboring_vertices.insert(triangle[j]);
                }
            }
        }

        // Convert set to vector to return unique neighboring vertices
        neighborList = std::vector<int>(neighboring_vertices.begin(), neighboring_vertices.end());

        // Iterate over each triangle that contains the vertex
        for (int triangle_index : triangles_containing_vertex->second) {
            const std::vector<int>& triangle = triangles[triangle_index];

            // Find the other vertices in the same triangle
            std::vector<int> other_vertices;
            for (int j = 0; j < 3; ++j) {
                if (triangle[j] != vertex_index) {
                    // Add vertex to neighborMap
                    for (int i = 0; i < neighborList.size(); i++) {
                        if (neighborList[i] == triangle[j]) {
                            neighborMap.push_back(i);
                            other_vertices.push_back(neighborList[i]);
                        }
                    }
                }
            }

            // If we have two other vertices in the triangle, calculate the area
            if (other_vertices.size() == 2) {
                int v1_idx = other_vertices[0];
                int v2_idx = other_vertices[1];

                // Calculate the signed area of the triangle formed by the current vertex and its two neighbors
                double area = calculate_signed_area(
                    source_points[vertex_index],  // Current vertex
                    source_points[v1_idx],        // First neighbor
                    source_points[v2_idx]         // Second neighbor
                );

                // If the area is negative, swap the two neighbor indices to correct the orientation
                if (area < 0.0f) {
                    std::swap(neighborMap[neighborMap.size() - 1], neighborMap[neighborMap.size() - 2]);
                }
            }
        }
    }
}

// Compute the angle between two vectors
double compute_angle(const point_t &v1, const point_t &v2) {
    double dot_product = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    double length_v1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    double length_v2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    return std::acos(dot_product / (length_v1 * length_v2));
}

// Find the normal vector for the given vertex
point_t Mesh::calculate_vertex_normal(std::vector<point_t> &points, int vertex_index) {
    point_t final_normal = {0.0, 0.0, 0.0};
    double total_weight = 0.0;

    // Find triangles connected to the vertex
    auto triangles_containing_vertex = vertex_to_triangles.find(vertex_index);
    if (triangles_containing_vertex != vertex_to_triangles.end()) {
        // Iterate over each triangle containing the vertex
        for (int triangle_index : triangles_containing_vertex->second) {
            const std::vector<int>& triangle = triangles[triangle_index];

            // Get the three vertices of the triangle
            point_t p0 = points[triangle[0]];
            point_t p1 = points[triangle[1]];
            point_t p2 = points[triangle[2]];

            // Determine which vertex in the triangle corresponds to the vertex_index
            point_t vertex, v1, v2;
            if (triangle[0] == vertex_index) {
                vertex = p0; v1 = p1; v2 = p2;
            } else if (triangle[1] == vertex_index) {
                vertex = p1; v1 = p2; v2 = p0;
            } else {
                vertex = p2; v1 = p0; v2 = p1;
            }

            // Compute two edges of the triangle
            point_t edge1 = {v1[0] - vertex[0], v1[1] - vertex[1], v1[2] - vertex[2]};
            point_t edge2 = {v2[0] - vertex[0], v2[1] - vertex[1], v2[2] - vertex[2]};

            // Compute the normal of the triangle (cross product of two edges)
            point_t normal = cross_product(edge1, edge2);
            normal = normalize(normal);

            // Compute the angle at the given vertex (between edge1 and edge2)
            double angle = compute_angle(edge1, edge2);

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
        return normalize(final_normal);
    }

    // If no triangles are connected, return a zero vector
    return {0.0, 0.0, 0.0};
}

// Function to calculate angle between two points with respect to a reference point
/*double calculateAngle(const point_t& a, const point_t& reference) {
    return std::atan2(a[1] - reference[1], a[0] - reference[0]);
}*/

// Build a dual cell from a given vertex
std::vector<point_t> Mesh::get_barycentric_dual_cell(int point, std::vector<std::vector<double>> &points) {
    std::vector<std::pair<int, int>> adjacent_edges;
    std::vector<int> adjacent_triangles;

    // Find adjacent edges and triangles
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>> adjacent_elements = find_adjacent_elements(point);

    adjacent_edges = adjacent_elements.first;
    adjacent_triangles = adjacent_elements.second;

    // Store dual cell vertices
    std::vector<point_t> dual_points;

    // Append triangle centroids to the dual cell vertices with angles
    for (int i = 0; i < adjacent_triangles.size(); i++) {
        int triangle_index = adjacent_triangles[i];

        const std::vector<int>& triangle = this->triangles[triangle_index];
        const point_t& p1 = points[triangle[0]];
        const point_t& p2 = points[triangle[1]];
        const point_t& p3 = points[triangle[2]];

        // Centroid of the triangle
        double centroid_x = (p1[0] + p2[0] + p3[0]) / 3.0;
        double centroid_y = (p1[1] + p2[1] + p3[1]) / 3.0;

        dual_points.push_back({centroid_x, centroid_y});
    }

    // Append edge centroids to the dual cell vertices with angles
    for (int i = 0; i < adjacent_edges.size(); i++) {
        std::pair<int, int> edge = adjacent_edges[i];

        const point_t& p1 = points[edge.first];
        const point_t& p2 = points[edge.second];

        // Centroid of the edge
        double centroid_x = (p1[0] + p2[0]) / 2.0;
        double centroid_y = (p1[1] + p2[1]) / 2.0;

        dual_points.push_back({centroid_x, centroid_y});
    }

    double epsilon = std::numeric_limits<double>::epsilon();

    std::vector<point_t> dual_points_copy;
    dual_points_copy.resize(dual_points.size());

    // add the vertex itself to the dual vertices if there are less than 4 adjacent triangles
    if (adjacent_triangles.size() <= 3) {
        for (int i=0; i<dual_points_copy.size(); i++) {
            dual_points_copy[i] = dual_points[i];
        }

        // add point that is slightly moved away from the other vertices instead of the point itself
        //point_t current_centroid = calculate_polygon_centroid(dual_points_copy);
        //point_t new_point = {current_centroid[0]*epsilon, current_centroid[1]*epsilon};
        //dual_points_copy.push_back(new_point);
        dual_points_copy.push_back(points[point]);

        std::sort(dual_points.begin(), dual_points.end(), [&](const point_t& a, const point_t& b) {
            double angle_a = std::atan2(a[1] - points[point][1], a[0] - points[point][0]);
            double angle_b = std::atan2(b[1] - points[point][1], b[0] - points[point][0]);
            return angle_a < angle_b;
        });

        // Create a vector of indices
        /*std::vector<int> indices(dual_points_copy.size());
        for (int i = 0; i < dual_points_copy.size(); ++i) {
            indices[i] = i;
        }

        // Sort the indices based on the angles
        std::sort(indices.begin(), indices.end(), [&](int a, int b) {
            double angle_a = std::atan2(dual_points_copy[a][1] - points[point][1], dual_points_copy[a][0] - points[point][0]);
            double angle_b = std::atan2(dual_points_copy[b][1] - points[point][1], dual_points_copy[b][0] - points[point][0]);
            return angle_a < angle_b;
        });

        for (int i = 0; i < dual_points_copy.size(); ++i) {
            //if (i == dual_points_copy.size() - 1) {
            //    dual_points[i] = points[point];
            //} else {
                dual_points[i] = dual_points_copy[indices[i]];
            //}
        }*/
    } else {
        // Sort triangles and edges based on angles with respect to the centroid
        std::sort(dual_points.begin(), dual_points.end(), [&](const point_t& a, const point_t& b) {
            double angle_a = std::atan2(a[1] - points[point][1], a[0] - points[point][0]);
            double angle_b = std::atan2(b[1] - points[point][1], b[0] - points[point][0]);
            return angle_a < angle_b;
        });
    }

    return dual_points;
}

// build barycentric dual mesh for the source mesh
void Mesh::build_source_dual_cells(std::vector<std::vector<point_t>> &cells) {
    for (int i=0; i<this->source_points.size(); i++) {
        std::vector<point_t> cell = get_barycentric_dual_cell(i, this->source_points);
        cells.push_back(cell);
    }
}

// build barycentric dual mesh for the target mesh
void Mesh::build_target_dual_cells(std::vector<std::vector<point_t>> &cells) {
    for (int i=0; i<this->target_points.size(); i++) {
        std::vector<point_t> cell = get_barycentric_dual_cell(i, this->target_points);
        cells.push_back(cell);
    }
}

// build barycentric dual mesh for the target mesh
void Mesh::build_circular_target_dual_cells(std::vector<std::vector<point_t>> &cells) {
    std::vector<point_t> temp_points = circular_transform(source_points);
    for (int i=0; i<temp_points.size(); i++) {
        std::vector<point_t> cell = get_barycentric_dual_cell(i, temp_points);
        cells.push_back(cell);
    }
}

// interpolate target mesh into a rectangular grid
std::vector<std::vector<double>> Mesh::interpolate_raster_target(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss) {
    build_target_bvh(5, 30);
    
    // Generate x and y vectors
    std::vector<double> x(res_x);
    std::vector<double> y(res_y);

    double epsilon = 1e-8;//std::numeric_limits<float>::epsilon();

    for (int i = 0; i < res_x; ++i) {
        //x[i] = ((static_cast<double>(i) + 1) / res_x) * width - (1 * width) / (res_x);
        x[i] = static_cast<double>(i) * (width - epsilon) / (res_x - 1) + 0.5 * epsilon;
        //x[i] = (static_cast<double>(i) + 1) * width / (res_x);
        //x[i] = x[i] - 0.000001 / res_x;
    }

    for (int i = 0; i < res_y; ++i) {
        //y[i] = static_cast<double>(i) * height / (res_y - 1);
        y[i] = static_cast<double>(i) * (height - epsilon) / (res_y - 1) + 0.5 * epsilon;
        //y[i] = y[i] - 0.000001 / res_y;
    }

    // Generate raster
    std::vector<std::vector<double>> raster;
    for (int i = 0; i < res_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < res_x; ++j) {
            point_t point = {x[j], y[i]};
            Hit hit;
            bool intersection = false;
            target_bvh->query(point, hit, intersection);
            if (intersection) {
                std::vector<double> vertex_values;
                vertex_values.reserve(3);
                for (int k = 0; k < 3; ++k)
                    vertex_values.push_back(errors[triangles[hit.face_id][k]]);
                double interpolation = 
                    vertex_values[0]*hit.barycentric_coords[0] + 
                    vertex_values[1]*hit.barycentric_coords[1] + 
                    vertex_values[2]*hit.barycentric_coords[2];
                row.push_back(interpolation);
                triangle_miss = false;
            } else {
                printf("interpolation miss!\r\n");
                printf("x: %f, y: %f\r\n", point[0], point[1]);
                exit(0);
                triangle_miss = true;
                row.push_back(NAN);
            }
        }
        raster.push_back(row);
    }

    return raster;
}

// interpolate target mesh into a rectangular grid
std::vector<std::vector<double>> Mesh::interpolate_raster_source(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss) {
    build_source_bvh(5, 30);
    
    // Generate x and y vectors
    std::vector<double> x(res_x);
    std::vector<double> y(res_y);

    double epsilon = 1e-8;//std::numeric_limits<float>::epsilon();

    for (int i = 0; i < res_x; ++i) {
        //x[i] = ((static_cast<double>(i) + 1) / res_x) * width - (1 * width) / (res_x);
        x[i] = static_cast<double>(i) * (width - epsilon) / (res_x - 1) + 0.5 * epsilon;
        //x[i] = (static_cast<double>(i) + 1) * width / (res_x);
        //x[i] = x[i] - 0.000001 / res_x;
    }

    for (int i = 0; i < res_y; ++i) {
        //y[i] = static_cast<double>(i) * height / (res_y - 1);
        y[i] = static_cast<double>(i) * (height - epsilon) / (res_y - 1) + 0.5 * epsilon;
        //y[i] = y[i] - 0.000001 / res_y;
    }

    // Generate raster
    std::vector<std::vector<double>> raster;
    for (int i = 0; i < res_y; ++i) {
        std::vector<double> row;
        for (int j = 0; j < res_x; ++j) {
            point_t point = {x[j], y[i]};
            Hit hit;
            bool intersection = false;
            source_bvh->query(point, hit, intersection);
            if (intersection) {
                std::vector<double> vertex_values;
                vertex_values.reserve(3);
                for (int k = 0; k < 3; ++k)
                    vertex_values.push_back(errors[triangles[hit.face_id][k]]);
                double interpolation = 
                    vertex_values[0]*hit.barycentric_coords[0] + 
                    vertex_values[1]*hit.barycentric_coords[1] + 
                    vertex_values[2]*hit.barycentric_coords[2];
                row.push_back(interpolation);
                triangle_miss = false;
            } else {
                printf("interpolation miss!\r\n");
                printf("x: %f, y: %f\r\n", point[0], point[1]);
                exit(0);
                triangle_miss = true;
                row.push_back(NAN);
            }
        }
        raster.push_back(row);
    }

    return raster;
}

// exports the inverted transport map as svg (mesh where its density distrbution is dependent on the image intensity)
std::vector<std::vector<double>> Mesh::calculate_inverted_transport_map() {
    build_target_bvh(5, 30);

    double epsilon = 1e-8;//std::numeric_limits<float>::epsilon();

    std::vector<std::vector<double>> inverted_points;
    for (int i=0; i<this->source_points.size(); ++i) {
        
        std::vector<double> point = {
            epsilon + this->source_points[i][0] * ((width - 2*epsilon) / width), 
            epsilon + this->source_points[i][1] * ((height - 2*epsilon) / height), 
            this->source_points[i][2]
        };

        Hit hit;
        bool intersection = false;
        target_bvh->query(point, hit, intersection);
        if (intersection) {
            std::vector<std::vector<double>> vertex_values;
            vertex_values.push_back(source_points[this->triangles[hit.face_id][0]]);
            vertex_values.push_back(source_points[this->triangles[hit.face_id][1]]);
            vertex_values.push_back(source_points[this->triangles[hit.face_id][2]]);
            
            double interpolation_x = 
                vertex_values[0][0]*hit.barycentric_coords[0] + 
                vertex_values[1][0]*hit.barycentric_coords[1] + 
                vertex_values[2][0]*hit.barycentric_coords[2];

            double interpolation_y = 
                vertex_values[0][1]*hit.barycentric_coords[0] + 
                vertex_values[1][1]*hit.barycentric_coords[1] + 
                vertex_values[2][1]*hit.barycentric_coords[2];

            int y = i / res_x;
            int x = i % res_x;

            if (x == 0 && y == 0) {
                inverted_points.push_back({0, 0});
            } else if (x == 0 && y == res_y - 1) {
                inverted_points.push_back({0, height});
            } else if (x == res_x - 1 && y == 0) {
                inverted_points.push_back({width, 0});
            } else if (x == res_x - 1 && y == res_y - 1) {
                inverted_points.push_back({width, height});
            } else if (x == 0 && (y != 0 && y != res_y - 1)) {
                inverted_points.push_back({0, interpolation_y});
            } else if (x == res_x - 1 && (y != 0 && y != res_y - 1)) {
                inverted_points.push_back({width, interpolation_y});
            } else if (y == 0 && (x != 0 && x != res_x - 1)) {
                inverted_points.push_back({interpolation_x, 0});
            } else if (y == res_y - 1 && (x != 0 && x != res_x - 1)) {
                inverted_points.push_back({interpolation_x, height});
            } else if (x != 0 && x != res_x - 1 && y != 0 && y != res_y - 1) {
                inverted_points.push_back({interpolation_x, interpolation_y});
            }

            //inverted_points.push_back({interpolation_x, interpolation_y});
        }
    }

    return inverted_points;
}

void Mesh::calculate_and_export_inverted_transport_map(std::string filename, double stroke_width) {
    std::vector<std::vector<double>> inverted_points = calculate_inverted_transport_map();
    export_grid_to_svg(inverted_points, this->width, this->height, this->res_x, this->res_y, filename, stroke_width);
}

// find the maximum delta_t given a triangle and the vertex velocities where the triangle will collapse
std::vector<double> find_t(const point_t& p1, const point_t& p2, const point_t& p3,
                              const point_t& dp1, const point_t& dp2, const point_t& dp3) {
    double x1 = p2[0] - p1[0], y1 = p2[1] - p1[1];
    double x2 = p3[0] - p1[0], y2 = p3[1] - p1[1];
    double u1 = dp2[0] - dp1[0], v1 = dp2[1] - dp1[1];
    double u2 = dp3[0] - dp1[0], v2 = dp3[1] - dp1[1];

    double a = u1 * v2 - u2 * v1;
    double b = x1 * v1 - y1 * u1 - x2 * v1 + y2 * u1;
    double c = x1 * y2 - x2 * y1;

    std::vector<double> result = {-123.0, -123.0};  // Initialize with invalid values

    if (a != 0) {
        double quotient = b * b - (4 * a) * c;
        if (quotient >= 0) {
            double d = std::sqrt(quotient);
            double t1 = (-b - d) / (2 * a);
            double t2 = (-b + d) / (2 * a);

            // Both t1 and t2 are valid
            result[0] = t1;
            result[1] = t2;
        }
    }

    return result;
}

// Function to find the minimum delta_t values for each triangle
double Mesh::find_min_delta_t(const std::vector<std::vector<double>>& velocities) {
    std::vector<double> min_t_values;

    //for (const auto& triangle : this->triangles) {
    for (int tri=0; tri<this->triangles.size(); tri++) {
        std::vector<std::vector<double>> t_values;
        // try every combination
        t_values.push_back(find_t(target_points[triangles[tri][0]], target_points[triangles[tri][1]], target_points[triangles[tri][2]], velocities[triangles[tri][0]], velocities[triangles[tri][1]], velocities[triangles[tri][2]]));
        t_values.push_back(find_t(target_points[triangles[tri][1]], target_points[triangles[tri][0]], target_points[triangles[tri][2]], velocities[triangles[tri][1]], velocities[triangles[tri][0]], velocities[triangles[tri][2]]));
        t_values.push_back(find_t(target_points[triangles[tri][2]], target_points[triangles[tri][0]], target_points[triangles[tri][1]], velocities[triangles[tri][2]], velocities[triangles[tri][0]], velocities[triangles[tri][1]]));
        t_values.push_back(find_t(target_points[triangles[tri][0]], target_points[triangles[tri][2]], target_points[triangles[tri][1]], velocities[triangles[tri][0]], velocities[triangles[tri][2]], velocities[triangles[tri][1]]));
        t_values.push_back(find_t(target_points[triangles[tri][2]], target_points[triangles[tri][1]], target_points[triangles[tri][0]], velocities[triangles[tri][2]], velocities[triangles[tri][1]], velocities[triangles[tri][0]]));
        t_values.push_back(find_t(target_points[triangles[tri][1]], target_points[triangles[tri][2]], target_points[triangles[tri][0]], velocities[triangles[tri][1]], velocities[triangles[tri][2]], velocities[triangles[tri][0]]));

        // Ignore negative or zero values
        std::vector<double> valid_t_values;
        for (int i=0; i<t_values.size(); i++) {
            //printf("delta_t[0] = %f, delta_t[1] = %f\r\n", t_values[i][0], t_values[i][1]);
            if (t_values[i][0] > 0) {
                valid_t_values.push_back(t_values[i][0]);
            }
            if (t_values[i][1] > 0) {
                valid_t_values.push_back(t_values[i][1]);
            }
        }

        if (!valid_t_values.empty()) {
            min_t_values.push_back(*std::min_element(valid_t_values.begin(), valid_t_values.end()));
        } else {
            min_t_values.push_back(std::numeric_limits<double>::infinity());
        }
    }

    // Calculate the minimum of the minimum delta_t values
    return *std::min_element(min_t_values.begin(), min_t_values.end());
}

// Function to update points based on velocities and minimum delta_t
double Mesh::step_grid(const std::vector<double>& dfx, const std::vector<double>& dfy, double step_size) {
    std::vector<std::vector<double>> velocities;

    // Populate velocities and delta_t
    for (int i = 0; i < target_points.size(); i++) {
        int y = i / res_x;
        int x = i % res_x;

        if (x == 0 && y == 0) {
            velocities.push_back({0, 0});
        } else if (x == 0 && y == res_y - 1) {
            velocities.push_back({0, 0});
        } else if (x == res_x - 1 && y == 0) {
            velocities.push_back({0, 0});
        } else if (x == res_x - 1 && y == res_y - 1) {
            velocities.push_back({0, 0});
        } else if (x == 0 && (y != 0 && y != res_y - 1)) {
            velocities.push_back({0, dfy[i]});
        } else if (x == res_x - 1 && (y != 0 && y != res_y - 1)) {
            velocities.push_back({0, dfy[i]});
        } else if (y == 0 && (x != 0 && x != res_x - 1)) {
            velocities.push_back({dfx[i], 0});
        } else if (y == res_y - 1 && (x != 0 && x != res_x - 1)) {
            velocities.push_back({dfx[i], 0});
        } else if (x != 0 && x != res_x - 1 && y != 0 && y != res_y - 1) {
            velocities.push_back({dfx[i], dfy[i]});
        }
    }

    // Apply regularization to the velocity field
    /*for (size_t i = 0; i < velocities.size(); ++i) {
        double regularization_term = 0.002 * std::sqrt(std::pow(velocities[i][0], 2) + std::pow(velocities[i][1], 2));
        velocities[i][0] -= regularization_term * velocities[i][0];
        velocities[i][1] -= regularization_term * velocities[i][1];
    }*/

    //double min_t = find_min_delta_t(velocities);
    double min_t = (width/res_x);
    //std::cout << "min_t = " << min_t << std::endl;

    // Move vertices along the gradient
    for (size_t i = 0; i < target_points.size(); ++i) {
        target_points[i][0] += velocities[i][0] * min_t * step_size;
        target_points[i][1] += velocities[i][1] * min_t * step_size;
    }

    return min_t;
}

void Mesh::laplacian_smoothing(std::vector<std::vector<double>> &points, double smoothing_factor) {
    std::vector<std::vector<double>> points_copy;
    for (int i = 0; i < points.size(); i++) {
        int y = i / res_x;
        int x = i % res_x;

        std::vector<double> new_point = points[y * res_x + x];

        if (x == 0 && y == 0) {
            points_copy.push_back(points[i]);
        } else if (x == 0 && y == res_y - 1) {
            points_copy.push_back(points[i]);
        } else if (x == res_x - 1 && y == 0) {
            points_copy.push_back(points[i]);
        } else if (x == res_x - 1 && y == res_y - 1) {
            points_copy.push_back(points[i]);
        } else if (x == 0 && (y != 0 && y != res_y - 1)) {
            new_point[1] = 0.0f;
            new_point[1] += points[y * res_x + (x + 1)][1];
            new_point[1] += points[(y + 1) * res_x + x][1];
            new_point[1] += points[(y - 1) * res_x + x][1];
            new_point[1] /= 3.0f;
            //points[y * res_x + x][0] = new_point[0] = 0.0f;
            points_copy.push_back(new_point);
        } else if (x == res_x - 1 && (y != 0 && y != res_y - 1)) {
            new_point[1] = 0.0f;
            new_point[1] += points[y * res_x + (x - 1)][1];
            new_point[1] += points[(y - 1) * res_x + x][1];
            new_point[1] += points[(y + 1) * res_x + x][1];
            new_point[1] /= 3.0f;
            //points[y * res_x + x][0] = new_point[0] = width;
            points_copy.push_back(new_point);
        } else if (y == 0 && (x != 0 && x != res_x - 1)) {
            new_point[0] = 0.0f;
            new_point[0] += points[y * res_x + (x - 1)][0];
            new_point[0] += points[y * res_x + (x + 1)][0];
            new_point[0] += points[(y + 1) * res_x + x][0];
            new_point[0] /= 3.0f;
            //points[y * res_x + x][1] = new_point[1] = 0.0f;
            points_copy.push_back(new_point);
        } else if (y == res_y - 1 && (x != 0 && x != res_x - 1)) {
            new_point[0] = 0.0f;
            new_point[0] += points[y * res_x + (x - 1)][0];
            new_point[0] += points[(y - 1) * res_x + x][0];
            new_point[0] += points[y * res_x + (x + 1)][0];
            new_point[0] /= 3.0f;
            //points[y * res_x + x][1] = new_point[1] = height;
            points_copy.push_back(new_point);
        } else if (x != 0 && x != res_x - 1 && y != 0 && y != res_y - 1) {
            new_point[0] = 0.0f;
            new_point[1] = 0.0f;

            new_point[0] += points[y * res_x + (x - 1)][0];
            new_point[0] += points[y * res_x + (x + 1)][0];
            new_point[0] += points[(y - 1) * res_x + x][0];
            new_point[0] += points[(y + 1) * res_x + x][0];
            new_point[1] += points[y * res_x + (x - 1)][1];
            new_point[1] += points[y * res_x + (x + 1)][1];
            new_point[1] += points[(y - 1) * res_x + x][1];
            new_point[1] += points[(y + 1) * res_x + x][1];

            new_point[0] /= 4.0f;
            new_point[1] /= 4.0f;

            points_copy.push_back(new_point);
        }
    }

    for (int i = 0; i < points.size(); i++) {
        points[i][0] = points_copy[i][0] * smoothing_factor + points[i][0] * (1.0f - smoothing_factor);
        points[i][1] = points_copy[i][1] * smoothing_factor + points[i][1] * (1.0f - smoothing_factor);
    }
}

void Mesh::export_paramererization_to_svg(std::string filename, double stroke_width) {
    export_grid_to_svg(this->target_points, this->width, this->height, this->res_x, this->res_y, filename, stroke_width);
}

/*std::vector<std::vector<double>> Mesh::calculate_refractive_normals(double focal_len, double refractive_index) {
    std::vector<double> x_normals;
    std::vector<double> y_normals;

    for (int i=0; i<this->target_points.size(); i++) {
        double dx = this->source_points[i][0] - this->target_points[i][0];
        double dy = this->source_points[i][1] - this->target_points[i][1];
        double dz = (this->source_points[i][2] + focal_len) + this->target_points[i][2];

        double k = refractive_index * std::sqrt(dx * dx + dy * dy + dz * dz) - dz;

        x_normals.push_back((1.0f / k) * dx);
        y_normals.push_back((1.0f / k) * dy);
    }

    return {x_normals, y_normals};
}*/

// calculate target vertex normals for refractive caustics
std::vector<std::vector<double>> Mesh::calculate_refractive_normals(double focal_len, double refractive_index) {
    std::vector<double> x_normals;
    std::vector<double> y_normals;
    std::vector<double> z_normals;

    // n = (t - µi) / ||(t - µi)||
    // where:
    // n = surface normal
    // t = transmitted ray normal
    // i = incident ray normal
    // µ = refraction coefficient

    //std::vector<double> point_src = {0, 0, -20.0f};

    std::vector<double> incident = {0.0f, 0.0f, 1.0f};

    for (int i=0; i<this->target_points.size(); i++) {
        std::vector<double> transmitted = {
            this->source_points[i][0] - this->target_points[i][0],
            this->source_points[i][1] - this->target_points[i][1],
            this->source_points[i][2] + this->target_points[i][2] + focal_len
        };

        //std::vector<double> incident = {0.0f, 0.0f, 0.0f};
        //incident[0] = this->target_points[i][0] - point_src[0];
        //incident[1] = this->target_points[i][1] - point_src[1];
        //incident[2] = this->target_points[i][2] - point_src[2];

        transmitted = normalize(transmitted);
        incident = normalize(incident);

        // t - µi
        double x_normal = transmitted[0] - incident[0] * refractive_index;
        double y_normal = transmitted[1] - incident[1] * refractive_index;
        double z_normal = transmitted[2] - incident[2] * refractive_index;

        // (t - µi) / ||(t - µi)||
        x_normals.push_back(x_normal / z_normal);
        y_normals.push_back(y_normal / z_normal);
        z_normals.push_back(z_normal / z_normal);
    }

    return {x_normals, y_normals, z_normals};
}

// calculate target vertex normals for refractive caustics
std::vector<std::vector<double>> Mesh::calculate_refractive_normals_uniform(double focal_len, double refractive_index) {
    std::vector<std::vector<double>> inverted_points = calculate_inverted_transport_map();
    
    std::vector<double> x_normals;
    std::vector<double> y_normals;
    std::vector<double> z_normals;

    // n = (t - µi) / ||(t - µi)||
    // where:
    // n = surface normal
    // t = transmitted ray normal
    // i = incident ray normal
    // µ = refraction coefficient

    //std::vector<double> point_src = {0, 0, -20.0f};

    std::vector<double> incident = {0.0f, 0.0f, 1.0f};

    for (int i=0; i<this->target_points.size(); i++) {
        std::vector<double> transmitted = {
            inverted_points[i][0] - this->source_points[i][0],
            inverted_points[i][1] - this->source_points[i][1],
            0 - this->source_points[i][2]  + focal_len
        };

        //std::vector<double> incident = {0.0f, 0.0f, 0.0f};
        //incident[0] = this->target_points[i][0] - point_src[0];
        //incident[1] = this->target_points[i][1] - point_src[1];
        //incident[2] = this->target_points[i][2] - point_src[2];

        transmitted = normalize(transmitted);
        incident = normalize(incident);

        // t - µi
        double x_normal = transmitted[0] - incident[0] * refractive_index;
        double y_normal = transmitted[1] - incident[1] * refractive_index;
        double z_normal = transmitted[2] - incident[2] * refractive_index;
        
        std::vector<double> normal = {x_normal, y_normal, z_normal};

        normal = normalize(normal);

        // (t - µi) / ||(t - µi)||
        x_normals.push_back(normal[0]);
        y_normals.push_back(normal[1]);
        z_normals.push_back(normal[2]);
    }

    return {x_normals, y_normals, z_normals};
}

double Mesh::set_source_heights(std::vector<double> heights) {
        // Find maximum height
    double max_h = 0;
    for (int i=0; i<heights.size(); i++) {
        double h = heights[i];

        if (max_h > h) {
            max_h = h;
        }
    }

    double update_sum = 0.0f;
    for (int i=0; i<heights.size(); i++) {
        heights[i] -= max_h;
        update_sum += (heights[i] - this->source_points[i][2]) * (heights[i] - this->source_points[i][2]);
        this->source_points[i][2] = heights[i];
    }
    return update_sum;
}

double Mesh::set_target_heights(std::vector<double> heights) {
    double update_sum = 0.0f;
    for (int i=0; i<heights.size(); i++) {
        update_sum += (heights[i] - this->target_points[i][2]) * (heights[i] - this->target_points[i][2]);
        this->target_points[i][2] = heights[i];
    }
    return update_sum;
}

void Mesh::save_solid_obj_source(double thickness, const std::string& filename) {
    save_solid_obj(this->source_points, this->source_points, this->triangles, thickness, this->width, this->height, this->res_x, this->res_y, filename);
}

void Mesh::save_solid_obj_target(double thickness, const std::string& filename) {
    save_solid_obj(this->target_points, this->source_points, this->triangles, thickness, this->width, this->height, this->res_x, this->res_y, filename);
}

void Mesh::get_vertex_neighbor_ids(int vertex_id, int &left_vertex, int &right_vertex, int &top_vertex, int &bottom_vertex) {
    int y = vertex_id / res_x;
    int x = vertex_id % res_x;

    if (x != 0) {
        left_vertex = (y) * res_x + (x - 1);
    } else {
        left_vertex = -1;
    }

    if (y != 0) {
        top_vertex = (y - 1) * res_x + (x);
    } else {
        top_vertex = -1;
    }

    if (x != res_x - 1) {
        right_vertex = (y) * res_x + (x + 1);
    } else {
        right_vertex = -1;
    }

    if (y != res_y - 1) {
        bottom_vertex = (y + 1) * res_x + (x);
    } else {
        bottom_vertex = -1;
    }
}

// Function to calculate the approximate vertex normal
/*std::vector<double> Mesh::calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index) {
    std::vector<double> avg_normal = {0.0, 0.0, 0.0}; // Initialize normal to zero vector
    
    int left_vtx = 0;
    int right_vtx = 0;
    int top_vtx = 0;
    int bot_vtx = 0;

    //printf("aa\r\n");
    
    get_vertex_neighbor_ids(vertex_index, left_vtx, right_vtx, top_vtx, bot_vtx);

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
    //double magnitude = sqrt(avg_normal[0] * avg_normal[0] + avg_normal[1] * avg_normal[1] + avg_normal[2] * avg_normal[2]);

    // Avoid division by zero
    //if (magnitude > 1e-12) {
    avg_normal[0] *= -1;
    avg_normal[1] *= -1;

    //printf("ae\r\n"); fflush(stdout);

    return avg_normal;
}*/