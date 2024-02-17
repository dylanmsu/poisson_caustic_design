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

    // Duplicate mesh points
    for (int i=0; i<this->target_points.size(); i++) {
        this->source_points.push_back(this->target_points[i]);
    }

    // Create instance of the bvh class used for interpolation
    target_bvh = new Bvh(triangles, target_points);
}

Mesh::~Mesh()
{
    delete(target_bvh);
}

// Build the BVH tree for the target mesh
void Mesh::build_bvh(int targetCellSize, int maxDepth) {
    target_bvh->build(targetCellSize, maxDepth);
}

// generates a structured triangulation used for the parameterization
void Mesh::generate_structured_mesh(int nx, int ny, double width, double height, std::vector<std::vector<int>> &triangles, std::vector<point_t> &points) {
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
std::vector<point_t> Mesh::circular_transform(std::vector<point_t> input_points) {
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
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)height / (double)width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
    
    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Draw polygons
    for (const auto& polygon : this->triangles) {
        std::vector<point_t> poly_points;
        for (const auto& point_idx : polygon) {
            poly_points.push_back(this->target_points[point_idx]);
        }

        std::string path_str = "M";
        for (std::size_t j = 0; j < poly_points.size(); ++j) {
            const auto& point = poly_points[j];
            path_str += std::to_string((point[0] / (double)width) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)height) * 1000.0f * ((double)height / (double)width));

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

// build mapping from vertices to connected triangles -> used for creating dual cells
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

// Function to calculate angle between two points with respect to a reference point
double calculateAngle(const point_t& a, const point_t& reference) {
    return std::atan2(a[1] - reference[1], a[0] - reference[0]);
}

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

// interpolate target mesh into a rectangular grid
std::vector<std::vector<double>> Mesh::interpolate_raster(const std::vector<double>& errors, int res_x, int res_y) {
    // Generate x and y vectors
    std::vector<double> x(res_x);
    std::vector<double> y(res_y);

    for (int i = 0; i < res_x; ++i)
        x[i] = static_cast<double>(i) * width / (res_x - 1);

    for (int i = 0; i < res_y; ++i)
        y[i] = static_cast<double>(i) * height / (res_y - 1);

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
            } else {
                row.push_back(1.0);
            }
        }
        raster.push_back(row);
    }

    return raster;
}

// exports the inverted transport map as svg (mesh where its density distrbution is dependent on the image intensity)
void Mesh::calculate_inverted_transport_map(std::string filename, double stroke_width) {
    std::vector<std::vector<double>> raster;
    for (int i=0; i<this->source_points.size(); ++i) {
        std::vector<double> point = this->source_points[i];
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

            raster.push_back({interpolation_x, interpolation_y});
        }
    }

    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)height / (double)width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    //int y = i / res_y;
    //int x = i % res_y;
    //int idx = y * this->res_x + x;

    for (int j = 0; j < this->res_y; j++) {
        std::string path_str = "M";
        for (int i = 0; i < this->res_x; i++) {
            int idx = j * this->res_x + i;

            const auto& point = raster[idx];
            path_str += std::to_string((point[0] / (double)width) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)height) * 1000.0f * ((double)height / (double)width));

            if (i < this->res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    for (int j = 0; j < this->res_y; j++) {
        std::string path_str = "M";
        for (int i = 0; i < this->res_x; i++) {
            int idx = i * this->res_x + j;

            const auto& point = raster[idx];
            path_str += std::to_string((point[0] / (double)width) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)height) * 1000.0f * ((double)height / (double)width));

            if (i < this->res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
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
        //t_values.push_back(find_t(target_points[triangles[tri][0]], target_points[triangles[tri][1]], target_points[triangles[tri][2]], velocities[triangles[tri][0]], velocities[triangles[tri][1]], velocities[triangles[tri][2]]));
        //t_values.push_back(find_t(target_points[triangles[tri][1]], target_points[triangles[tri][0]], target_points[triangles[tri][2]], velocities[triangles[tri][1]], velocities[triangles[tri][0]], velocities[triangles[tri][2]]));
        //t_values.push_back(find_t(target_points[triangles[tri][2]], target_points[triangles[tri][0]], target_points[triangles[tri][1]], velocities[triangles[tri][2]], velocities[triangles[tri][0]], velocities[triangles[tri][1]]));
        //t_values.push_back(find_t(target_points[triangles[tri][0]], target_points[triangles[tri][2]], target_points[triangles[tri][1]], velocities[triangles[tri][0]], velocities[triangles[tri][2]], velocities[triangles[tri][1]]));
        t_values.push_back(find_t(target_points[triangles[tri][2]], target_points[triangles[tri][1]], target_points[triangles[tri][0]], velocities[triangles[tri][2]], velocities[triangles[tri][1]], velocities[triangles[tri][0]]));
        //t_values.push_back(find_t(target_points[triangles[tri][1]], target_points[triangles[tri][2]], target_points[triangles[tri][0]], velocities[triangles[tri][1]], velocities[triangles[tri][2]], velocities[triangles[tri][0]]));

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
        //printf("\r\n");

        std::vector<std::vector<double>> triangle;
        triangle.push_back(target_points[triangles[tri][0]]);
        triangle.push_back(target_points[triangles[tri][1]]);
        triangle.push_back(target_points[triangles[tri][2]]);
        if (calculate_polygon_area_vec(triangle) <= 0.0f) {
            printf("negative area!\r\n");
        }

        //std::cout << valid_t_values.size() << std::endl;

        //std::cout << "e" << std::endl;

        if (!valid_t_values.empty()) {
            min_t_values.push_back(*std::min_element(valid_t_values.begin(), valid_t_values.end()));
        } else {
            min_t_values.push_back(std::numeric_limits<double>::infinity());
        }
    }

    //std::cout << "f" << std::endl;

    // Calculate the minimum of the minimum delta_t values
    return *std::min_element(min_t_values.begin(), min_t_values.end());
}

// Function to update points based on velocities and minimum delta_t
double Mesh::step_grid(const std::vector<double>& dfx, const std::vector<double>& dfy, double step_size) {
    std::vector<std::vector<double>> velocities;

    printf("populated velocities\r\n");

    // Populate velocities and delta_t
    for (int i = 0; i < target_points.size(); i++) {
        int y = i / res_y;
        int x = i % res_y;

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

    // Apply regularization to the velocity field if needed
    // for (size_t i = 0; i < velocities.size(); ++i) {
    //     double regularization_term = regularization_coeff * std::sqrt(std::pow(velocities[i][0], 2) + std::pow(velocities[i][1], 2));
    //     velocities[i][0] -= regularization_term * velocities[i][0];
    //     velocities[i][1] -= regularization_term * velocities[i][1];
    // }

    printf("finding min_delta_t\r\n");

    double min_t = find_min_delta_t(velocities);
    //std::cout << "min_t = " << min_t << std::endl;

    // Move vertices along the gradient
    for (size_t i = 0; i < target_points.size(); ++i) {
        target_points[i][0] += velocities[i][0] * min_t * step_size;
        target_points[i][1] += velocities[i][1] * min_t * step_size;
    }

    printf("stepped\r\n");

    return min_t;
}

void Mesh::export_paramererization_to_svg(std::string filename, double stroke_width) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)height / (double)width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    svg_file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    //int y = i / res_y;
    //int x = i % res_y;
    //int idx = y * this->res_x + x;

    for (int j = 0; j < this->res_y; j++) {
        std::string path_str = "M";
        for (int i = 0; i < this->res_x; i++) {
            int idx = j * this->res_x + i;

            const auto& point = target_points[idx];
            path_str += std::to_string((point[0] / (double)width) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)height) * 1000.0f * ((double)height / (double)width));

            if (i < this->res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    for (int j = 0; j < this->res_y; j++) {
        std::string path_str = "M";
        for (int i = 0; i < this->res_x; i++) {
            int idx = i * this->res_x + j;

            const auto& point = target_points[idx];
            path_str += std::to_string((point[0] / (double)width) * 1000.0f) + "," +
                        std::to_string((point[1] / (double)height) * 1000.0f * ((double)height / (double)width));

            if (i < this->res_x - 1)
                path_str += "L";
        }
        svg_file << "<path d=\"" << path_str << "\" fill=\"none\" stroke=\"black\" stroke-width=\"" << stroke_width << "\"/>\n";
    }

    // Write SVG footer
    svg_file << "</svg>\n";
    svg_file.close();
}

std::vector<std::vector<double>> Mesh::calculate_refractive_normals(double focal_len, double refractive_index) {
    std::vector<double> x_normals;
    std::vector<double> y_normals;

    for (int i=0; i<this->target_points.size(); i++) {
        double dx = this->source_points[i][0] - this->target_points[i][0];
        double dy = this->source_points[i][1] - this->target_points[i][1];
        double dz = this->source_points[i][2] - this->target_points[i][2];

        double k = refractive_index * std::sqrt((dx * dx + dy * dy) + (focal_len - dz) * (focal_len - dz)) - (focal_len - dz);

        x_normals.push_back((1.0f / k) * dx);
        y_normals.push_back((1.0f / k) * dy);
    }

    return {x_normals, y_normals};
}

void Mesh::set_source_heights(std::vector<double> heights) {
    for (int i=0; i<heights.size(); i++) {
        this->source_points[i][2] = -heights[i];
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

void Mesh::save_solid_obj(double thickness, const std::string& filename) {
    int num_points = this->source_points.size();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Find maximum height
    double min_h = 0;
    for (int i=0; i<num_points; i++) {
        double h = this->source_points[i][2];

        if (min_h > h) {
            min_h = h;
        }
    }

    file << "# Generated by the software algorithm written by Dylan Missuwe" << "\n";
    file << "# The algorithm used to create the lens is based on the paper " 
        << "\"Poisson-Based Continuous Surface Generation for Goal-Based Caustics\" " 
        << "by Yue et al (2014)" << "\n";

    // Curved mesh verts on the bottom
    for (const auto& point : this->source_points) {
        file << "v " << point[0] << " " << this->width - point[1] << " " << -(point[2] - min_h) << "\n";
    }

    // Flat mesh verts on the bottom
    for (const auto& point : this->source_points) {
        file << "v " << point[0] << " " << this->width - point[1] << " " << -thickness << "\n";
    }

    // Curved mesh triangles on the top
    for (const auto& triangle : this->triangles) {
        file << "f " << triangle[0] + 1 << " " << triangle[2] + 1 << " " << triangle[1] + 1 << "\n";
    }

    // Flat mesh triangles on the bottom
    for (const auto& triangle : this->triangles) {
        file << "f " << triangle[0] + num_points + 1 << " " << triangle[1] + num_points + 1 << " " << triangle[2] + num_points + 1 << "\n";
    }

    // Generate triangles connecting top and bottom mesh
    std::vector<int> perimeter_verts;
    find_perimeter_vertices(this->res_x, this->res_y, perimeter_verts);

    for (size_t i = 0; i < perimeter_verts.size(); ++i) {
        int top_idx = perimeter_verts[i];
        int bottom_idx = perimeter_verts[i] + num_points;
        int next_top_idx = perimeter_verts[(i + 1) % perimeter_verts.size()];
        int next_bottom_idx = perimeter_verts[(i + 1) % perimeter_verts.size()] + num_points;


        file << "f " << top_idx + 1 << " " << next_bottom_idx + 1 << " " << bottom_idx + 1 << "\n";
        file << "f " << top_idx + 1 << " " << next_top_idx + 1 << " " << next_bottom_idx + 1 << "\n";
    }

    //std::cout << "Exported model \"" << filename << "\"" << std::endl;
}
