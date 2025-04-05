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

    build_adjacency_lookups();

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

// find triangles, edges, and neighboring vertices connected to a specific vertex by index
std::tuple<std::vector<std::pair<int, int>>, std::vector<int>, std::vector<int>> Mesh::find_adjacent_elements(int vertex_index) {
    std::unordered_set<std::pair<int, int>, HashPair> adjacent_edges;
    std::unordered_set<int> adjacent_triangles;
    std::unordered_set<int> neighboring_vertices;

    // Find triangles containing the vertex
    auto triangles_containing_vertex = vertex_to_triangles.find(vertex_index);
    if (triangles_containing_vertex != vertex_to_triangles.end()) {
        for (int triangle_index : triangles_containing_vertex->second) {
            adjacent_triangles.insert(triangle_index);
            const std::vector<int>& triangle = triangles[triangle_index];

            // Find edges and neighboring vertices directly connected to the vertex
            for (int j = 0; j < 3; ++j) {
                int v1 = triangle[j];
                int v2 = triangle[(j + 1) % 3];

                // Add edge if it involves the vertex
                if (vertex_index == v1 || vertex_index == v2) {
                    adjacent_edges.insert(std::make_pair(std::min(v1, v2), std::max(v1, v2)));
                }

                // Add neighboring vertex (other vertex of the edge if it's not the input vertex)
                if (v1 == vertex_index) {
                    neighboring_vertices.insert(v2);
                } else if (v2 == vertex_index) {
                    neighboring_vertices.insert(v1);
                }
            }
        }
    }

    // Convert sets to vectors
    std::vector<std::pair<int, int>> adjacent_edges_vector(adjacent_edges.begin(), adjacent_edges.end());
    std::vector<int> adjacent_triangles_vector(adjacent_triangles.begin(), adjacent_triangles.end());
    std::vector<int> neighboring_vertices_vector(neighboring_vertices.begin(), neighboring_vertices.end());

    return std::make_tuple(adjacent_edges_vector, adjacent_triangles_vector, neighboring_vertices_vector);
}

bool Mesh::is_boundary_vertex(int vertex_index, std::vector<std::pair<int, int>>& boundary_edges) {
    std::unordered_map<std::pair<int, int>, int, HashPair> edge_triangle_count;
    for (int triangle_index : vertex_adjecent_triangles[vertex_index]) {
        const std::vector<int>& triangle = this->triangles[triangle_index];
        for (int j = 0; j < 3; ++j) {
            int v1 = triangle[j];
            int v2 = triangle[(j + 1) % 3];
            std::pair<int, int> edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            edge_triangle_count[edge]++;
        }
    }

    bool is_boundary = false;
    for (const auto& edge : vertex_adjecent_edges[vertex_index]) {
        if (edge_triangle_count[edge] == 1) { // Boundary edge
            boundary_edges.push_back(edge);
            is_boundary = true;
        }
    }

    return is_boundary;
}

void Mesh::build_adjacency_lookups() {
    for (int i = 0; i < source_points.size(); i++)
    {
        auto [adjacent_edges, adjacent_triangles, neighboring_vertices] = find_adjacent_elements(i);
        vertex_adjecent_edges.push_back(adjacent_edges);
        vertex_adjecent_triangles.push_back(adjacent_triangles);
        vertex_adjecent_vertices.push_back(neighboring_vertices);
    }

    for (int i = 0; i < source_points.size(); i++)
    {
        std::vector<std::pair<int, int>> boundary_edges;
        bool is_boundary = is_boundary_vertex(i, boundary_edges);
        vertex_is_boundary.push_back(is_boundary);
    }
    
    
}

std::vector<double> cross_v(std::vector<double> v1, std::vector<double> v2){
    std::vector<double> result(3);
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return result;
}

double dot_v(std::vector<double> a, std::vector<double> b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

std::vector<double> mult_v(double a, std::vector<double> b) {
    return {a*b[0], a*b[1], a*b[2]};
}

std::vector<double> add_v(std::vector<double> a, std::vector<double> b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::vector<double> sub_v(std::vector<double> a, std::vector<double> b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double magnitude_v(std::vector<double> a) {
    return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

double cot(const std::vector<double>& a, const std::vector<double>& b) {
    auto cross_product = cross_v(a, b);
    double cross_magnitude = magnitude_v(cross_product);

    if (cross_magnitude < 1e-12) {
        //throw std::invalid_argument("Vectors are parallel or one is a zero vector, cotangent undefined.");
        cross_magnitude = 1e-12;
    }

    return dot_v(a, b) / cross_magnitude;
}

std::vector<double> Mesh::compute_laplacian(int i) {
    std::vector<double> laplacian(vertex_adjecent_vertices[i].size(), 0.0f);

    for (int j_index = 0; j_index < vertex_adjecent_vertices[i].size(); ++j_index) {
        int j = vertex_adjecent_vertices[i][j_index];

        // Find triangles shared between `i` and `j`
        std::vector<int> shared_triangles;
        for (int triangle : vertex_adjecent_triangles[i]) {
          // Check if `j` is one of the vertices in this triangle
          const auto& vertices = triangles[triangle];
          if (std::find(vertices.begin(), vertices.end(), j) != vertices.end()) {
              shared_triangles.push_back(triangle);
          }
        }

        // Handle cases based on the number of shared triangles
        if (shared_triangles.size() == 2) {
            // Interior case: Two triangles are connected
            std::vector<int> k_vertices;
            for (int triangle : shared_triangles) {
                for (int vertex : this->triangles[triangle]) {
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

            edge1 = sub_v(this->source_points[k1], this->source_points[j]);
            edge2 = sub_v(this->source_points[k1], this->source_points[i]);
            double cot_k1 = cot(edge1, edge2);

            edge1 = sub_v(this->source_points[k2], this->source_points[j]);
            edge2 = sub_v(this->source_points[k2], this->source_points[i]);
            double cot_k2 = cot(edge1, edge2);

            laplacian[j_index] += cot_k1 * 0.5;
            laplacian[j_index] += cot_k2 * 0.5;

            //std::cout << "k1=" << k1 << ", k2=" << k2 << std::endl;

        } else if (shared_triangles.size() == 1) {
            // Boundary case: Only one triangle is connected
            int triangle = shared_triangles[0];
            int k = -1;

            // Find the single `k` vertex
            for (int vertex : this->triangles[triangle]) {
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

            edge1 = sub_v(this->source_points[k], this->source_points[j]);
            edge2 = sub_v(this->source_points[k], this->source_points[i]);
            double cot_k = cot(edge1, edge2);

            laplacian[j_index] += cot_k;

            //std::cout << "k=" << k << std::endl;

        } else {
            throw std::runtime_error("No shared triangles between i and j; invalid mesh or disconnected vertex.");
        }
    }

    return laplacian;
}

void Mesh::calculate_vertex_laplacians() {
    build_adjacency_lookups();
    vertex_laplacians.clear();
    
    for (int i = 0; i < this->source_points.size(); i++)
    {
        vertex_laplacians.push_back(compute_laplacian(i));
    }
}


// Function to calculate angle between two points with respect to a reference point
double calculateAngle(const point_t& a, const point_t& reference) {
    return std::atan2(a[1] - reference[1], a[0] - reference[0]);
}

point_t edge_centroid(point_t &point_i, point_t &point_j) {
    point_t centroid(3);
    for (int i = 0; i < 3; i++)
    {
        centroid[i] = (point_i[i] + point_j[i]) / 2.0f;
    }

    return centroid;
}

point_t triangle_centroid(point_t &point_i, point_t &point_j, point_t &point_k) {
    std::vector<double> centroid(3);
    for (int i = 0; i < 3; i++)
    {
        centroid[i] = (point_i[i] + point_j[i] + point_k[i]) / 3.0f;
    }

    return centroid;
}

polygon_t Mesh::get_triangle_quad(int vertex_idx_i, int triangle_idx, std::vector<point_t>& points) {
    std::vector<int> triangle = this->triangles[triangle_idx];
    int vertex_idx_j = -1;
    int vertex_idx_k = -1;

    for (int i = 0; i < 3; i++) {
        if (triangle[i] == vertex_idx_i) {
            vertex_idx_j = triangle[(i + 1) % 3];
            vertex_idx_k = triangle[(i + 2) % 3];
            break;
        }
    }

    point_t vertex_i = points[vertex_idx_i];
    point_t vertex_j = points[vertex_idx_j];
    point_t vertex_k = points[vertex_idx_k];

    polygon_t quad;

    quad.push_back(vertex_i);
    quad.push_back(edge_centroid(vertex_i, vertex_j));
    quad.push_back(triangle_centroid(vertex_i, vertex_j, vertex_k));
    quad.push_back(edge_centroid(vertex_i, vertex_k));

    return quad;
}

std::vector<polygon_t> Mesh::get_partitioned_barycentric_dual_cell(int v_point, std::vector<std::vector<double>>& points) {
    std::vector<polygon_t> cell;

    for (int i = 0; i < vertex_adjecent_triangles[v_point].size(); i++)
    {
        int adjacent_triangle_idx = vertex_adjecent_triangles[v_point][i];

        cell.push_back(get_triangle_quad(v_point, adjacent_triangle_idx, points));
    }
    
    return cell;
}

void Mesh::smoothMeshUntilNoFlops(std::vector<std::vector<int>> &triangles, std::vector<point_t> &points, const std::vector<std::vector<double>> &laplacian, const std::vector<std::vector<int>> &vertex_adjacent_vertices) {
    const int maxIterations = 100; // Limit iterations to prevent infinite loops
    const double epsilon = 1e-6;  // Small threshold for detecting flopped triangles
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        bool hasFloppedTriangles = false;
        std::vector<point_t> newPoints = points;
        
        for (size_t i = 0; i < points.size(); ++i) {
            // Skip boundary points
            if (vertex_is_boundary[i]) continue;
            
            point_t displacement = {0.0, 0.0, 0.0};
            double weightSum = 0.0;
            
            for (size_t j = 0; j < vertex_adjacent_vertices[i].size(); ++j) {
                int neighbor = vertex_adjacent_vertices[i][j];
                double weight = laplacian[i][j];
                
                displacement[0] += weight * (points[neighbor][0] - points[i][0]);
                displacement[1] += weight * (points[neighbor][1] - points[i][1]);
                displacement[2] += weight * (points[neighbor][2] - points[i][2]);
                weightSum += weight;
            }
            
            if (weightSum > 0.0) {
                newPoints[i][0] += 0.5 * (displacement[0] / weightSum);
                newPoints[i][1] += 0.5 * (displacement[1] / weightSum);
                newPoints[i][2] += 0.5 * (displacement[2] / weightSum);
            }
        }
        
        // Check if any triangles are still flopped
        for (const auto &tri : triangles) {

            std::vector<std::vector<double>> triangle;
            triangle.push_back(points[tri[0]]);
            triangle.push_back(points[tri[1]]);
            triangle.push_back(points[tri[2]]);

            double eps = 1e-12;

            double area = calculate_polygon_area_vec(triangle);
            if (area < epsilon) {
                hasFloppedTriangles = true;
                break;
            }
        }
        
        points = newPoints;
        if (!hasFloppedTriangles) break;
    }
}

// build barycentric dual mesh for the source mesh
void Mesh::build_source_partitioned_dual_cells(std::vector<std::vector<polygon_t>> &cells) {
    for (int i=0; i<this->source_points.size(); i++) {
        std::vector<polygon_t> cell = get_partitioned_barycentric_dual_cell(i, this->source_points);
        cells.push_back(cell);
    }
}

// build barycentric dual mesh for the target mesh
void Mesh::build_target_partitioned_dual_cells(std::vector<std::vector<polygon_t>> &cells) {
    for (int i=0; i<this->target_points.size(); i++) {
        std::vector<polygon_t> cell = get_partitioned_barycentric_dual_cell(i, this->target_points);
        cells.push_back(cell);
    }
}

std::vector<std::vector<double>> Mesh::interpolate_raster(const std::vector<double>& errors, std::vector<std::vector<double>>& points, std::vector<std::vector<int>> &triangles, int res_x, int res_y, double width, double height, bool &triangle_miss) {
    Bvh *bvh = new Bvh(triangles, points);
    bvh->build(5, 30);
    
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
            std::vector<Hit> hits;
            bool intersection = false;
            bvh->query(point, hits, intersection);

            double interpolation_sum = 0.0;

            for (int k = 0; k < hits.size(); k++)
            {
                std::vector<double> vertex_values;
                for (int l = 0; l < 3; ++l) {
                    int value_index = triangles[hits[k].face_id][l];
                    double value = errors[value_index];
                    vertex_values.push_back(value);
                }
                    
                double interpolation = 
                    vertex_values[0]*hits[k].barycentric_coords[0] + 
                    vertex_values[1]*hits[k].barycentric_coords[1] + 
                    vertex_values[2]*hits[k].barycentric_coords[2];
                
                interpolation_sum += interpolation;
                triangle_miss = false;
            }
            
            if (intersection) {
                row.push_back(interpolation_sum);
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
std::vector<std::vector<double>> Mesh::interpolate_raster_target(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss) {
    return interpolate_raster(errors, this->target_points, this->triangles, res_x, res_y, this->width, this->height, triangle_miss);
}

// interpolate target mesh into a rectangular grid
std::vector<std::vector<double>> Mesh::interpolate_raster_source(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss) {
    return interpolate_raster(errors, this->source_points, this->triangles, res_x, res_y, this->width, this->height, triangle_miss);
}

// exports the inverted transport map as svg (mesh where its density distrbution is dependent on the image intensity)
std::vector<point_t> Mesh::calculate_inverted_transport_map() {
    build_target_bvh(5, 30);

    double epsilon = 1e-8;//std::numeric_limits<float>::epsilon();

    std::vector<point_t> inverted_points;
    for (int i=0; i<this->source_points.size(); ++i) {
        
        point_t point = {
            epsilon + this->source_points[i][0] * ((width - 2*epsilon) / width), 
            epsilon + this->source_points[i][1] * ((height - 2*epsilon) / height), 
            this->source_points[i][2]
        };

        std::vector<Hit> hits;
        bool intersection = false;
        target_bvh->query(point, hits, intersection);
        if (intersection) {
            std::vector<point_t> vertex_values;
            vertex_values.push_back(source_points[this->triangles[hits[0].face_id][0]]);
            vertex_values.push_back(source_points[this->triangles[hits[0].face_id][1]]);
            vertex_values.push_back(source_points[this->triangles[hits[0].face_id][2]]);
            
            double interpolation_x = 
                vertex_values[0][0]*hits[0].barycentric_coords[0] + 
                vertex_values[1][0]*hits[0].barycentric_coords[1] + 
                vertex_values[2][0]*hits[0].barycentric_coords[2];

            double interpolation_y = 
                vertex_values[0][1]*hits[0].barycentric_coords[0] + 
                vertex_values[1][1]*hits[0].barycentric_coords[1] + 
                vertex_values[2][1]*hits[0].barycentric_coords[2];

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
    std::vector<point_t> inverted_points = calculate_inverted_transport_map();
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
double Mesh::find_min_delta_t(std::vector<std::vector<double>>& velocities) {
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
            if (t_values[i][0] > 0 && t_values[i][0] < 10e3) {
                valid_t_values.push_back(t_values[i][0]);
            }
            if (t_values[i][1] > 0 && t_values[i][1] < 10e3) {
                valid_t_values.push_back(t_values[i][1]);
            }
        }

        if (!valid_t_values.empty()) {
            min_t_values.push_back(*std::min_element(valid_t_values.begin(), valid_t_values.end()));
        } else {
            min_t_values.push_back(1);
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

void Mesh::laplacian_smoothing(std::vector<point_t> &points, double smoothing_factor) {
    std::vector<point_t> points_copy;
    for (int i = 0; i < points.size(); i++) {
        int y = i / res_x;
        int x = i % res_x;

        point_t new_point = points[y * res_x + x];

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

        // (t - µi) / ||(t - µi)||
        x_normals.push_back(x_normal / z_normal);
        y_normals.push_back(y_normal / z_normal);
        z_normals.push_back(z_normal / z_normal);
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
