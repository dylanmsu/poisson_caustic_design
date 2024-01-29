#include "mesh.h"

Mesh::Mesh(double width, double height, int res_x, int res_y)
{
    this->width = width;
    this->height = height;
    this->res_x = res_x;
    this->res_y = res_y;

    generate_structured_mesh(res_x, res_y, width, height, this->triangles, this->target_points);
    build_vertex_to_triangles();

    // duplicate points
    for (int i=0; i<this->target_points.size(); i++) {
        this->source_points.push_back(this->target_points[i]);
    }
}

Mesh::~Mesh()
{
}

void Mesh::generate_structured_mesh(int nx, int ny, double width, double height, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &points) {
    // Creating x and y vectors
    std::vector<double> x(nx);
    std::vector<double> y(ny);

    for (int i = 0; i < nx; ++i)
        x[i] = static_cast<double>(i) * width / static_cast<double>(nx - 1);

    for (int i = 0; i < ny; ++i)
        y[i] = static_cast<double>(i) * height / static_cast<double>(ny - 1);

    // Creating xv and yv matrices using meshgrid
    std::vector<std::vector<double>> xv(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> yv(ny, std::vector<double>(nx, 0.0));

    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            xv[i][j] = x[j];
            yv[i][j] = y[i];
        }
    }

    // Creating points
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            points.push_back({xv[i][j], yv[i][j], 0.0});
        }
    }

    // Create triangles
    for (int i = 0; i < ny - 1; ++i) {
        for (int j = 0; j < nx - 1; ++j) {
            int idx = i * nx + j;
            triangles.push_back({idx, idx + 1, idx + nx});
            triangles.push_back({idx + nx, idx + 1, idx + nx + 1});
        }
    }
}

void Mesh::export_to_svg(std::string filename, double stroke_width) {
    std::ofstream svg_file(filename, std::ios::out);
    if (!svg_file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    // Write SVG header
    svg_file << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
    svg_file << "<svg width=\"1000\" height=\"" << 1000.0f * ((double)height / (double)width) << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

    // Draw polygons
    for (const auto& polygon : this->triangles) {
        std::vector<std::vector<double>> poly_points;
        for (const auto& point_idx : polygon) {
            poly_points.push_back(this->target_points[point_idx]);
        }

        std::string path_str = "M";
        for (std::size_t j = 0; j < poly_points.size(); ++j) {
            const auto& point = poly_points[j];
            path_str += std::to_string((point[0] / (double)width) * 1000.0f) + "," +
                        std::to_string((1 - point[1] / (double)height) * 1000.0f * ((double)height / (double)width));

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

polygon_t Mesh::get_barycentric_dual_cell(int point) {
    std::vector<std::pair<int, int>> adjacent_edges;
    std::vector<int> adjacent_triangles;

    // Find adjacent edges and triangles
    std::pair<std::vector<std::pair<int, int>>, std::vector<int>> adjacent_elements = find_adjacent_elements(point);

    adjacent_edges = adjacent_elements.first;
    adjacent_triangles = adjacent_elements.second;

    // Calculate dual points
    polygon_t dual_points;

    for (int i=0; i<adjacent_triangles.size(); i++) {
        int triangle_index = adjacent_triangles[i];

        const std::vector<int>& triangle = this->triangles[triangle_index];
        const std::vector<double>& p1 = this->target_points[triangle[0]];
        const std::vector<double>& p2 = this->target_points[triangle[1]];
        const std::vector<double>& p3 = this->target_points[triangle[2]];

        // Centroid of the triangle
        double centroid_x = (p1[0] + p2[0] + p3[0]) / 3.0;
        double centroid_y = (p1[1] + p2[1] + p3[1]) / 3.0;

        dual_points.push_back({centroid_x, centroid_y});
    }

    for (int i=0; i<adjacent_edges.size(); i++) {
        std::pair<int, int> edge = adjacent_edges[i];

        const std::vector<double>& p1 = this->target_points[edge.first];
        const std::vector<double>& p2 = this->target_points[edge.second];

        // Centroid of the edge
        double centroid_x = (p1[0] + p2[0]) / 2.0;
        double centroid_y = (p1[1] + p2[1]) / 2.0;

        dual_points.push_back({centroid_x, centroid_y});
    }

    if (adjacent_triangles.size() <= 3) {
        dual_points.push_back(this->target_points[point]);
    }

    // Order points based on angles with respect to the centroid
    double centroid_x = 0.0;
    double centroid_y = 0.0;
    for (const auto& point : dual_points) {
        centroid_x += point[0];
        centroid_y += point[1];
    }
    centroid_x /= dual_points.size();
    centroid_y /= dual_points.size();

    std::vector<double> angles;
    for (const auto& point : dual_points) {
        double angle = std::atan2(point[1] - centroid_y, point[0] - centroid_x);
        angles.push_back(angle);
    }

    std::vector<size_t> order(dual_points.size());
    std::iota(order.begin(), order.end(), 0);

    // Sort points based on angles
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return angles[a] < angles[b];
    });

    polygon_t polygon;
    for (size_t index : order) {
        polygon.push_back(dual_points[index]);
    }

    return polygon;
}

std::vector<polygon_t> Mesh::build_dual_cells() {
    std::vector<polygon_t> cells;

    for (int i=0; i<this->target_points.size(); i++)
    {
        polygon_t cell = get_barycentric_dual_cell(i);
        cells.push_back(cell);
    }

    return cells;
}


