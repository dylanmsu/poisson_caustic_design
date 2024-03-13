#ifndef MESH_H
#define MESH_H

#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "bvh.h"
#include "utils.h"

typedef std::vector<double> point_t;

struct HashPair {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Using XOR (^) to combine hashes
        return h1 ^ h2;
    }
};

class Mesh {
    private:
        std::vector<std::vector<int>> triangles;

        std::unordered_map<int, std::vector<int>> vertex_to_triangles;

        // stores the BVH tree for the target mesh
        Bvh *target_bvh;

        void generate_structured_mesh(int nx, int ny, double width, double height, std::vector<std::vector<int>> &triangles, std::vector<point_t> &points);
        void build_vertex_to_triangles();
        std::pair<std::vector<std::pair<int, int>>, std::vector<int>> find_adjacent_elements(int vertex_index);

        double find_min_delta_t(const std::vector<std::vector<double>>& velocities);

    public:
        Mesh(double width, double height, int res_x, int res_y);
        ~Mesh();

        std::vector<point_t> source_points;
        std::vector<point_t> target_points;

        double width;
        double height;

        int res_x;
        int res_y;

        void export_to_svg(std::string filename, double stroke_width);
        void export_paramererization_to_svg(std::string filename, double stroke_width);

        std::vector<std::vector<double>> get_barycentric_dual_cell(int point, std::vector<std::vector<double>> &points);
        void build_target_dual_cells(std::vector<std::vector<point_t>> &cells);
        void build_source_dual_cells(std::vector<std::vector<point_t>> &cells);

        std::vector<std::vector<double>> interpolate_raster(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss);

        double step_grid(const std::vector<double>& dfx, const std::vector<double>& dfy, double step_size);

        void build_bvh(int targetCellSize, int maxDepth);

        std::vector<std::vector<double>> calculate_refractive_normals(double focal_len, double refractive_index);

        void set_source_heights(std::vector<double> heights);
        double set_target_heights(std::vector<double> heights);

        void save_solid_obj_target(double thickness, const std::string& filename);
        void save_solid_obj_source(double thickness, const std::string& filename);

        void calculate_inverted_transport_map(std::string filename, double stroke_width);

        void build_target_dual_cells_circ(std::vector<std::vector<point_t>> &cells);

        std::vector<point_t> circular_transform(std::vector<std::vector<double>> &input_points);

        void build_circular_target_dual_cells(std::vector<std::vector<point_t>> &cells);

        void laplacian_smoothing(std::vector<std::vector<double>> &points, double smoothing_factor);

        void get_vertex_neighbor_ids(int vertex_id, int &left_vertex, int &right_vertex, int &top_vertex, int &bottom_vertex);
};

#endif
