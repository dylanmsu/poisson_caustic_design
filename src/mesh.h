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
typedef std::vector<point_t> polygon_t;

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
        Bvh *source_bvh;

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

        int find_shared_triangle(int v, const std::pair<int, int>& e1, const std::pair<int, int>& e2);

        void export_to_svg(std::string filename, double stroke_width);
        void export_paramererization_to_svg(std::string filename, double stroke_width);

        polygon_t get_triangle_quad(int vertex_idx, int triangle_idx, std::vector<point_t>& points);
        polygon_t get_barycentric_dual_cell(int point, std::vector<point_t> &points);
        std::vector<polygon_t> get_partitioned_barycentric_dual_cell(int v_point, std::vector<point_t>& points);
        
        void build_target_dual_cells(std::vector<polygon_t> &cells);
        void build_source_dual_cells(std::vector<polygon_t> &cells);

        void build_target_partitioned_dual_cells(std::vector<std::vector<polygon_t>> &cells);
        void build_source_partitioned_dual_cells(std::vector<std::vector<polygon_t>> &cells);

        std::vector<std::vector<double>> interpolate_raster_target(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss);
        std::vector<std::vector<double>> interpolate_raster_source(const std::vector<double>& errors, int res_x, int res_y, bool &triangle_miss);

        double step_grid(const std::vector<double>& dfx, const std::vector<double>& dfy, double step_size);

        void build_target_bvh(int targetCellSize, int maxDepth);
        void build_source_bvh(int targetCellSize, int maxDepth);

        std::vector<std::vector<double>> calculate_refractive_normals(double focal_len, double refractive_index);
        std::vector<std::vector<double>> calculate_refractive_normals_uniform(double focal_len, double refractive_index);

        double set_source_heights(std::vector<double> heights);
        double set_target_heights(std::vector<double> heights);

        void save_solid_obj_target(double thickness, const std::string& filename);
        void save_solid_obj_source(double thickness, const std::string& filename);
        
        std::vector<point_t> calculate_inverted_transport_map();
        void calculate_and_export_inverted_transport_map(std::string filename, double stroke_width);

        void build_target_dual_cells_circ(std::vector<polygon_t> &cells);

        std::vector<point_t> circular_transform(std::vector<point_t> &input_points);

        void build_circular_target_dual_cells(std::vector<polygon_t> &cells);

        void laplacian_smoothing(std::vector<point_t> &points, double smoothing_factor);

        void get_vertex_neighbor_ids(int vertex_id, int &left_vertex, int &right_vertex, int &top_vertex, int &bottom_vertex);
};

#endif
