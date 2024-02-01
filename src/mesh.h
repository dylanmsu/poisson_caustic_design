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

typedef std::vector<std::vector<double>> polygon_t;

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
        double width;
        double height;
        int res_x;
        int res_y;

        std::vector<std::vector<int>> triangles;
        std::vector<std::vector<double>> source_points;
        std::vector<std::vector<double>> target_points;

        std::unordered_map<int, std::vector<int>> vertex_to_triangles;

        void generate_structured_mesh(int nx, int ny, double width, double height, std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &points);
        void build_vertex_to_triangles();
        std::pair<std::vector<std::pair<int, int>>, std::vector<int>> find_adjacent_elements(int vertex_index);

        double find_min_delta_t(const std::vector<std::vector<double>>& velocities);

        Bvh *bvh;

    public:
        Mesh(double width, double height, int res_x, int res_y);
        ~Mesh();

        void export_to_svg(std::string filename, double stroke_width);
        polygon_t get_barycentric_dual_cell(int point, std::vector<std::vector<double>> &points);
        void build_target_dual_cells(std::vector<polygon_t> &cells);
        void build_source_dual_cells(std::vector<polygon_t> &cells);

        std::vector<std::vector<double>> interpolate_raster(const std::vector<double>& errors, int res_x, int res_y);

        void step_grid(const std::vector<double>& dfx, const std::vector<double>& dfy, double step_size);

        void build_bvh(int targetCellSize, int maxDepth);
};

#endif
