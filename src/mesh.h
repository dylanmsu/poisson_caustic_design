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

    public:
        Mesh(double width, double height, int res_x, int res_y);
        ~Mesh();

        void export_to_svg(std::string filename, double stroke_width);
        polygon_t get_barycentric_dual_cell(int point);
        std::vector<polygon_t> build_dual_cells();
};

#endif
