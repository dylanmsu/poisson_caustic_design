#ifndef BVH_H
#define BVH_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "polygon_utils.h"

struct Node {
    double bbox_min_x;
    double bbox_min_y;
    double bbox_max_x;
    double bbox_max_y;
    int first_child_id;
    int first_face_id;
    int nb_faces;
    bool is_leaf;
};

struct Hit {
    int face_id;
    double barycentric_coords[3];
};

class Bvh
{
    private:
        std::vector<std::vector<double>> centroids;
        std::vector<double> sorted_triangle_ids;
        std::vector<Node> nodes;

        std::vector<std::vector<int>> &triangles;
        std::vector<std::vector<double>> &points;

        void buildNode(int nodeId, int start, int end, int level, int targetCellSize, int maxDepth);
        int split(int start, int end, int dim, float split_value);
        void intersectNode(int nodeId, std::vector<double> &point, Hit &hit, bool &found);

    public:
        Bvh(std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &points);
        ~Bvh();

        void build(int targetCellSize, int maxDepth);
        void query(std::vector<double> point, Hit &hit, bool &intersection_found);
};

#endif
