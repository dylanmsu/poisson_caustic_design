#include "bvh.h"

// the implementation was heavily based on the work of Gael Guennebaud at https://github.com/ggael/otmap

// Copyright (C) 2016-2018 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

Bvh::Bvh(std::vector<std::vector<int>> &triangles, std::vector<std::vector<double>> &points)
    :triangles(triangles), points(points)
{
}

Bvh::~Bvh()
{
}

void Bvh::build(int targetCellSize, int maxDepth) {
    // Reset tree data
    centroids.clear();
    sorted_triangle_ids.clear();
    nodes.clear();

    // Sdd the first empty node
    nodes.resize(1);
    
    // Set size beforehand
    centroids.resize(triangles.size());
    sorted_triangle_ids.resize(triangles.size());

    // Initialize the tree data with centroids and triangle indices
    for (int i=0; i<triangles.size(); i++) {
        std::vector<std::vector<double>> triangle;
        triangle.resize(triangles[i].size());

        for (int j=0; j<triangle.size(); j++) {
            triangle[j] = points[triangles[i][j]];
        }

        centroids[i] = calculate_polygon_centroid(triangle);
        sorted_triangle_ids[i] = i;
    }

    // Build the BHV-tree recursively 
    buildNode(0, 0, triangles.size(), 0, targetCellSize, maxDepth);
}

int Bvh::split(int start, int end, int dim, float split_value)
{
    int left = start, right = end - 1;

    while (left < right)
    {
        // Find the first on the left
        while (left < end && centroids[left][dim] < split_value)
            left += 1;

        // Find the first on the right
        while (right >= start && centroids[right][dim] >= split_value)
            right -= 1;

        if (left >= right)
            break;

        // Swap centroids and corresponding triangle IDs
        std::swap(centroids[left], centroids[right]);
        std::swap(sorted_triangle_ids[left], sorted_triangle_ids[right]);

        ++left;
        --right;
    }

    return centroids[left][dim] <= split_value ? std::min(end, left + 1) : left;
}

void Bvh::buildNode(int nodeId, int start, int end, int level, int targetCellSize, int maxDepth)
{
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (int i = start; i < end; i++) {
        for (int j = 0; j < triangles[sorted_triangle_ids[i]].size(); j++) {
            int point_index = triangles[sorted_triangle_ids[i]][j];

            double x = points[point_index][0];
            double y = points[point_index][1];

            min_x = std::fmin(min_x, x);
            min_y = std::fmin(min_y, y);
            max_x = std::fmax(max_x, x);
            max_y = std::fmax(max_y, y);
        }
    }

    double epsilon = std::numeric_limits<double>::epsilon();

    // Enlarge the bounding box by epsilon
    min_x -= 0.5 * epsilon;
    min_y -= 0.5 * epsilon;
    max_x += 0.5 * epsilon;
    max_y += 0.5 * epsilon;

    // Now, min_x, min_y, max_x, and max_y represent the enlarged bounding box
    nodes[nodeId].bbox_min_x = min_x;
    nodes[nodeId].bbox_min_y = min_y;
    nodes[nodeId].bbox_max_x = max_x;
    nodes[nodeId].bbox_max_y = max_y;

    // stopping criteria
    if(end-start <= targetCellSize || level>=maxDepth)
    {
        // we got a leaf !
        nodes[nodeId].is_leaf = true;
        nodes[nodeId].first_face_id = start;
        nodes[nodeId].nb_faces = std::max(0,end-start);
        return;
    }
    nodes[nodeId].is_leaf = false;

    int dim = 0;
    if ((max_x - min_x) < (max_y - min_y)) {
        dim = 1;
    }

    // Split at the middle
    double split_value = 0;
    if (dim == 0) {
        split_value = 0.5f * (nodes[nodeId].bbox_max_x + nodes[nodeId].bbox_min_x);
    } else {
        split_value = 0.5f * (nodes[nodeId].bbox_max_y + nodes[nodeId].bbox_min_y);
    }

    //std::cout << "level=" << level << " dim=" << dim << ", val=" << split_value << std::endl;

    // Sort the faces according to the split plane
    int mid_id = split(start, end, dim, split_value);

    // second stopping criteria
    if(mid_id==start || mid_id==end)
    {
        // no improvement
        nodes[nodeId].is_leaf = true;
        nodes[nodeId].first_face_id = start;
        nodes[nodeId].nb_faces = std::max(0,end-start);
        return;
    }

    // create the children
    int child_id = nodes.size();
    nodes[nodeId].first_child_id = nodes.size();
    nodes.resize(nodes.size()+2);
    // node is not a valid reference anymore !

    buildNode(child_id  , start, mid_id, level+1, targetCellSize, maxDepth);
    buildNode(child_id+1, mid_id, end, level+1, targetCellSize, maxDepth);
}

std::vector<double> get_barycentric_coordinates(const std::vector<double>& t0, const std::vector<double>& t1, const std::vector<double>& t2, const std::vector<double>& point) {
    // Calculate the vectors of the triangle sides
    std::vector<double> v0(2), v1(2), v2(2);
    v0[0] = t2[0] - t0[0];
    v0[1] = t2[1] - t0[1];
    v1[0] = t1[0] - t0[0];
    v1[1] = t1[1] - t0[1];
    v2[0] = point[0] - t0[0];
    v2[1] = point[1] - t0[1];

    // Calculate dot products
    double dot00 = v0[0] * v0[0] + v0[1] * v0[1];
    double dot01 = v0[0] * v1[0] + v0[1] * v1[1];
    double dot02 = v0[0] * v2[0] + v0[1] * v2[1];
    double dot11 = v1[0] * v1[0] + v1[1] * v1[1];
    double dot12 = v1[0] * v2[0] + v1[1] * v2[1];

    // Check for division by zero
    double denom = dot00 * dot11 - dot01 * dot01;

    // Calculate barycentric coordinates
    double inv_denom = 1 / denom;
    double u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
    double v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
    double w = 1.0f - u - v;

    // Return barycentric coordinates
    return {u, v, w};
}

bool is_inside_bbox(const Node &node, const std::vector<double> &point) {
    return node.bbox_min_x <= point[0] && point[0] <= node.bbox_max_x && node.bbox_min_y <= point[1] && point[1] <= node.bbox_max_y;
}

void Bvh::intersectNode(int nodeId, std::vector<double> &point, Hit &hit, bool &found)
{
    if(nodes[nodeId].is_leaf)
    {
        if (is_inside_bbox(nodes[nodeId], point)) {
            for(int i=nodes[nodeId].first_face_id; i<nodes[nodeId].first_face_id + nodes[nodeId].nb_faces; ++i) {
                std::vector<double> bary_coord = get_barycentric_coordinates(points[triangles[sorted_triangle_ids[i]][2]], points[triangles[sorted_triangle_ids[i]][1]], points[triangles[sorted_triangle_ids[i]][0]], point);
                double eps = 1e-10;

                if ((bary_coord[0] >= -eps && bary_coord[1] >= -eps) && ((bary_coord[0] + bary_coord[1]) <= 1.0f + eps)) {
                    hit.barycentric_coords[0] = bary_coord[0];
                    hit.barycentric_coords[1] = bary_coord[1];
                    hit.barycentric_coords[2] = bary_coord[2];
                    hit.face_id = sorted_triangle_ids[i];
                    found = true;
                    return;
                }
            }
        } else {
            //std::cout << "point is not inside bbox" << std::endl;
        }
    }
    else
    {
        int child_id1 = nodes[nodeId].first_child_id;
        int child_id2 = nodes[nodeId].first_child_id+1;

        if (nodes[child_id1].bbox_min_x <= point[0] && point[0] <= nodes[child_id1].bbox_max_x &&
            nodes[child_id1].bbox_min_y <= point[1] && point[1] <= nodes[child_id1].bbox_max_y)
        {
        //if(nodes_[child_id1].box.contains(target)) {
            intersectNode(child_id1, point, hit, found);
            if(found) return;
        }
        
        if (nodes[child_id2].bbox_min_x <= point[0] && point[0] <= nodes[child_id2].bbox_max_x &&
            nodes[child_id2].bbox_min_y <= point[1] && point[1] <= nodes[child_id2].bbox_max_y)
        {
        //if(nodes_[child_id2].box.contains(target)) {
            intersectNode(child_id2, point, hit, found);
            if(found) return;
        }
    }
}

void Bvh::query(std::vector<double> point, Hit &hit, bool &intersection_found)
{
    if (nodes[0].bbox_min_x <= point[0] && point[0] <= nodes[0].bbox_max_x &&
        nodes[0].bbox_min_y <= point[1] && point[1] <= nodes[0].bbox_max_y)
    {
        intersectNode(0, point, hit, intersection_found);
    } else {
        printf("outside main bbox\r\n");
    }
}