#ifndef _CAUSTIC_DESIGN_H
#define _CAUSTIC_DESIGN_H

#include "mesh.h"
#include "solver.h"

class Caustic_design
{
private:
    /* data */

public:
    Caustic_design(/* args */);
    ~Caustic_design();

    Mesh *mesh;
    std::vector<std::vector<double>> phi;
    std::vector<double> errors;
    std::vector<std::vector<std::vector<double>>> target_cells;
    std::vector<std::vector<std::vector<double>>> source_cells;
    std::vector<double> target_areas;
    std::vector<std::vector<double>> pixels;
    std::vector<std::vector<double>> raster;
    std::vector<std::vector<std::vector<double>>> gradient;
    std::vector<std::vector<double>> h;
    std::vector<std::vector<double>> divergence;
    std::vector<std::vector<double>> norm_x;
    std::vector<std::vector<double>> norm_y;
    std::vector<std::vector<double>> vertex_gradient;
    std::vector<std::vector<double>> normals;

    int mesh_res_x;
    int mesh_res_y;

    int resolution_x;
    int resolution_y;

    double width;
    double height;

    double focal_l;
    double thickness;
    int nthreads;

    std::vector<double> calculate_vertex_normal(std::vector<std::vector<double>> &points, int vertex_index);

    double perform_transport_iteration();

    void perform_height_map_iteration(int itr);

    void initialize_solvers(std::vector<std::vector<double>> image);

    void set_mesh_resolution(int width, int heigth);
    void set_domain_resolution(int width, int heigth);
    void set_mesh_size(double width, double heigth);
    void set_lens_focal_length(double focal_length);
    void set_lens_thickness(double thickness);
    void set_solver_max_threads(int n_threads);

    void save_solid_obj_target(const std::string& filename);
    void save_solid_obj_source(const std::string& filename);

    void export_paramererization_to_svg(const std::string& filename, double line_width);

    void export_inverted_transport_map(std::string filename, double stroke_width);
};

#endif // _CAUSTIC_DESIGN_H
