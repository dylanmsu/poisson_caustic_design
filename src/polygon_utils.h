#ifndef POLYGON_UTILS_H
#define POLYGON_UTILS_H

#include <math.h>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x, y;
} vec_t;

typedef struct {
    int len, alloc;
    vec_t *v;
} poly_t;

typedef struct {
    double xmin, xmax, ymin, ymax;
} bbox_t;

poly_t *poly_clip(const poly_t *sub, const poly_t *clip);
void calculate_bounding_box(poly_t *polygon, double *min_x, double *min_y, double *max_x, double *max_y);

double calculate_polygon_area(const poly_t *polygon);
double calculate_polygon_area_vec(const std::vector<std::vector<double>> input_polygon);
double calculate_partitioned_cell_area(const std::vector<std::vector<std::vector<double>>> input_polygons);

double integrate_intensity(poly_t **polygons, double *intensities, int num_polygons);
void append_polygon(poly_t *polygons, int *num_polygons, poly_t *polygon);
void append_intensity(double *intensities, int *num_intensities, double intensity);

std::vector<double> calculate_polygon_centroid(std::vector<std::vector<double>> vertices);

poly_t *poly_new();
void poly_free(poly_t *p);
void poly_append(poly_t *p, vec_t *v);

double integrate_cell_intensities(std::vector<std::vector<double>> &image, poly_t *polygon, int image_w, int image_h, double width);
void integrate_cell_gradient(std::vector<std::vector<double>> &grad_x, std::vector<std::vector<double>> &grad_y, std::vector<std::vector<double>> &input_polygon, int grad_w, int grad_h, double width, double &interp_x, double &interp_y);

std::vector<std::vector<double>> integrate_cell_gradients(std::vector<std::vector<std::vector<double>>> &gradient, std::vector<std::vector<std::vector<double>>> &input_polygons, int image_w, int image_h, double width, double height);

std::vector<double> integrate_grid_into_cells(std::vector<std::vector<double>> &image, std::vector<std::vector<std::vector<double>>> &input_polygons, int image_w, int image_h, double width, double height);

std::vector<double> get_target_areas(std::vector<std::vector<double>> &image, std::vector<std::vector<std::vector<double>>> &input_polygons, int image_w, int image_h, double width, double height);
std::vector<double> get_source_areas(std::vector<std::vector<std::vector<double>>> &input_polygons);

std::vector<double> get_partitioned_source_areas(std::vector<std::vector<std::vector<std::vector<double>>>> &input_polygons);

#ifdef __cplusplus
}
#endif

#endif // POLYGON_UTILS_H