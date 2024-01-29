// Source for the polygon clipping algorithm
// https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <algorithm>

#include "polygon_utils.h"

inline double dot(vec_t *a, vec_t *b)
{
	return a->x * b->x + a->y * b->y;
}

inline double cross(vec_t *a, vec_t *b)
{
	return a->x * b->y - a->y * b->x;
}

inline vec_t *vsub(vec_t *a, vec_t *b, vec_t *res)
{
	res->x = a->x - b->x;
	res->y = a->y - b->y;
	return res;
}

/* tells if vec c lies on the left side of directed edge a->b
 * 1 if left, -1 if right, 0 if colinear
 */
int left_of(vec_t *a, vec_t *b, vec_t *c)
{
	vec_t tmp1, tmp2;
	double x;
	vsub(b, a, &tmp1);
	vsub(c, b, &tmp2);
	x = cross(&tmp1, &tmp2);
	return x < 0 ? -1 : x > 0;
}

int line_sect(vec_t *x0, vec_t *x1, vec_t *y0, vec_t *y1, vec_t *res)
{
	vec_t dx, dy, d;
	vsub(x1, x0, &dx);
	vsub(y1, y0, &dy);
	vsub(x0, y0, &d);
	/* x0 + a dx = y0 + b dy ->
	   x0 X dx = y0 X dx + b dy X dx ->
	   b = (x0 - y0) X dx / (dy X dx) */
	double dyx = cross(&dy, &dx);
	if (!dyx) return 0;
	dyx = cross(&d, &dx) / dyx;
	if (dyx <= 0 || dyx >= 1) return 0;

	res->x = y0->x + dyx * dy.x;
	res->y = y0->y + dyx * dy.y;
	return 1;
}

poly_t *poly_new()
{
	return (poly_t *)calloc(1, sizeof(poly_t));
}

void poly_free(poly_t *p)
{
	free(p->v);
	free(p);
}

void poly_append(poly_t *p, vec_t *v)
{
	if (p->len >= p->alloc) {
		p->alloc *= 2;
		if (!p->alloc) p->alloc = 4;
		p->v = (vec_t *)realloc(p->v, sizeof(vec_t) * p->alloc);
	}
	p->v[p->len++] = *v;
}

/* this works only if all of the following are true:
 *   1. poly has no colinear edges;
 *   2. poly has no duplicate vertices;
 *   3. poly has at least three vertices;
 *   4. poly is convex (implying 3).
*/
int poly_winding(const poly_t *p)
{
	return left_of(p->v, p->v + 1, p->v + 2);
}

void poly_edge_clip(const poly_t *sub, vec_t *x0, vec_t *x1, int left, poly_t *res)
{
	int i, side0, side1;
	vec_t tmp;
	vec_t *v0 = sub->v + sub->len - 1, *v1;
	res->len = 0;

	side0 = left_of(x0, x1, v0);
	if (side0 != -left) poly_append(res, v0);

	for (i = 0; i < sub->len; i++) {
		v1 = sub->v + i;
		side1 = left_of(x0, x1, v1);
		if (side0 + side1 == 0 && side0)
			/* last point and current straddle the edge */
			if (line_sect(x0, x1, v0, v1, &tmp))
				poly_append(res, &tmp);
		if (i == sub->len - 1) break;
		if (side1 != -left) poly_append(res, v1);
		v0 = v1;
		side0 = side1;
	}
}

poly_t *poly_clip(const poly_t *sub, const poly_t *clip)
{
	int i;
	poly_t *p1 = poly_new(), *p2 = poly_new(), *tmp;

	int dir = poly_winding(clip);
	poly_edge_clip(sub, clip->v + clip->len - 1, clip->v, dir, p2);
	for (i = 0; i < clip->len - 1; i++) {
		tmp = p2; p2 = p1; p1 = tmp;
		if(p1->len == 0) {
			p2->len = 0;
			break;
		}
		poly_edge_clip(p1, clip->v + i, clip->v + i + 1, dir, p2);
	}

	poly_free(p1);
	return p2;
}

bbox_t calculateBoundingBox(const poly_t *polygon) {
    bbox_t boundingBox;

    if (polygon->len == 0) {
        // Handle empty polygon
        boundingBox.xmin = boundingBox.xmax = boundingBox.ymin = boundingBox.ymax = 0.0;
        return boundingBox;
    }

    // Initialize bounding box with the first vertex
    boundingBox.xmin = boundingBox.xmax = polygon->v[0].x;
    boundingBox.ymin = boundingBox.ymax = polygon->v[0].y;

    // Iterate through the remaining vertices to update the bounding box
    for (int i = 1; i < polygon->len; ++i) {
        double x = polygon->v[i].x;
        double y = polygon->v[i].y;

        if (x < boundingBox.xmin) {
            boundingBox.xmin = x;
        } else if (x > boundingBox.xmax) {
            boundingBox.xmax = x;
        }

        if (y < boundingBox.ymin) {
            boundingBox.ymin = y;
        } else if (y > boundingBox.ymax) {
            boundingBox.ymax = y;
        }
    }

    return boundingBox;
}

double calculate_polygon_area(const poly_t *polygon) {
    if (polygon == NULL || polygon->len < 3) {
        return 0.0;
    }

	int n = polygon->len;
    double area = 0.0;

    for (int i = 0; i < n; i++) {
		double x1 = polygon->v[i].x;
		double y1 = polygon->v[i].y;
		double x2 = polygon->v[(i + 1) % n].x;
		double y2 = polygon->v[(i + 1) % n].y;
        area += (x1 * y2) - (x2 * y1);
    }

    // Take the absolute value and divide by 2
    area = 0.5 * (area);

    return area;
}

double integrate_intensity(poly_t **polygons, double *intensities, int num_polygons) {
    double intensity = 0.0;

    for (int i = 0; i < num_polygons; i++) {
        double area = calculate_polygon_area(polygons[i]);
        intensity += area * intensities[i];
    }

    return intensity;
}

void calculate_bounding_box(poly_t *polygon, double *min_x, double *min_y, double *max_x, double *max_y) {
    // Extract x and y coordinates
    double *x_coords = (double *)malloc(polygon->len * sizeof(double));
    double *y_coords = (double *)malloc(polygon->len * sizeof(double));
    
    for (int i = 0; i < polygon->len; i++) {
        x_coords[i] = polygon->v[i].x;
        y_coords[i] = polygon->v[i].y;
    }

    // Calculate bounding box coordinates
    *min_x = x_coords[0];
    *min_y = y_coords[0];
    *max_x = x_coords[0];
    *max_y = y_coords[0];

    for (int i = 1; i < polygon->len; i++) {
        if (x_coords[i] < *min_x) {
            *min_x = x_coords[i];
        }
        if (y_coords[i] < *min_y) {
            *min_y = y_coords[i];
        }
        if (x_coords[i] > *max_x) {
            *max_x = x_coords[i];
        }
        if (y_coords[i] > *max_y) {
            *max_y = y_coords[i];
        }
    }

    free(x_coords);
    free(y_coords);
}

void reverseVertices(poly_t *polygon) {
    if (polygon->len <= 1) {
        // Nothing to reverse for polygons with 0 or 1 vertices
        return;
    }

    // Allocate temporary memory to store reversed vertices
    vec_t *reversed = (vec_t *)malloc(polygon->len * sizeof(vec_t));

    if (reversed == NULL) {
        // Handle memory allocation failure
        // (You might want to add error handling based on your specific requirements)
        return;
    }

    // Copy vertices in reverse order
    for (int i = 0; i < polygon->len; i++) {
        reversed[i] = polygon->v[polygon->len - 1 - i];
    }

    // Replace the original vertices with the reversed ones
    free(polygon->v); // Free the original memory
    polygon->v = reversed;
}

double integrate_cell_intensities(double *image, poly_t *polygon, int image_w, int image_h, double width) {
	bbox_t bbox = calculateBoundingBox(polygon);

	vec_t square_vertices[4];

	double intensity = 0.0;

	double px_side_length = width / ((double)image_w);

	for (int y = std::fmax(floor(bbox.ymin/px_side_length), 0.0f); y < std::fmin(ceil(bbox.ymax/px_side_length), image_h); ++y) {
		for (int x = std::fmax(floor(bbox.xmin/px_side_length), 0.0f); x < std::fmin(ceil(bbox.xmax/px_side_length), image_w); x++) {
            double center[2] = {(double)x + 0.5, (double)y + 0.5};

			center[0] *= px_side_length;
			center[1] *= px_side_length;

			square_vertices[0].x = center[0] - px_side_length / 2.0f;
			square_vertices[0].y = center[1] - px_side_length / 2.0f;

			square_vertices[1].x = center[0] - px_side_length / 2.0f;
			square_vertices[1].y = center[1] + px_side_length / 2.0f;

			square_vertices[2].x = center[0] + px_side_length / 2.0f;
			square_vertices[2].y = center[1] + px_side_length / 2.0f;

			square_vertices[3].x = center[0] + px_side_length / 2.0f;
			square_vertices[3].y = center[1] - px_side_length / 2.0f;

			poly_t *pixel_poly = poly_new();
			poly_append(pixel_poly, &(square_vertices[0]));
			poly_append(pixel_poly, &(square_vertices[1]));
			poly_append(pixel_poly, &(square_vertices[2]));
			poly_append(pixel_poly, &(square_vertices[3]));

			poly_t *result = poly_clip(polygon, pixel_poly);

			intensity += calculate_polygon_area(result) * image[(image_h - y - 1) * image_w + x];

			poly_free(pixel_poly);
			poly_free(result);
        }
    }

    return intensity;
}

void get_target_areas(double *image, poly_t **polygons, int num_polygons, int image_w, int image_h, double width, double height, double *target_areas) {
	double sum_target_area = 0.0f;

	for (int i=0; i<num_polygons; i++) {
		target_areas[i] = integrate_cell_intensities(image, polygons[i], image_w, image_h, width);
		sum_target_area += target_areas[i];
	}

	double scaling_factor = (width * height) / sum_target_area;

	for (int i=0; i<num_polygons; i++) {
		target_areas[i] *= scaling_factor;
	}
}

void integrate_cell_gradient(double *grad_x, double *grad_y, poly_t *polygon, int grad_w, int grad_h, double width, double *interp_x, double *interp_y) {
	bbox_t bbox = calculateBoundingBox(polygon);

	vec_t square_vertices[4];

	double intensity = 0.0;

	for (int y = fmax(floor((bbox.ymin) * (((double)grad_w)/width)), 0); y < fmin(ceil((bbox.ymax) * (((double)grad_w)/width)), grad_h); ++y) {
		for (int x = fmax(floor((bbox.xmin) * (((double)grad_w)/width)), 0); x < fmin(ceil((bbox.xmax) * (((double)grad_w)/width)), grad_w); x++) {
            double center[2] = {(double)x + 0.5, (double)y + 0.5};
            double side_length = width / ((double)grad_w);

			center[0] *= side_length;
			center[1] *= side_length;

			square_vertices[0].x = center[0] + side_length / 2.0f;
			square_vertices[0].y = center[1] - side_length / 2.0f;

			square_vertices[1].x = center[0] + side_length / 2.0f;
			square_vertices[1].y = center[1] + side_length / 2.0f;

			square_vertices[2].x = center[0] - side_length / 2.0f;
			square_vertices[2].y = center[1] + side_length / 2.0f;

			square_vertices[3].x = center[0] - side_length / 2.0f;
			square_vertices[3].y = center[1] - side_length / 2.0f;

			poly_t *pixel_poly = poly_new();
			poly_append(pixel_poly, &(square_vertices[0]));
			poly_append(pixel_poly, &(square_vertices[1]));
			poly_append(pixel_poly, &(square_vertices[2]));
			poly_append(pixel_poly, &(square_vertices[3]));

			poly_t *result = poly_clip(polygon, pixel_poly);

			double area = calculate_polygon_area(result);
			*interp_x += area * grad_x[(grad_h - y - 1) * grad_w + x];
			*interp_y += area * grad_y[(grad_h - y - 1) * grad_w + x];

			poly_free(pixel_poly);
			poly_free(result);
        }
    }
}