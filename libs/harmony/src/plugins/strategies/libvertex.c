/*
 * Copyright 2003-2016 Jeffrey K. Hollingsworth
 *
 * This file is part of Active Harmony.
 *
 * Active Harmony is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Active Harmony is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
 */
#define _XOPEN_SOURCE 500 // Needed for M_PI.

#include "libvertex.h"
#include "session-core.h"
#include "hspace.h"
#include "hpoint.h"
#include "hperf.h"

#include <stdlib.h> // For free(), realloc(), and NULL.
#include <string.h> // For memcpy() and memset().
#include <math.h>   // For fabs(), nextafter(), HUGE_VAL, and NAN.

/*
 * Internal helper function prototypes.
 */
static int       copy_vertex(vertex_t* dest, const vertex_t* src,
                             unsigned count);
static double    fraction(const hval_t* val, const hrange_t* range);
static int       inbounds(const vertex_t* base, unsigned count,
                          const hspace_t* space);
static double    l1_norm(const vertex_t* a, const vertex_t* b);
static double    l2_norm(const vertex_t* a, const vertex_t* b);
static void      pair_by_index(unsigned n, unsigned index, unsigned pair[]);
static int       rotate_simplex(simplex_t* simplex, int dimensions);
static unsigned* shuffle(unsigned size);
static int       transform(const vertex_t* base, unsigned count,
                           const vertex_t* origin, double coefficient,
                           vertex_t* result);
static void      unit_simplex(simplex_t* simplex, unsigned dimensions);
static int       validate_simplex(simplex_t* simplex, const hspace_t* space);
static hval_t    value(double fraction, const hrange_t* range);

/*
 * Vertex structure management implementation.
 */
int vertex_init(vertex_t* vertex, int newlen)
{
    if (vertex->len != newlen) {
        double* newbuf = realloc(vertex->term, newlen * sizeof(*newbuf));
        if (!newbuf)
            return -1;

        vertex->term = newbuf;
        vertex->len  = newlen;
    }
    return 0;
}

int vertex_copy(vertex_t* dest, const vertex_t* src)
{
    return copy_vertex(dest, src, 1);
}

void vertex_fini(vertex_t* vertex)
{
    hperf_fini(&vertex->perf);
    free(vertex->term);
}

/*
 * Vertex initialization implementation.
 */
int vertex_center(vertex_t* vertex, const hspace_t* space)
{
    if (vertex->len != space->len) {
        if (vertex_init(vertex, space->len) != 0)
            return -1;
    }

    for (int i = 0; i < space->len; ++i) {
        if (hrange_finite(&space->dim[i])) {
            vertex->term[i] = 0.5;
        }
        else {
            vertex->term[i]  = space->dim[i].bounds.r.max;
            vertex->term[i] -= space->dim[i].bounds.r.min;
            vertex->term[i] /= 2;
            vertex->term[i] += space->dim[i].bounds.r.min;
        }
    }

    hperf_reset(&vertex->perf);
    return 0;
}

int vertex_maximum(vertex_t* vertex, const hspace_t* space)
{
    if (vertex->len != space->len) {
        if (vertex_init(vertex, space->len) != 0)
            return -1;
    }

    for (int i = 0; i < space->len; ++i) {
        if (hrange_finite(&space->dim[i]))
            vertex->term[i] = nextafter(1.0, 0.0);
        else
            vertex->term[i] = space->dim[i].bounds.r.max;
    }

    hperf_reset(&vertex->perf);
    return 0;
}

int vertex_minimum(vertex_t* vertex, const hspace_t* space)
{
    if (vertex->len != space->len) {
        if (vertex_init(vertex, space->len) != 0)
            return -1;
    }

    for (int i = 0; i < space->len; ++i) {
        if (hrange_finite(&space->dim[i]))
            vertex->term[i] = 0.0;
        else
            vertex->term[i] = space->dim[i].bounds.r.min;
    }

    hperf_reset(&vertex->perf);
    return 0;
}

int vertex_parse(vertex_t* vertex, const hspace_t* space, const char* buf)
{
    hpoint_t point = hpoint_zero;
    int retval = 0;

    if (hpoint_parse(&point, buf, space) != 0)
        goto error;

    if (hpoint_align(&point, space) != 0)
        goto error;

    if (vertex_set(vertex, space, &point) != 0)
        goto error;

    goto cleanup;

  error:
    retval = -1;

  cleanup:
    hpoint_fini(&point);
    return retval;
}

int vertex_random(vertex_t* vertex, const hspace_t* space, double radius)
{
    if (vertex->len != space->len) {
        if (vertex_init(vertex, space->len) != 0)
            return -1;
    }

    for (int i = 0; i < space->len; ++i) {
        double base, range;
        if (hrange_finite(&space->dim[i])) {
            base  = 0.0;
            range = 1.0;
        }
        else {
            range_real_t* bounds = &space->dim[i].bounds.r;
            base  = bounds->min;
            range = nextafter(bounds->max - bounds->min, HUGE_VAL);
        }

        double fraction = search_drand48();
        vertex->term[i]  = range * (      radius) * fraction; // Target window.
        vertex->term[i] += range * (1.0 - radius) / 2.0;     // Excluded frame.
        vertex->term[i] += base;
    }

    hperf_reset(&vertex->perf);
    return 0;
}

int vertex_set(vertex_t* vertex, const hspace_t* space, const hpoint_t* point)
{
    if (vertex->len != space->len) {
        if (vertex_init(vertex, space->len) != 0)
            return -1;
    }

    for (int i = 0; i < space->len; ++i) {
        if (hrange_finite(&space->dim[i]))
            vertex->term[i] = fraction(&point->term[i], &space->dim[i]);
        else
            vertex->term[i] = point->term[i].value.r;
    }

    vertex->id = point->id;
    hperf_reset(&vertex->perf);
    return 0;
}

/*
 * Vertex utility function implementation.
 */
int vertex_inbounds(const vertex_t* vertex, const hspace_t* space)
{
    return inbounds(vertex, 1, space);
}

double vertex_norm(const vertex_t* a, const vertex_t* b, vertex_norm_t norm)
{
    switch (norm) {
    case VERTEX_NORM_L1: return l1_norm(a, b);
    case VERTEX_NORM_L2: return l2_norm(a, b);
    default:             return NAN;
    }
}

int vertex_point(const vertex_t* vertex, const hspace_t* space,
                 hpoint_t* point)
{
    if (point->cap < space->len) {
        if (hpoint_init(point, space->len) != 0)
            return -1;
    }

    for (int i = 0; i < space->len; ++i) {
        if (hrange_finite(&space->dim[i])) {
            point->term[i] = value(vertex->term[i], &space->dim[i]);
        }
        else {
            point->term[i].type    = HVAL_REAL;
            point->term[i].value.r = vertex->term[i];
        }
    }

    point->len = space->len;
    point->id  = vertex->id;
    return 0;
}

int vertex_transform(const vertex_t* base, const vertex_t* origin,
                     double coefficient, vertex_t* result)
{
    return transform(base, 1, origin, coefficient, result);
}

/*
 * Simplex structure management implementation.
 */
int simplex_init(simplex_t* simplex, unsigned dimensions)
{
    unsigned newlen = dimensions + 1;
    if (simplex->len != newlen) {
        vertex_t* newbuf = realloc(simplex->vertex, newlen * sizeof(*newbuf));
        if (!newbuf)
            return -1;

        unsigned oldlen = simplex->len;
        memset(&newbuf[ oldlen ], 0, (newlen - oldlen) * sizeof(*newbuf));
        simplex->vertex = newbuf;
        simplex->len    = newlen;
    }

    for (int i = 0; i < newlen; ++i) {
        if (simplex->vertex[i].len != dimensions) {
            if (vertex_init(&simplex->vertex[i], dimensions) != 0)
                return -1;
        }
    }
    return 0;
}

int simplex_copy(simplex_t* dest, const simplex_t* src)
{
    if (dest->len != src->len) {
        if (simplex_init(dest, src->len - 1) != 0)
            return -1;
    }
    return copy_vertex(dest->vertex, src->vertex, src->len);
}

void simplex_fini(simplex_t* simplex)
{
    for (unsigned i = 0; i < simplex->len; ++i)
        vertex_fini(&simplex->vertex[i]);
    free(simplex->vertex);
}

/*
 * Simplex initialization implementation.
 */
int simplex_set(simplex_t* simplex, const hspace_t* space,
                const vertex_t* base, double radius)
{
    int dimensions = space->len;

    if (simplex->len != dimensions + 1) {
        if (simplex_init(simplex, dimensions) != 0)
            return -1;
    }

    // Generate a unit simplex.
    unit_simplex(simplex, dimensions);

    // Rotate it randomly.
    rotate_simplex(simplex, dimensions);

    radius /= 2.0;
    for (unsigned j = 0; j < dimensions; ++j) {
        double range;
        if (hrange_finite(&space->dim[j]))
            range = 1.0;
        else
            range = space->dim[j].bounds.r.max - space->dim[j].bounds.r.min;

        for (unsigned i = 0; i <= dimensions; ++i) {
            simplex->vertex[i].term[j] *= range * radius;
            simplex->vertex[i].term[j] += base->term[j];
        }
    }

    // Move the simplex into the space if necessary.
    if (validate_simplex(simplex, space) != 0)
        return -1;

    return 0;
}

/*
 * Simplex utility function implementation.
 */

/*
 * General centroid formula for vertex points x1 through xN:
 *     sum(x1, x2, x3, ... , xN) / N
 *
 * This function will ignore any vertices with id == 0.
 */
int simplex_centroid(const simplex_t* simplex, vertex_t* centroid)
{
    int total = 0;
    int vert_len = simplex->vertex[0].len;
    int perf_len = simplex->vertex[0].perf.len;

    if (centroid->len != vert_len) {
        if (vertex_init(centroid, vert_len) != 0)
            return -1;
    }

    if (centroid->perf.len != perf_len) {
        if (hperf_init(&centroid->perf, perf_len) != 0)
            return -1;
    }

    memset(centroid->term,     0, vert_len * sizeof(*centroid->term));
    memset(centroid->perf.obj, 0, perf_len * sizeof(*centroid->perf.obj));

    for (unsigned i = 0; i < simplex->len; ++i) {
        if (!simplex->vertex[i].id)
            continue;

        for (int j = 0; j < vert_len; ++j)
            centroid->term[j] += simplex->vertex[i].term[j];

        for (int j = 0; j < perf_len; ++j)
            centroid->perf.obj[j] += simplex->vertex[i].perf.obj[j];

        ++total;
    }

    for (int j = 0; j < vert_len; ++j)
        centroid->term[j] /= total;

    for (int j = 0; j < perf_len; ++j)
        centroid->perf.obj[j] /= total;

    return 0;
}

int simplex_collapsed(const simplex_t* simplex, const hspace_t* space)
{
    for (unsigned j = 0; j < space->len; ++j) {
        hval_t value_0, value_i;

        value_0 = value(simplex->vertex[0].term[j], &space->dim[j]);
        for (int i = 1; i < simplex->len; ++i) {
            value_i = value(simplex->vertex[i].term[j], &space->dim[j]);

            if (!hval_eq(&value_0, &value_i))
                return 0;
        }
    }
    return -1;
}

int simplex_inbounds(const simplex_t* simplex, const hspace_t* space)
{
    return inbounds(simplex->vertex, simplex->len, space);
}

int simplex_transform(const simplex_t* base, const vertex_t* origin,
                      double coefficient, simplex_t* result)
{
    if (result->len != base->len) {
        if (simplex_init(result, base->len - 1) != 0)
            return -1;
    }

    return transform(base->vertex, base->len, origin, coefficient,
                     result->vertex);
}

/*
 * Internal helper function implementation.
 */
int copy_vertex(vertex_t* dest, const vertex_t* src, unsigned count)
{
    for (unsigned i = 0; i < count; ++i) {
        if (dest[i].len != src[i].len) {
            if (vertex_init(&dest[i], src[i].len) != 0)
                return -1;
        }

        if (hperf_copy(&dest[i].perf, &src[i].perf) != 0) {
            return -1;
        }

        memcpy(dest[i].term, src[i].term, src[i].len * sizeof(*src[i].term));
        dest[i].id = src[i].id;
    }
    return 0;
}

/*
 * Convert a raw value into a number between [0, 1).
 */
double fraction(const hval_t* val, const hrange_t* range)
{
    double fraction;

    fraction  = 0.5;
    fraction += hrange_index(range, val);
    fraction /= hrange_limit(range);

    return fraction;
}

/*
 * Check if a vertex falls within the bounds of a search space.
 */
int inbounds(const vertex_t* vertex, unsigned count, const hspace_t* space)
{
    if (!vertex)
        return 0;

    for (int j = 0; j < space->len; ++j) {
        double min, max;
        if (hrange_finite(&space->dim[j])) {
            min = 0.0;
            max = nextafter(1.0, 0.0);
        }
        else {
            min = space->dim[j].bounds.r.min;
            max = space->dim[j].bounds.r.max;
        }

        for (unsigned i = 0; i < count; ++i) {
            if (vertex[i].term[j] < min || max < vertex[i].term[j])
                return 0;
        }
    }
    return 1;
}

double l1_norm(const vertex_t* a, const vertex_t* b)
{
    double sum = 0.0;
    for (int i = 0; i < a->len; ++i)
        sum += fabs(a->term[i] - b->term[i]);

    return sum;
}

double l2_norm(const vertex_t* a, const vertex_t* b)
{
    double sum = 0.0;
    for (int i = 0; i < a->len; ++i) {
        double offset = a->term[i] - b->term[i];
        sum += offset * offset;
    }
    return sqrt(sum);
}

/*
 * Given an index (i), returns the i'th lexicographical combination of
 * choosing two elements from a set of size N.
 */
void pair_by_index(unsigned n, unsigned index, unsigned pair[])
{
    int limit = 0;
    unsigned first = 0;

    ++index;
    do limit += n - ++first;
    while (limit < index);

    pair[0] = first - 1;
    pair[1] = n + index - limit - 1;
}

/*
 * Randomly rotate a simplex in N dimensions about the origin.
 */
int rotate_simplex(simplex_t* simplex, int dimensions)
{
    unsigned combos = (dimensions * (dimensions - 1)) / 2;

    // Produce a random order of (N choose 2) integers.
    unsigned* order = shuffle(combos);
    if (!order)
        return -1;

    for (unsigned i = 0; i < combos; ++i) {
        // Produce a random angle for each pair of dimensions.
        double theta = search_drand48() * (2 * M_PI);

        for (unsigned j = 0; j <= dimensions; ++j) {
            unsigned pair[2];
            pair_by_index(dimensions, order[i], pair);

            double first  = simplex->vertex[j].term[ pair[0] ];
            double second = simplex->vertex[j].term[ pair[1] ];
            simplex->vertex[j].term[ pair[0] ] = (first  * cos(theta) -
                                                  second * sin(theta));
            simplex->vertex[j].term[ pair[1] ] = (first  * sin(theta) +
                                                  second * cos(theta));
        }
    }

    free(order);
    return 0;
}

/*
 * Generate a random permutation of N indexes using the Knuth shuffle.
 * Returns a newly allocated heap array which must be freed.
 */
unsigned* shuffle(unsigned size)
{
    unsigned* retval = malloc(size * sizeof(*retval));
    if (!retval)
        return NULL;

    for (unsigned i = 0; i < size; ++i) {
        unsigned j = search_lrand48() % (i + 1);
        retval[i] = retval[j];
        retval[j] = i;
    }
    return retval;
}

/*
 * Translate each base point by a coefficient of its distance from an
 * origin point.  For example, a coefficient of 0.0 produces the base
 * point, and a coefficient of -1.0 produces the origin.
 */
int transform(const vertex_t* base, unsigned count, const vertex_t* origin,
              double coefficient, vertex_t* result)
{
    int dimensions = origin->len;

    for (unsigned i = 0; i < count; ++i) {
        if (base[i].len != dimensions)
            return -1;

        if (result[i].len != dimensions) {
            if (vertex_init(&result[i], dimensions) != 0)
                return -1;
        }

        for (int j = 0; j < dimensions; ++j) {
            double offset = (base[i].term[j] - origin->term[j]) * coefficient;
            result[i].term[j] = base[i].term[j] + offset;
        }
    }
    return 0;
}

/*
 * Write a regular N-dimensional simplex (whose vertex-to-origin
 * distance is 1) into the given simplex_t structure.
 */
void unit_simplex(simplex_t* simplex, unsigned dimensions)
{
    simplex->vertex[0].term[0] = 1.0;
    for (unsigned i = 0; i < dimensions; ++i) {
        for (unsigned j = 0; j < i; ++j)
            simplex->vertex[j].term[i] = 0.0;

        if (i > 0) {
            double a = simplex->vertex[i-1].term[i-1];
            double b = simplex->vertex[i  ].term[i-1];
            simplex->vertex[i].term[i] = sqrt(a * a - b * b);
        }

        double portion = simplex->vertex[i].term[i] / (dimensions - i);
        for (unsigned j = i + 1; j <= dimensions; ++j)
            simplex->vertex[j].term[i] = -portion;
    }
}

/*
 * Check that a simplex is within the bounds of its search space,
 * shifting the entire simplex as necessary.  If shifting the simplex
 * will not resolve the issue (i.e., the simplex is larger than the
 * search space), return -1.
 *
 * This function only considers the first D+1 vertices, where D is the
 * number of search space dimensions.
 */
int validate_simplex(simplex_t* simplex, const hspace_t* space)
{
    // Process each dimension in isolation.
    for (unsigned i = 0; i < space->len; ++i) {
        // Find the lowest and highest value for the i'th term.
        double low = simplex->vertex[0].term[i];
        double high = simplex->vertex[0].term[i];
        for (unsigned j = 1; j <= space->len; ++j) {
            if (low > simplex->vertex[j].term[i])
                low = simplex->vertex[j].term[i];
            if (high < simplex->vertex[j].term[i])
                high = simplex->vertex[j].term[i];
        }

        // Verify that simplex is not larger than the search space.
        double min;
        double max;
        if (hrange_finite(&space->dim[i])) {
            min = 0.0;
            max = nextafter(1.0, 0.0);
        }
        else {
            min = space->dim[i].bounds.r.min;
            max = space->dim[i].bounds.r.max;
        }

        // Fail if simplex is larger than the search space.
        if ((high - low) >= (max - min))
            return -1;

        // Shift all vertices up, if necessary.
        if (low < min) {
            double shift = min - low;
            for (unsigned j = 0; j < simplex->len; ++j)
                simplex->vertex[j].term[i] += shift;
        }

        // Shift all vertices down, if necessary.
        if (high >= max) {
            double shift = high - max;
            for (unsigned j = 0; j < simplex->len; ++j)
                simplex->vertex[j].term[i] -= shift;
        }
    }
    return 0;
}

/*
 * Convert a fractional value [0, 1) to a raw value.
 */
hval_t value(double fraction, const hrange_t* range)
{
    unsigned long index;

    fraction *= hrange_limit(range);
    index = (unsigned long) fraction;

    return hrange_value(range, index);
}
