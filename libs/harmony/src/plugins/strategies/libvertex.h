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

#ifndef __LIBVERTEX_H__
#define __LIBVERTEX_H__

#include "hspace.h"
#include "hpoint.h"
#include "hperf.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Vertex structure for geometric manipulation of hpoints.
 */
typedef struct vertex {
    unsigned id;
    double*  term;
    int      len;
    hperf_t  perf;
} vertex_t;

typedef enum vertex_norm {
    VERTEX_NORM_UNKNOWN = 0,
    VERTEX_NORM_L1,
    VERTEX_NORM_L2,

    VERTEX_NORM_MAX
} vertex_norm_t;

/*
 * Vertex structure management interface.
 */
int  vertex_init(vertex_t* vertex, int newlen);
int  vertex_copy(vertex_t* dest, const vertex_t* src);
void vertex_fini(vertex_t* vertex);

/*
 * Vertex initialization interface.
 */
int vertex_center(vertex_t* vertex, const hspace_t* space);
int vertex_maximum(vertex_t* vertex, const hspace_t* space);
int vertex_minimum(vertex_t* vertex, const hspace_t* space);
int vertex_parse(vertex_t* vertex, const hspace_t* space, const char* buf);
int vertex_random(vertex_t* vertex, const hspace_t* space, double radius);
int vertex_set(vertex_t* vertex, const hspace_t* space, const hpoint_t* point);

/*
 * Vertex utility function interface.
 */
int    vertex_inbounds(const vertex_t* vertex, const hspace_t* space);
double vertex_norm(const vertex_t* a, const vertex_t* b, vertex_norm_t norm);
int    vertex_point(const vertex_t* vertex, const hspace_t* space,
                    hpoint_t* point);
int    vertex_transform(const vertex_t* base, const vertex_t* origin,
                        double coefficient, vertex_t* result);

/*
 * Simplex structure to support operations on vertex groups.
 */
typedef struct simplex {
    vertex_t* vertex;
    int       len;
} simplex_t;

/*
 * Simplex structure management interface.
 */
int  simplex_init(simplex_t* simplex, unsigned dimensions);
int  simplex_copy(simplex_t* dest, const simplex_t* src);
void simplex_fini(simplex_t* simplex);

/*
 * Simplex initialization interface.
 */
int simplex_set(simplex_t* simplex, const hspace_t* space,
                const vertex_t* base, double radius);

/*
 * Simplex utility function interface.
 */
int simplex_centroid(const simplex_t* simplex, vertex_t* centroid);
int simplex_collapsed(const simplex_t* simplex, const hspace_t* space);
int simplex_inbounds(const simplex_t* simplex, const hspace_t* space);
int simplex_transform(const simplex_t* base, const vertex_t* origin,
                      double coefficient, simplex_t* result);

#ifdef __cplusplus
}
#endif

#endif
