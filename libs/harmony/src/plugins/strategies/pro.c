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

/**
 * \page pro Parallel Rank Order (pro.so)
 *
 * This search strategy uses a simplex-based method similar to the
 * Nelder-Mead algorithm.  It improves upon the Nelder-Mead algorithm
 * by allowing the simultaneous search of all simplex points at each
 * step of the algorithm.  As such, it is ideal for a parallel search
 * utilizing multiple nodes, for instance when integrated in OpenMP or
 * MPI programs.
 */

#include "hstrategy.h"
#include "session-core.h"
#include "hcfg.h"
#include "hspace.h"
#include "hpoint.h"
#include "hperf.h"
#include "hutil.h"
#include "libvertex.h"

#include <string.h> // For strcmp().
#include <math.h>   // For isnan().

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_INIT_POINT, NULL,
      "Centroid point used to initialize the search simplex.  If this key "
      "is left undefined, the simplex will be initialized in the center of "
      "the search space." },
    { CFGKEY_INIT_RADIUS, "0.50",
      "Size of the initial simplex, specified as a fraction of the total "
      "search space radius." },
    { CFGKEY_REJECT_METHOD, "penalty",
      "How to choose a replacement when dealing with rejected points. "
      "    penalty: Use this method if the chance of point rejection is "
      "relatively low. It applies an infinite penalty factor for invalid "
      "points, allowing the strategy to select a sensible next point.  "
      "However, if the entire simplex is comprised of invalid points, an "
      "infinite loop of rejected points may occur.\n"
      "    random: Use this method if the chance of point rejection is "
      "high.  It reduces the risk of infinitely selecting invalid points "
      "at the cost of increasing the risk of deforming the simplex." },
    { CFGKEY_REFLECT, "1.0",
      "Multiplicative coefficient for simplex reflection step." },
    { CFGKEY_EXPAND, "2.0",
      "Multiplicative coefficient for simplex expansion step." },
    { CFGKEY_SHRINK, "0.5",
      "Multiplicative coefficient for simplex shrink step." },
    { CFGKEY_FVAL_TOL, "0.0001",
      "Convergence test succeeds if difference between all vertex "
      "performance values fall below this value." },
    { CFGKEY_SIZE_TOL, "0.005",
      "Convergence test succeeds if the simplex radius becomes smaller "
      "than this percentage of the total search space.  Simplex radius "
      "is measured from centroid to furthest vertex.  Total search space "
      "is measured from minimum to maximum point." },
    { NULL }
};

typedef enum reject_method {
    REJECT_METHOD_UNKNOWN = 0,
    REJECT_METHOD_PENALTY,
    REJECT_METHOD_RANDOM,

    REJECT_METHOD_MAX
} reject_method_t;

typedef enum simplex_state {
    SIMPLEX_STATE_UNKNONW = 0,
    SIMPLEX_STATE_INIT,
    SIMPLEX_STATE_REFLECT,
    SIMPLEX_STATE_EXPAND_ONE,
    SIMPLEX_STATE_EXPAND_ALL,
    SIMPLEX_STATE_SHRINK,
    SIMPLEX_STATE_CONVERGED,

    SIMPLEX_STATE_MAX
} simplex_state_t;

/*
 * Structure to hold data for an individual PRO search instance.
 */
struct hplugin_data {
    hspace_t* space;
    hpoint_t  best;
    hperf_t   best_perf;

    // Search options.
    vertex_t        init_point;
    double          init_radius;
    reject_method_t reject_type;

    double reflect_val;
    double expand_val;
    double shrink_val;
    double fval_tol;
    double size_tol;

    // Search state.
    simplex_t       simplex;
    simplex_t       reflect;
    simplex_t       expand;
    vertex_t        centroid;
    simplex_state_t state;

    simplex_t* next;
    int        best_base;
    int        best_test;
    int        best_reflect;
    int        next_id;
    int        send_idx;
    int        reported;
};

/*
 * Internal helper function prototypes.
 */
static int config_strategy(hplugin_data_t* data);
static int pro_algorithm(hplugin_data_t* data);
static int pro_next_state(hplugin_data_t* data);
static int pro_next_simplex(hplugin_data_t* data);
static int check_convergence(hplugin_data_t* data);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* strategy_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    retval->next_id = 1;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search task.
 */
int strategy_init(hplugin_data_t* data, hspace_t* space)
{
    if (data->space != space) {
        if (simplex_init(&data->simplex, space->len) != 0) {
            search_error("Could not initialize base simplex");
            return -1;
        }

        if (simplex_init(&data->reflect, space->len) != 0) {
            search_error("Could not initialize reflection simplex");
            return -1;
        }

        if (simplex_init(&data->expand, space->len) != 0) {
            search_error("Could not initialize expansion simplex");
            return -1;
        }
        data->space = space;
    }

    if (config_strategy(data) != 0)
        return -1;

    if (simplex_set(&data->simplex, space,
                    &data->init_point, data->init_radius) != 0)
    {
        search_error("Could not generate initial simplex");
        return -1;
    }

    if (search_setcfg(CFGKEY_CONVERGED, "0") != 0) {
        search_error("Could not set " CFGKEY_CONVERGED " config variable");
        return -1;
    }

    data->send_idx = 0;
    data->state = SIMPLEX_STATE_INIT;
    if (pro_next_simplex(data) != 0) {
        search_error("Could not initiate the simplex");
        return -1;
    }

    data->reported = 0;
    return 0;
}

/*
 * Generate a new candidate configuration point.
 */
int strategy_generate(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    if (data->send_idx > data->space->len ||
        data->state == SIMPLEX_STATE_CONVERGED)
    {
        flow->status = HFLOW_WAIT;
        return 0;
    }

    data->next->vertex[data->send_idx].id = data->next_id;
    if (vertex_point(&data->next->vertex[data->send_idx],
                     data->space, point) != 0)
    {
        search_error("Could not copy point during generate");
        return -1;
    }
    ++data->next_id;
    ++data->send_idx;

    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Regenerate a point deemed invalid by a later plug-in.
 */
int strategy_rejected(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    int reject_idx;
    hpoint_t* hint = &flow->point;

    // Find the rejected vertex.
    for (reject_idx = 0; reject_idx <= data->space->len; ++reject_idx) {
        if (data->next->vertex[reject_idx].id == point->id)
            break;
    }
    if (reject_idx > data->space->len) {
        search_error("Could not find rejected point");
        return -1;
    }

    if (hint && hint->id) {
        // Update our state to include the hint point.
        if (vertex_set(&data->next->vertex[reject_idx],
                       data->space, hint) != 0)
        {
            search_error("Could not copy hint into simplex during reject");
            return -1;
        }

        // Return the hint point.
        if (hpoint_copy(point, hint) != 0) {
            search_error("Could not return hint during reject");
            return -1;
        }
    }
    else {
        if (data->reject_type == REJECT_METHOD_PENALTY) {
            // Apply an infinite penalty to the invalid point and
            // allow the algorithm to determine the next point to try.
            //
            hperf_reset(&data->next->vertex[reject_idx].perf);
            ++data->reported;

            if (data->reported > data->space->len) {
                if (pro_algorithm(data) != 0) {
                    search_error("PRO algorithm failure");
                    return -1;
                }
                data->reported = 0;
                data->send_idx = 0;
            }

            if (data->send_idx > data->space->len) {
                flow->status = HFLOW_WAIT;
                return 0;
            }

            data->next->vertex[data->send_idx].id = data->next_id;
            if (vertex_point(&data->next->vertex[data->send_idx],
                             data->space, point) != 0)
            {
                search_error("Could not convert vertex during reject");
                return -1;
            }

            ++data->next_id;
            ++data->send_idx;
        }
        else if (data->reject_type == REJECT_METHOD_RANDOM) {
            // Replace the rejected point with a random point.
            if (vertex_random(&data->next->vertex[reject_idx],
                              data->space, 1.0) != 0)
            {
                search_error("Could not make random point during reject");
                return -1;
            }

            if (vertex_point(&data->next->vertex[reject_idx],
                             data->space, point) != 0)
            {
                search_error("Could not convert vertex during reject");
                return -1;
            }
        }
    }
    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Analyze the observed performance for this configuration point.
 */
int strategy_analyze(hplugin_data_t* data, htrial_t* trial)
{
    int report_idx;
    for (report_idx = 0; report_idx <= data->space->len; ++report_idx) {
        if (data->next->vertex[report_idx].id == trial->point.id)
            break;
    }
    if (report_idx > data->space->len) {
        // Ignore rouge vertex reports.
        return 0;
    }

    ++data->reported;
    hperf_copy(&data->next->vertex[report_idx].perf, &trial->perf);
    if (hperf_cmp(&data->next->vertex[report_idx].perf,
                  &data->next->vertex[data->best_test].perf) < 0)
        data->best_test = report_idx;

    if (data->reported > data->space->len) {
        if (pro_algorithm(data) != 0) {
            search_error("PRO algorithm failure");
            return -1;
        }
        data->reported = 0;
        data->send_idx = 0;
    }

    // Update the best performing point, if necessary.
    if (!data->best.id || hperf_cmp(&data->best_perf, &trial->perf) > 0) {
        if (hperf_copy(&data->best_perf, &trial->perf) != 0) {
            search_error("Could not store best performance during analyze");
            return -1;
        }

        if (hpoint_copy(&data->best, &trial->point) != 0) {
            search_error("Could not copy best point during analyze");
            return -1;
        }
    }
    return 0;
}

/*
 * Return the best performing point thus far in the search.
 */
int strategy_best(hplugin_data_t* data, hpoint_t* point)
{
    if (hpoint_copy(point, &data->best) != 0) {
        search_error("Could not copy best point during strategy_best()");
        return -1;
    }
    return 0;
}

/*
 * Free memory associated with this search task.
 */
int strategy_fini(hplugin_data_t* data)
{
    vertex_fini(&data->centroid);
    simplex_fini(&data->expand);
    simplex_fini(&data->reflect);
    simplex_fini(&data->simplex);
    vertex_fini(&data->init_point);
    hperf_fini(&data->best_perf);
    hpoint_fini(&data->best);

    free(data);
    return 0;
}

/*
 * Internal helper function implementation.
 */

int config_strategy(hplugin_data_t* data)
{
    const char* cfgstr;
    double cfgval;

    cfgstr = hcfg_get(search_cfg, CFGKEY_INIT_POINT);
    if (cfgstr) {
        if (vertex_parse(&data->init_point, data->space, cfgstr) != 0) {
            search_error("Could not convert initial point to vertex");
            return -1;
        }
    }
    else {
        if (vertex_center(&data->init_point, data->space) != 0) {
            search_error("Could not create central vertex");
            return -1;
        }
    }

    cfgval = hcfg_real(search_cfg, CFGKEY_INIT_RADIUS);
    if (isnan(cfgval) || cfgval <= 0 || cfgval > 1) {
        search_error("Configuration key " CFGKEY_INIT_RADIUS
                     " must be between 0.0 and 1.0 (exclusive)");
        return -1;
    }
    data->init_radius = cfgval;

    cfgstr = hcfg_get(search_cfg, CFGKEY_REJECT_METHOD);
    if (cfgstr) {
        if (strcmp(cfgstr, "penalty") == 0) {
            data->reject_type = REJECT_METHOD_PENALTY;
        }
        else if (strcmp(cfgstr, "random") == 0) {
            data->reject_type = REJECT_METHOD_RANDOM;
        }
        else {
            search_error("Invalid value for "
                         CFGKEY_REJECT_METHOD " configuration key");
            return -1;
        }
    }

    cfgval = hcfg_real(search_cfg, CFGKEY_REFLECT);
    if (isnan(cfgval) || cfgval <= 0.0) {
        search_error("Configuration key " CFGKEY_REFLECT
                     " must be positive");
        return -1;
    }
    data->reflect_val = cfgval;

    cfgval = hcfg_real(search_cfg, CFGKEY_EXPAND);
    if (isnan(cfgval) || cfgval <= data->reflect_val) {
        search_error("Configuration key " CFGKEY_EXPAND
                     " must be greater than the reflect coefficient");
        return -1;
    }
    data->expand_val = cfgval;

    cfgval = hcfg_real(search_cfg, CFGKEY_SHRINK);
    if (isnan(cfgval) || cfgval <= 0.0 || cfgval >= 1.0) {
        search_error("Configuration key " CFGKEY_SHRINK
                     " must be between 0.0 and 1.0 (exclusive)");
        return -1;
    }
    data->shrink_val = cfgval;

    cfgval = hcfg_real(search_cfg, CFGKEY_FVAL_TOL);
    if (isnan(cfgval)) {
        search_error("Configuration key " CFGKEY_FVAL_TOL " is invalid");
        return -1;
    }
    data->fval_tol = hcfg_real(search_cfg, CFGKEY_FVAL_TOL);

    cfgval = hcfg_real(search_cfg, CFGKEY_SIZE_TOL);
    if (isnan(cfgval) || cfgval <= 0.0 || cfgval >= 1.0) {
        search_error("Configuration key " CFGKEY_SIZE_TOL
                     " must be between 0.0 and 1.0 (exclusive)");
        return -1;
    }
    data->size_tol = cfgval;

    // Use the first two verticies of the base simplex as temporaries
    // to calculate the size tolerance.
    vertex_t* vertex = data->simplex.vertex;
    if (vertex_minimum(&vertex[0], data->space) != 0 ||
        vertex_maximum(&vertex[1], data->space) != 0)
        return -1;

    data->size_tol *= vertex_norm(&vertex[0], &vertex[1], VERTEX_NORM_L2);
    return 0;
}

int pro_algorithm(hplugin_data_t* data)
{
    do {
        if (data->state == SIMPLEX_STATE_CONVERGED)
            break;

        if (pro_next_state(data) != 0)
            return -1;

        if (data->state == SIMPLEX_STATE_REFLECT) {
            if (check_convergence(data) != 0)
                return -1;
        }

        if (pro_next_simplex(data) != 0)
            return -1;

    } while (!simplex_inbounds(data->next, data->space));

    return 0;
}

int pro_next_state(hplugin_data_t* data)
{
    switch (data->state) {
    case SIMPLEX_STATE_INIT:
    case SIMPLEX_STATE_SHRINK:
        // Simply accept the candidate simplex and prepare to reflect.
        data->best_base = data->best_test;
        data->state = SIMPLEX_STATE_REFLECT;
        break;

    case SIMPLEX_STATE_REFLECT:
        if (hperf_cmp(&data->reflect.vertex[data->best_test].perf,
                      &data->simplex.vertex[data->best_base].perf) < 0)
        {
            // Reflected simplex has best known performance.
            // Attempt a trial expansion.
            //
            data->best_reflect = data->best_test;
            data->state = SIMPLEX_STATE_EXPAND_ONE;
        }
        else {
            // Reflected simplex does not improve performance.
            // Shrink the simplex instead.
            //
            data->state = SIMPLEX_STATE_SHRINK;
        }
        break;

    case SIMPLEX_STATE_EXPAND_ONE:
        if (hperf_cmp(&data->expand.vertex[0].perf,
                      &data->reflect.vertex[data->best_reflect].perf) < 0)
        {
            // Trial expansion has found the best known vertex thus far.
            // We are now free to expand the entire reflected simplex.
            //
            data->state = SIMPLEX_STATE_EXPAND_ALL;
        }
        else {
            // Expanded test vertex does not improve performance.
            // Accept the original (unexpanded) reflected simplex and
            // attempt another reflection.
            //
            if (simplex_copy(&data->simplex, &data->reflect) != 0) {
                search_error("Could not copy reflection simplex"
                             " for next state");
                return -1;
            }
            data->best_base = data->best_reflect;
            data->state = SIMPLEX_STATE_REFLECT;
        }
        break;

    case SIMPLEX_STATE_EXPAND_ALL:
        if (simplex_inbounds(&data->expand, data->space)) {
            // If the entire expanded simplex is valid (in bounds),
            // accept it as the reference simplex.
            //
            if (simplex_copy(&data->simplex, &data->expand) != 0) {
                search_error("Could not copy expanded simplex"
                             " for next state");
                return -1;
            }
            data->best_base = data->best_test;
        }
        else {
            // Otherwise, accept the original (unexpanded) reflected
            // simplex as the reference simplex.
            //
            if (simplex_copy(&data->simplex, &data->reflect) != 0) {
                search_error("Could not copy reflected simplex"
                             " for next state");
                return -1;
            }
            data->best_base = data->best_reflect;
        }

        // Either way, test reflection next.
        data->state = SIMPLEX_STATE_REFLECT;
        break;

    default:
        return -1;
    }
    return 0;
}

int pro_next_simplex(hplugin_data_t* data)
{
    switch (data->state) {
    case SIMPLEX_STATE_INIT:
        // Bootstrap the process by testing the reference simplex.
        data->next = &data->simplex;
        break;

    case SIMPLEX_STATE_REFLECT:
        // Each simplex vertex is translated to a position opposite
        // its origin (best performing) vertex.  This corresponds to a
        // coefficient that is -1 * (1.0 + <reflect coefficient>).
        //
        simplex_transform(&data->simplex,
                          &data->simplex.vertex[data->best_base],
                          -(1.0 + data->reflect_val), &data->reflect);
        data->next = &data->reflect;
        break;

    case SIMPLEX_STATE_EXPAND_ONE:
        // Next simplex should have one vertex extending the best.
        // And the rest should be copies of the best known vertex.
        //
        vertex_transform(&data->simplex.vertex[data->best_reflect],
                         &data->simplex.vertex[data->best_base],
                         -(1.0 + data->expand_val), &data->expand.vertex[0]);

        for (int i = 1; i < data->expand.len; ++i)
            vertex_copy(&data->expand.vertex[i],
                        &data->simplex.vertex[data->best_base]);

        data->next = &data->expand;
        break;

    case SIMPLEX_STATE_EXPAND_ALL:
        // Expand all original simplex vertices away from the best
        // known vertex thus far.
        //
        simplex_transform(&data->simplex,
                          &data->simplex.vertex[data->best_base],
                          -(1.0 + data->expand_val), &data->expand);
        data->next = &data->expand;
        break;

    case SIMPLEX_STATE_SHRINK:
        // Shrink all original simplex vertices towards the best
        // known vertex thus far.
        //
        simplex_transform(&data->simplex,
                          &data->simplex.vertex[data->best_base],
                          -data->shrink_val, &data->simplex);
        data->next = &data->simplex;
        break;

    case SIMPLEX_STATE_CONVERGED:
        // Simplex has converged.  Nothing to do.
        // In the future, we may consider new search at this point.
        //
        break;

    default:
        return -1;
    }
    return 0;
}

int check_convergence(hplugin_data_t* data)
{
    double avg_perf;
    double fval_err;
    double size_max;

    if (simplex_centroid(&data->simplex, &data->centroid) != 0)
        return -1;

    if (simplex_collapsed(&data->simplex, data->space))
        goto converged;

    avg_perf = hperf_unify(&data->centroid.perf);
    fval_err = 0.0;
    for (int i = 0; i < data->simplex.len; ++i) {
        double vert_perf = hperf_unify(&data->simplex.vertex[i].perf);
        fval_err += ((vert_perf - avg_perf) * (vert_perf - avg_perf));
    }
    fval_err /= data->simplex.len;

    size_max = 0.0;
    for (int i = 0; i < data->simplex.len; ++i) {
        double dist = vertex_norm(&data->simplex.vertex[i], &data->centroid,
                                  VERTEX_NORM_L2);
        if (size_max < dist)
            size_max = dist;
    }

    if (fval_err < data->fval_tol && size_max < data->size_tol)
        goto converged;

    return 0;

  converged:
    data->state = SIMPLEX_STATE_CONVERGED;
    search_setcfg(CFGKEY_CONVERGED, "1");
    return 0;
}
