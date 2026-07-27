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
 * \page nm Nelder-Mead (nm.so)
 *
 * This search strategy uses a simplex-based method to estimate the
 * relative slope of a search space without calculating gradients.  It
 * functions by evaluating the performance for each point of the
 * simplex, and systematically replacing the worst performing point
 * with a reflection, expansion, or contraction in relation to the
 * simplex centroid.  In some cases, the entire simplex may also be
 * shrunken.
 *
 * \note Due to the nature of the underlying algorithm, this strategy
 * is best suited for serial tuning tasks.  It often waits on a single
 * performance report before a new point may be generated.
 *
 * For details of the algorithm, see:
 * > Nelder, John A.; R. Mead (1965). "A simplex method for function
 * >  minimization". Computer Journal 7: 308–313. doi:10.1093/comjnl/7.4.308
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
    { CFGKEY_CONTRACT, "0.5",
      "Multiplicative coefficient for simplex contraction step." },
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
    SIMPLEX_STATE_EXPAND,
    SIMPLEX_STATE_CONTRACT,
    SIMPLEX_STATE_SHRINK,
    SIMPLEX_STATE_CONVERGED,

    SIMPLEX_STATE_MAX
} simplex_state_t;

/*
 * Structure to hold data for an individual Nelder-Mead search instance.
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
    double contract_val;
    double shrink_val;
    double fval_tol;
    double size_tol;

    // Search state.
    vertex_t        centroid;
    vertex_t        reflect;
    vertex_t        expand;
    vertex_t        contract;
    simplex_t       simplex;
    simplex_state_t state;

    vertex_t* next;
    int index_best;
    int index_worst;
    int index_curr; // For INIT or SHRINK.
    int next_id;
};

/*
 * Internal helper function prototypes.
 */
static void check_convergence(hplugin_data_t* data);
static int  config_strategy(hplugin_data_t* data);
static int  nm_algorithm(hplugin_data_t* data);
static int  nm_state_transition(hplugin_data_t* data);
static int  nm_next_vertex(hplugin_data_t* data);
static int  update_centroid(hplugin_data_t* data);

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
    data->space = space;
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

    data->index_curr = 0;
    data->state = SIMPLEX_STATE_INIT;
    if (nm_next_vertex(data) != 0) {
        search_error("Could not initiate test vertex");
        return -1;
    }

    return 0;
}

/*
 * Generate a new candidate configuration point.
 */
int strategy_generate(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    if (data->next->id == data->next_id) {
        flow->status = HFLOW_WAIT;
        return 0;
    }

    data->next->id = data->next_id;
    if (vertex_point(data->next, data->space, point) != 0) {
        search_error("Could not make point from vertex during generate");
        return -1;
    }

    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Regenerate a point deemed invalid by a later plug-in.
 */
int strategy_rejected(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    hpoint_t* hint = &flow->point;

    if (hint->id) {
        // Update our state to include the hint point.
        hint->id = point->id;
        if (vertex_set(data->next, data->space, hint) != 0) {
            search_error("Could not copy hint into simplex during reject");
            return -1;
        }

        if (hpoint_copy(point, hint) != 0) {
            search_error("Could not return hint during reject");
            return -1;
        }
    }
    else if (data->reject_type == REJECT_METHOD_PENALTY) {
        // Apply an infinite penalty to the invalid point and
        // allow the algorithm to determine the next point to try.
        //
        hperf_reset(&data->next->perf);
        if (nm_algorithm(data) != 0) {
            search_error("Nelder-Mead algorithm failure");
            return -1;
        }

        data->next->id = data->next_id;
        if (vertex_point(data->next, data->space, point) != 0) {
            search_error("Could not copy next point during reject");
            return -1;
        }
    }
    else if (data->reject_type == REJECT_METHOD_RANDOM) {
        // Replace the rejected point with a random point.
        if (vertex_random(data->next, data->space, 1.0) != 0) {
            search_error("Could not randomize point during reject");
            return -1;
        }

        data->next->id = data->next_id;
        if (vertex_point(data->next, data->space, point) != 0) {
            search_error("Could not copy random point during reject");
            return -1;
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
    if (trial->point.id != data->next->id)
        return 0;

    if (hperf_copy(&data->next->perf, &trial->perf) != 0) {
        search_error("Could not copy performance to vertex");
        return -1;
    }

    if (nm_algorithm(data) != 0) {
        search_error("Nelder-Mead algorithm failure");
        return -1;
    }

    // Update the best performing point, if necessary.
    if (hperf_cmp(&data->best_perf, &trial->perf) > 0) {
        if (hperf_copy(&data->best_perf, &trial->perf) != 0) {
            search_error("Could not store best performance");
            return -1;
        }

        if (hpoint_copy(&data->best, &trial->point) != 0) {
            search_error("Could not copy best point during analyze");
            return -1;
        }
    }

    if (data->state != SIMPLEX_STATE_CONVERGED)
        ++data->next_id;

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
    simplex_fini(&data->simplex);
    vertex_fini(&data->contract);
    vertex_fini(&data->expand);
    vertex_fini(&data->reflect);
    vertex_fini(&data->centroid);
    vertex_fini(&data->init_point);
    hperf_fini(&data->best_perf);
    hpoint_fini(&data->best);

    free(data);
    return 0;
}

/*
 * Internal helper function implementations.
 */

void check_convergence(hplugin_data_t* data)
{
    double fval_err, size_max;
    double avg_perf = hperf_unify(&data->centroid.perf);

    if (simplex_collapsed(&data->simplex, data->space))
        goto converged;

    fval_err = 0.0;
    for (int i = 0; i < data->simplex.len; ++i) {
        double point_perf = hperf_unify(&data->simplex.vertex[i].perf);
        fval_err += ((point_perf - avg_perf) * (point_perf - avg_perf));
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

    return;

  converged:
    data->state = SIMPLEX_STATE_CONVERGED;
    search_setcfg(CFGKEY_CONVERGED, "1");
}

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
    data->expand_val = hcfg_real(search_cfg, CFGKEY_EXPAND);

    cfgval = hcfg_real(search_cfg, CFGKEY_CONTRACT);
    if (isnan(cfgval) || cfgval <= 0.0 || cfgval >= 1.0) {
        search_error("Configuration key " CFGKEY_CONTRACT
                     " must be between 0.0 and 1.0 (exclusive)");
        return -1;
    }
    data->contract_val = cfgval;

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
    data->fval_tol = cfgval;

    cfgval = hcfg_real(search_cfg, CFGKEY_SIZE_TOL);
    if (isnan(cfgval) || cfgval <= 0.0 || cfgval >= 1.0) {
        search_error("Configuration key " CFGKEY_SIZE_TOL
                     " must be between 0.0 and 1.0 (exclusive)");
        return -1;
    }
    data->size_tol = cfgval;

    // Use the expand and reflect vertex variables as temporaries to
    // calculate the size tolerance.
    if (vertex_minimum(&data->expand, data->space) != 0 ||
        vertex_maximum(&data->reflect, data->space) != 0)
        return -1;

    data->size_tol *= vertex_norm(&data->expand, &data->reflect,
                                  VERTEX_NORM_L2);
    return 0;
}

int nm_algorithm(hplugin_data_t* data)
{
    do {
        if (data->state == SIMPLEX_STATE_CONVERGED)
            break;

        if (nm_state_transition(data) != 0)
            return -1;

        if (data->state == SIMPLEX_STATE_REFLECT) {
            if (update_centroid(data) != 0)
                return -1;

            check_convergence(data);
        }

        if (nm_next_vertex(data) != 0)
            return -1;

    } while (!vertex_inbounds(data->next, data->space));

    return 0;
}

int nm_state_transition(hplugin_data_t* data)
{
    switch (data->state) {
    case SIMPLEX_STATE_INIT:
    case SIMPLEX_STATE_SHRINK:
        // Simplex vertex performance value.
        if (++data->index_curr == data->space->len + 1)
            data->state = SIMPLEX_STATE_REFLECT;

        break;

    case SIMPLEX_STATE_REFLECT:
        if (hperf_cmp(&data->reflect.perf,
                      &data->simplex.vertex[data->index_best].perf) < 0)
        {
            // Reflected point performs better than all simplex points.
            // Attempt expansion.
            //
            data->state = SIMPLEX_STATE_EXPAND;
        }
        else if (hperf_cmp(&data->reflect.perf,
                           &data->simplex.vertex[data->index_worst].perf) < 0)
        {
            // Reflected point performs better than worst simplex point.
            // Replace the worst simplex point with reflected point
            // and attempt reflection again.
            //
            if (vertex_copy(&data->simplex.vertex[data->index_worst],
                            &data->reflect) != 0)
                return -1;
        }
        else {
            // Reflected point does not improve the current simplex.
            // Attempt contraction.
            //
            data->state = SIMPLEX_STATE_CONTRACT;
        }
        break;

    case SIMPLEX_STATE_EXPAND:
        if (hperf_cmp(&data->expand.perf, &data->reflect.perf) < 0) {
            // Expanded point performs even better than reflected point.
            // Replace the worst simplex point with the expanded point
            // and attempt reflection again.
            //
            if (vertex_copy(&data->simplex.vertex[data->index_worst],
                            &data->expand) != 0)
                return -1;
        }
        else {
            // Expanded point did not result in improved performance.
            // Replace the worst simplex point with the original
            // reflected point and attempt reflection again.
            //
            if (vertex_copy(&data->simplex.vertex[data->index_worst],
                            &data->reflect) != 0)
                return -1;
        }
        data->state = SIMPLEX_STATE_REFLECT;
        break;

    case SIMPLEX_STATE_CONTRACT:
        if (hperf_cmp(&data->contract.perf,
                      &data->simplex.vertex[data->index_worst].perf) < 0)
        {
            // Contracted point performs better than the worst simplex point.
            //
            // Replace the worst simplex point with contracted point
            // and attempt reflection.
            //
            if (vertex_copy(&data->simplex.vertex[data->index_worst],
                            &data->contract) != 0)
                return -1;

            data->state = SIMPLEX_STATE_REFLECT;
        }
        else {
            // Contracted test vertex has worst known performance.
            // Shrink the entire simplex towards the best point.
            //
            data->index_curr = -1; // Indicates the beginning of SHRINK.
            data->state = SIMPLEX_STATE_SHRINK;
        }
        break;

    default:
        return -1;
    }
    return 0;
}

int nm_next_vertex(hplugin_data_t* data)
{
    switch (data->state) {
    case SIMPLEX_STATE_INIT:
        // Test individual vertices of the initial simplex.
        data->next = &data->simplex.vertex[data->index_curr];
        break;

    case SIMPLEX_STATE_REFLECT:
        // Test a vertex reflected from the worst performing vertex
        // through the centroid point.
        //
        if (vertex_transform(&data->centroid,
                             &data->simplex.vertex[data->index_worst],
                             data->reflect_val, &data->reflect) != 0)
            return -1;

        data->next = &data->reflect;
        break;

    case SIMPLEX_STATE_EXPAND:
        // Test a vertex that expands the reflected vertex even
        // further from the the centroid point.
        //
        if (vertex_transform(&data->centroid,
                             &data->simplex.vertex[data->index_worst],
                             data->expand_val, &data->expand) != 0)
            return -1;

        data->next = &data->expand;
        break;

    case SIMPLEX_STATE_CONTRACT:
        // Test a vertex contracted from the worst performing vertex
        // towards the centroid point.
        //
        if (vertex_transform(&data->simplex.vertex[data->index_worst],
                             &data->centroid,
                             -data->contract_val, &data->contract) != 0)
            return -1;

        data->next = &data->contract;
        break;

    case SIMPLEX_STATE_SHRINK:
        if (data->index_curr == -1) {
            // Shrink the entire simplex towards the best known vertex
            // thus far.
            //
            if (simplex_transform(&data->simplex,
                                  &data->simplex.vertex[data->index_best],
                                  -data->shrink_val, &data->simplex) != 0)
                return -1;

            data->index_curr = 0;
        }

        // Test individual vertices of the initial simplex.
        data->next = &data->simplex.vertex[data->index_curr];
        break;

    case SIMPLEX_STATE_CONVERGED:
        // Simplex has converged.  Nothing to do.
        // In the future, we may consider new search at this point.
        //
        data->next = &data->simplex.vertex[data->index_best];
        break;

    default:
        return -1;
    }

    hperf_reset(&data->next->perf);
    return 0;
}

int update_centroid(hplugin_data_t* data)
{
    data->index_best = 0;
    data->index_worst = 0;

    for (int i = 1; i < data->simplex.len; ++i) {
        if (hperf_cmp(&data->simplex.vertex[i].perf,
                      &data->simplex.vertex[data->index_best].perf) < 0)
            data->index_best = i;

        if (hperf_cmp(&data->simplex.vertex[i].perf,
                      &data->simplex.vertex[data->index_worst].perf) > 0)
            data->index_worst = i;
    }

    unsigned stashed_id = data->simplex.vertex[data->index_worst].id;
    data->simplex.vertex[data->index_worst].id = 0;
    if (simplex_centroid(&data->simplex, &data->centroid) != 0)
        return -1;

    data->simplex.vertex[data->index_worst].id = stashed_id;
    return 0;
}
