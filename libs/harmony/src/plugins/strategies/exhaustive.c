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
 * \page exhaustive Exhaustive (exhaustive.so)
 *
 * This search strategy starts with the minimum-value point (i.e.,
 * using the minimum value for each tuning variable), and incrementing
 * the tuning variables like an odometer until the maximum-value point
 * is reached.  This strategy is guaranteed to visit all points within
 * a search space.
 *
 * It is mainly used as a basis of comparison for more intelligent
 * search strategies.
 */

#include "hstrategy.h"
#include "session-core.h"
#include "hutil.h"

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_PASSES, "1",
      "Number of passes through the search space before the search "
      "is considered converged." },
    { CFGKEY_INIT_POINT, NULL,
      "Initial point begin testing from." },
    { NULL }
};

/*
 * Index of an individual hpoint_t term.  The index may be represented
 * as a double, in the case of non-finite value ranges.
 */
typedef union unit {
    unsigned long index;
    double        value;
} unit_u;

/*
 * Structure to hold data for an individual exhaustive search instance.
 */
struct hplugin_data {
    int       space_id;
    hspace_t* space;
    hpoint_t  best;
    double    best_perf;

    unit_u*  head;
    unit_u*  next;
    unsigned next_id;
    unit_u*  wrap;

    int remaining_passes;
    int final_id;
    int outstanding_points;
    int final_point_received;
};

/*
 * Internal helper function prototypes.
 */
static int  config_strategy(hplugin_data_t* data);
static void increment(hplugin_data_t* data);
static int  make_next_point(hplugin_data_t* data, hpoint_t* point);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* strategy_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    retval->best_perf = HUGE_VAL;
    retval->next_id = 1;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search task.
 */
int strategy_init(hplugin_data_t* data, hspace_t* space)
{
    if (data->space_id != space->id) {
        free(data->head);
        data->head = malloc(space->len * sizeof(*data->head));
        if (!data->head) {
            search_error("Could not allocate head index array");
            return -1;
        }

        free(data->next);
        data->next = malloc(space->len * sizeof(*data->next));
        if (!data->next) {
            search_error("Could not allocate next index array");
            return -1;
        }

        free(data->wrap);
        data->wrap = malloc(space->len * sizeof(*data->wrap));
        if (!data->wrap) {
            search_error("Could not allocate wrap index array");
            return -1;
        }
        data->space = space;
    }

    if (config_strategy(data) != 0)
        return -1;

    // Determine each search dimension's upper limit.
    for (int i = 0; i < space->len; ++i) {
        if (hrange_finite(&space->dim[i]))
            data->wrap[i].index = hrange_limit(&space->dim[i]);
        else
            data->wrap[i].value = space->dim[i].bounds.r.max;
    }

    // Initialize the next point.
    memcpy(data->next, data->head, space->len * sizeof(*data->next));

    if (search_setcfg(CFGKEY_CONVERGED, "0") != 0) {
        search_error("Could not set " CFGKEY_CONVERGED " config variable");
        return -1;
    }
    return 0;
}

/*
 * Generate a new candidate configuration.
 */
int strategy_generate(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    if (data->remaining_passes > 0) {
        if (make_next_point(data, point) != 0) {
            search_error("Could not make point from index during generate");
            return -1;
        }
        point->id = data->next_id;

        increment(data);
        ++data->next_id;
    }
    else {
        if (hpoint_copy(point, &data->best) != 0) {
            search_error("Could not copy best point during generation");
            return -1;
        }
    }

    // Every time we send out a point that's before the final point,
    // increment the numebr of points we're waiting for results from.
    if (!data->final_id || data->next_id <= data->final_id)
        ++data->outstanding_points;

    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Regenerate a point deemed invalid by a later plug-in.
 */
int strategy_rejected(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    if (flow->point.id) {
        hpoint_t* hint = &flow->point;

        hint->id = point->id;
        if (hpoint_copy(point, hint) != 0) {
            search_error("Could not copy hint during reject");
            return -1;
        }
    }
    else {
        if (make_next_point(data, point) != 0) {
            search_error("Could not make point from index during reject");
            return -1;
        }
        increment(data);
    }

    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Analyze the observed performance for this configuration point.
 */
int strategy_analyze(hplugin_data_t* data, htrial_t* trial)
{
    // Function local variables.
    double perf = hperf_unify(&trial->perf);

    if (data->best_perf > perf) {
        data->best_perf = perf;
        if (hpoint_copy(&data->best, &trial->point) != 0) {
            search_error("Could not copy best point during analyze");
            return -1;
        }
    }

    if (trial->point.id == data->final_id) {
        if (search_setcfg(CFGKEY_CONVERGED, "1") != 0) {
            search_error("Could not set convergence status");
            return -1;
        }
    }

    // Decrement the number of points we're waiting for when we get a
    // point back that was generated before the final point.
    if (!data->final_id || trial->point.id <= data->final_id)
        --data->outstanding_points;

    if (trial->point.id == data->final_id)
        data->final_point_received = 1;

    // Converged when the final point has been received, and there are
    // no outstanding points.
    if (data->outstanding_points <= 0 && data->final_point_received) {
        if (search_setcfg(CFGKEY_CONVERGED, "1") != 0) {
            search_error("Could not set convergence status");
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
        search_error("Could not copy best point during request for best");
        return -1;
    }
    return 0;
}

/*
 * Free memory associated with this search task.
 */
int strategy_fini(hplugin_data_t* data)
{
    free(data->wrap);
    free(data->next);
    free(data->head);
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

    data->remaining_passes = hcfg_int(search_cfg, CFGKEY_PASSES);
    if (data->remaining_passes < 0) {
        search_error("Invalid value for " CFGKEY_PASSES);
        return -1;
    }

    cfgstr = hcfg_get(search_cfg, CFGKEY_INIT_POINT);
    if (cfgstr) {
        hpoint_t init;

        if (hpoint_parse(&init, cfgstr, data->space) != 0) {
            search_error("Error parsing point from " CFGKEY_INIT_POINT);
            return -1;
        }

        if (!hpoint_align(&init, data->space) != 0) {
            search_error("Could not align initial point to search space");
            return -1;
        }

        for (int i = 0; i < data->space->len; ++i) {
            if (hrange_finite(&data->space->dim[i]))
                data->head[i].index = hrange_index(&data->space->dim[i],
                                                   &init.term[i]);
            else
                data->head[i].value = init.term[i].value.r;
        }
        hpoint_fini(&init);
    }
    else {
        memset(data->head, 0, data->space->len * sizeof(*data->head));
    }
    return 0;
}

void increment(hplugin_data_t* data)
{
    if (data->remaining_passes <= 0)
        return;

    for (int i = 0; i < data->space->len; ++i) {
        if (hrange_finite(&data->space->dim[i])) {
            ++data->next[i].index;
            if (data->next[i].index == data->wrap[i].index) {
                data->next[i].index = 0;
                continue; // Overflow detected.
            }
        }
        else {
            double nextval = nextafter(data->next[i].value, HUGE_VAL);
            if (!(data->next[i].value < nextval)) {
                data->next[i].value = data->space->dim[i].bounds.r.min;
                continue; // Overflow detected.
            }
            data->next[i].value = nextval;
        }
        return; // No overflow detected.  Exit function.
    }

    // All values overflowed.
    if (--data->remaining_passes <= 0)
        data->final_id = data->next_id;
}

int make_next_point(hplugin_data_t* data, hpoint_t* point)
{
    if (point->cap < data->space->len) {
        if (hpoint_init(point, data->space->len) != 0)
            return -1;
    }

    for (int i = 0; i < data->space->len; ++i) {
        if (hrange_finite(&data->space->dim[i]))
            point->term[i] = hrange_value(&data->space->dim[i],
                                          data->next[i].index);
        else {
            point->term[i].type = HVAL_REAL;
            point->term[i].value.r = data->next[i].value;
        }
    }

    point->len = data->space->len;
    return 0;
}
