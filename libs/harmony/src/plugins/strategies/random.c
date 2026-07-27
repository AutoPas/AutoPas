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
 * \page random Random (random.so)
 *
 * This search strategy generates random points within the search
 * space.  Using* a pseudo-random method, a value is selected for each
 * tuning variable according to its defined bounds.  This search will
 * never reach a converged state.
 *
 * It is mainly used as a basis of comparison for more intelligent
 * search strategies.
 */

#include "hstrategy.h"
#include "session-core.h"
#include "hcfg.h"
#include "hspace.h"
#include "hpoint.h"
#include "hperf.h"
#include "hutil.h"

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_INIT_POINT, NULL, "Initial point begin testing from." },
    { NULL }
};

/*
 * Structure to hold data for an individual exhaustive search instance.
 */
struct hplugin_data {
    int       space_id;
    hspace_t* space;
    hpoint_t  best;
    double    best_perf;

    hpoint_t  next;
};

/*
 * Internal helper function prototypes.
 */
static int  config_strategy(hplugin_data_t* data);
static void randomize(hplugin_data_t* data, hpoint_t* point);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* strategy_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    retval->best_perf = HUGE_VAL;
    retval->next.id = 1;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search task.
 */
int strategy_init(hplugin_data_t* data, hspace_t* space)
{
    if (data->space_id != space->id) {
        if (hpoint_init(&data->next, space->len) != 0) {
            search_error("Could not initialize point structure");
            return -1;
        }
        data->next.len = space->len;
        data->space = space;
    }

    if (config_strategy(data) != 0)
        return -1;

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
    if (hpoint_copy(point, &data->next) != 0) {
        search_error("Could not copy point during generation");
        return -1;
    }

    // Prepare a new random vertex for the next call to strategy_generate.
    randomize(data, &data->next);
    ++data->next.id;

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
            search_error("Could not copy hint point during reject");
            return -1;
        }
    }
    else {
        randomize(data, point);
    }

    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Analyze the observed performance for this configuration point.
 */
int strategy_analyze(hplugin_data_t* data, htrial_t* trial)
{
    double perf = hperf_unify(&trial->perf);

    if (data->best_perf > perf) {
        data->best_perf = perf;
        if (hpoint_copy(&data->best, &trial->point) != 0) {
            search_error("Could not copy best point");
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
        search_error("Could not copy best point out of strategy");
        return -1;
    }
    return 0;
}

/*
 * Free memory associated with this search task.
 */
int strategy_fini(hplugin_data_t* data)
{
    hpoint_fini(&data->next);
    hpoint_fini(&data->best);

    free(data);
    return 0;
}

/*
 * Internal helper function implementation.
 */

int config_strategy(hplugin_data_t* data)
{
    const char* cfgval = hcfg_get(search_cfg, CFGKEY_INIT_POINT);
    if (cfgval) {
        if (hpoint_parse(&data->next, cfgval, data->space) != 0) {
            search_error("Error parsing point from " CFGKEY_INIT_POINT);
            return -1;
        }

        if (hpoint_align(&data->next, data->space) != 0) {
            search_error("Could not align initial point to search space");
            return -1;
        }
    }
    else {
        randomize(data, &data->next);
    }
    return 0;
}

void randomize(hplugin_data_t* data, hpoint_t* point)
{
    for (int i = 0; i < data->space->len; ++i)
        point->term[i] = hrange_random(&data->space->dim[i], search_drand48());
}
