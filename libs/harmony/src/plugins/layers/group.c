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
 * \page group Input Grouping (group.so)
 *
 * This processing layer allows search space input variables
 * (dimensions) to be iteratively search in groups.  The value for all
 * input variables not in the current search group remain constant.
 *
 * Input variables are specified via a zero-based index, determined by
 * the order they are defined in the session search space.  For instance,
 * the first variable defined via harmony_int(), harmony_real(), or
 * harmony_enum() can be referenced via the index 0.
 *
 * ### Example ###
 * Given a tuning session with nine input variables, the following
 * group specification:
 *
 *     (0,1,2),(3,4,5,6),(7,8)
 *
 * would instruct the layer to search the first three variables until
 * the search strategy converges.  Input values for the other six
 * variables remain constant during this search.
 *
 * Upon convergence, the search is restarted and allowed to
 * investigate the next four variables.  The first three variables are
 * forced to use the best discovered values from the prior search.
 * This pattern continues until all groups have been searched.
 *
 * Grouping should be beneficial for search spaces whose input
 * variables are relatively independent with respect to the reported
 * performance.
 */

#include "hlayer.h"
#include "session-core.h"
#include "hspace.h"
#include "hpoint.h"
#include "hutil.h"
#include "hcfg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>

/*
 * Name used to identify this plugin layer.
 * All Harmony plugin layers must define this variable.
 */
const char hplugin_name[] = "group";

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_GROUP_LIST, NULL, "List of input parameter indexes." },
    { NULL }
};

typedef struct group_def {
    int* idx;
    int  idx_len;
} group_def_t;

/*
 * Structure to hold all data needed by an individual search instance.
 *
 * To support multiple parallel search instances, no global variables
 * should be defined or used in this plug-in layer.  They should
 * instead be defined as a part of this structure.
 */
struct hplugin_data {
    char internal_restart_req;
    int cap_max;
    hval_t* locked_val;
    hval_t* hint_val;

    group_def_t* glist;
    int glist_len, glist_cap;
    int glist_curr;

    hpoint_t best;
};

/*
 * Internal helper function prototypes.
 */
static int parse_group(hplugin_data_t* data, const char* buf);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* group_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search task.
 */
int group_init(hplugin_data_t* data, hspace_t* space)
{
    const char* ptr;

    if (data->internal_restart_req) {
        // Ignore our own requests to re-initialize.
        return 0;
    }

    ptr = hcfg_get(search_cfg, CFGKEY_GROUP_LIST);
    if (!ptr) {
        search_error(CFGKEY_GROUP_LIST
                     " configuration variable must be defined");
        return -1;
    }

    // The maximum group size is the number of input ranges.
    data->cap_max = space->len;

    data->locked_val = calloc(data->cap_max, sizeof(*data->locked_val));
    if (!data->locked_val) {
        search_error("Could not allocate memory for locked value list");
        return -1;
    }

    data->hint_val = calloc(data->cap_max, sizeof(*data->hint_val));
    if (!data->hint_val) {
        search_error("Could not allocate memory for hint value list");
        return -1;
    }

    if (hpoint_init(&data->best, data->cap_max) != 0) {
        search_error("Could not initialize best point container");
        return -1;
    }

    data->glist = NULL;
    data->glist_len = data->glist_cap = 0;
    return parse_group(data, ptr);
}

int group_setcfg(hplugin_data_t* data, const char* key, const char* val)
{
    int retval = 0;

    if (strcmp(key, CFGKEY_CONVERGED) == 0 && val && *val == '1') {
        search_best(&data->best);

        // Update locked values with converged group.
        for (int i = 0; i < data->glist[data->glist_curr].idx_len; ++i) {
            int idx = data->glist[data->glist_curr].idx[i];
            data->locked_val[idx] = data->best.term[idx];
        }

        if (++data->glist_curr < data->glist_len) {
            data->internal_restart_req = 1;
            retval = search_restart();
            data->internal_restart_req = 0;
        }
    }
    return retval;
}

int group_generate(hplugin_data_t* data, hflow_t* flow, htrial_t* trial)
{
    int ptlen = data->cap_max * sizeof(*data->hint_val);

    if (data->glist_curr < data->glist_len) {
        // Initialize locked values, if needed.
        for (int i = 0; i < data->cap_max; ++i) {
            if (data->locked_val[i].type == HVAL_UNKNOWN)
                data->locked_val[i] = trial->point.term[i];
        }

        // Base hint values on locked values.
        memcpy(data->hint_val, data->locked_val, ptlen);

        // Allow current group values from the trial point to pass through.
        for (int i = 0; i < data->glist[data->glist_curr].idx_len; ++i) {
            int idx = data->glist[data->glist_curr].idx[i];
            data->hint_val[idx] = trial->point.term[idx];
        }

        // If any values differ from our hint, reject it.
        if (memcmp(data->hint_val, trial->point.term, ptlen) != 0) {
            if (hpoint_init(&flow->point, data->cap_max) != 0) {
                search_error("Could not initialize rejection hint point");
                return -1;
            }
            memcpy((void *) flow->point.term, data->hint_val, ptlen);

            flow->status = HFLOW_REJECT;
            return 0;
        }
    }

    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Free memory associated with this search task.
 */
int group_fini(hplugin_data_t* data)
{
    if (!data->internal_restart_req) {
        for (int i = 0; i < data->glist_len; ++i)
            free(data->glist[i].idx);
        free(data->glist);

        hpoint_fini(&data->best);
        free(data->locked_val);
        free(data->hint_val);
        free(data);
    }
    return 0;
}

/*
 * Internal helper function implementation.
 */

int parse_group(hplugin_data_t* data, const char* buf)
{
    int i, j, count, success;

    int* list = malloc(data->cap_max * sizeof(int));
    char* seen = calloc(data->cap_max, sizeof(char));

    if (!list || !seen) {
        search_error("Error allocating memory for group parsing function");
        return -1;
    }

    success = 0;
    for (i = 0; *buf != '\0'; ++i) {
        while (isspace(*buf)) ++buf;
        if (*buf != '(') break;
        ++buf;

        if (i >= data->glist_cap - 1) {
            if (array_grow(&data->glist, &data->glist_cap,
                           sizeof(*data->glist)) != 0)
            {
                search_error("Error allocating group definition list");
                goto cleanup;
            }
        }

        for (j = 0; *buf != ')'; ++j) {
            if (j >= data->cap_max) {
                search_error("Too many indexes in group specification");
                goto cleanup;
            }
            if (sscanf(buf, " %d %n", &list[j], &count) < 1) {
                search_error("Error parsing group specification");
                goto cleanup;
            }
            if (list[j] < 0 || list[j] >= data->cap_max) {
                search_error("Group specification member out of bounds");
                goto cleanup;
            }
            seen[ list[j] ] = 1;
            buf += count;
            if (*buf == ',') ++buf;
        }
        ++buf;

        // Allocate memory for the parsed index list.
        data->glist[i].idx = malloc(j * sizeof(int));
        if (!data->glist[i].idx) {
            search_error("Error allocating memory for group index list");
            goto cleanup;
        }
        // Copy index from list into glist.
        memcpy(data->glist[i].idx, list, j * sizeof(int));
        data->glist[i].idx_len = j;

        while (isspace(*buf) || *buf == ',') ++buf;
    }
    if (*buf != '\0') {
        search_error("Error parsing group specification");
        goto cleanup;
    }
    data->glist_len = i;

    // Produce a final group of all unseen input indexes.
    for (i = 0, j = 0; j < data->cap_max; ++j) {
        if (!seen[j])
            list[i++] = j;
    }

    // Only add the final group if necessary.
    if (i) {
        data->glist[data->glist_len].idx = list;
        data->glist[data->glist_len].idx_len = i;
        ++data->glist_len;
    }
    else {
        free(list);
    }
    success = 1;

  cleanup:
    free(seen);

    if (!success) {
        free(list);
        return -1;
    }

    return 0;
}
