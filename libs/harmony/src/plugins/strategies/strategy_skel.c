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

/*
 * All Active Harmony strategy plug-in implementations must begin by
 * including this header.
 */
#include "hstrategy.h"

/*
 * Configuration variables used in this plugin.
 *
 * Any entries in this structure will automatically be loaded into the
 * search environment.
 */
const hcfg_info_t hplugin_keyinfo[] = {
};

/*
 * Structure to hold all data needed by an individual search instance.
 *
 * To support multiple parallel search instances, no global variables
 * should be defined or used in this plug-in layer.  They should
 * instead be defined as a part of this structure.
 */
struct hplugin_data {
};

/*
 * The following functions are used for plug-in setup.
 *
 * They are called once along with each search task that uses this
 * plug-in.
 */

/*
 * Allocate plug-in state data for a new search instance.  This
 * function is invoked prior to any other plug-in interface function.
 *
 * Upon success, this function should return an address the plug-in
 * can use to track an individual search instance.  Otherwise, this
 * function should return NULL.
 */
hplugin_data_t* strategy_alloc(void)
{
    return NULL;
}

/*
 * The following functions are required.
 *
 * Active Harmony will not recognize shared objects as search
 * strategies unless these functions exist.
 */

/*
 * Generate a new candidate configuration point.
 *
 * Params:
 *   data  - Search instance private data.
 *   flow  - Inform plug-in manager of search state.
 *   point - New candidate point to be sent to a client.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, this routine should return 0, and the contents of the
 * "flow" variable should be set appropriately.
 */
int strategy_generate(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Regenerate a rejected candidate configuration point.
 *
 * Params:
 *   data  - Search instance private data.
 *   flow  - Inform plug-in manager of search state. Also may contain a
 *           search hint from the rejecting layer.
 *   point - Rejected candidate point to regenerate.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, this routine should return 0, and the contents of the
 * "flow" variable should be set appropriately.
 */
int strategy_rejected(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Analyze the observed performance associated with a configuration point.
 *
 * Params:
 *   data  - Search instance private data.
 *   trial - Observed point/performance pair.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, returning 0 indicates success.
 */
int strategy_analyze(hplugin_data_t* data, htrial_t* trial)
{
    return 0;
}

/*
 * The following functions are optional.
 *
 * They will be invoked at the appropriate time if and only if they
 * exist.
 */

/*
 * Initialize (or re-initialize) data for this search instance.
 *
 * Param:
 *   data  - Search instance private data.
 *   space - Details of the parameter space (dimensions, bounds, etc.).
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, returning 0 indicates success.
 */
int strategy_init(hplugin_data_t* data, hspace_t* space)
{
    return 0;
}

/*
 * Invoked when a client joins the tuning session.
 *
 * Params:
 *   data   - Search instance private data.
 *   client - Uniquely identifying string for the new client.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, returning 0 indicates success.
 */
int strategy_join(hplugin_data_t* data, const char* client)
{
    return 0;
}

/*
 * Invoked when a client writes to the configuration system.
 *
 * Params:
 *   data - Search instance private data.
 *   key  - Configuration key to be modified.
 *   val  - New value for configuration key.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, returning 0 indicates success.
 */
int strategy_setcfg(hplugin_data_t* data, const char* key, const char* val)
{
    return 0;
}

/*
 * Invoked on session exit.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, returning 0 indicates success.
 */
int strategy_fini(hplugin_data_t* data)
{
    return 0;
}
