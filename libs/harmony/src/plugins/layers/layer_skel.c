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
 * \page refname Layer Name (sharedobject.so)
 *
 * Documentation written in this comment block will automatically be
 * processed by Doxygen and combined into the user document.
 */

/*
 * All Active Harmony layer plug-in implementations must begin by
 * including this header.
 */
#include "hlayer.h"

#include "session-core.h"
#include "hspace.h"
#include "hpoint.h"
#include "hcfg.h"

/*
 * Name used to identify this plugin.
 * All Harmony plugins must define this variable.
 */
const char hplugin_name[] = "<name>";

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { "VARNAME", "Default Value", "Text description." },
    { NULL }
};

/*
 * Structure to hold all data needed by an individual search instance.
 *
 * To support multiple parallel search instances, no global variables
 * should be defined or used in this plug-in layer.  They should
 * instead be defined as a part of this structure.
 */
struct hplugin_data {
    // Per-search variables go here.
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
hplugin_data_t* <name>_alloc(void)
{
    return NULL;
}

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
int <name>_init(hplugin_data_t* data, hspace_t* space)
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
int <name>_join(hplugin_data_t* data, const char* client)
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
int <name>_setcfg(hplugin_data_t* data, const char* key, const char* val)
{
    return 0;
}

/*
 * Invoked after search driver generates a candidate point, but before
 * that point is returned to the client.
 *
 * Params:
 *   data  - Search instance private data.
 *   flow  - Controls how the plug-in manager will process the point.
 *   trial - Candidate point, value, and other auxiliary trial information.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, this routine should return 0, and the contents of the
 * "flow" variable should be set appropriately.
 */
int <name>_generate(hplugin_data_t* data, hflow_t* flow, htrial_t* trial)
{
    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Invoked after the client reports a performance value, but before it
 * is processed by the search strategy.
 *
 * Params:
 *   data  - Search instance private data.
 *   flow  - Controls how the plug-in manager will process the point.
 *   trial - Candidate point, value, and other auxiliary trial information.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, this routine should return 0, and the contents of the
 * "flow" variable should be set appropriately.
 */
int <name>_analyze(hplugin_data_t* data, hflow_t* flow, htrial_t* trial)
{
    flow->status = HFLOW_ACCEPT;
    return 0;
}

/*
 * Invoked after the tuning session completes.
 *
 * Params:
 *   data - Search instance private data.
 *
 * Upon error, this function should call search_error() with a
 * human-readable string explaining the problem and return -1.
 * Otherwise, returning 0 indicates success.
 */
int <name>_fini(hplugin_data_t* data)
{
    return 0;
}
