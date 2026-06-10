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

#ifndef __HPLUGIN_H__
#define __HPLUGIN_H__

#include "session-core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum hplugin_type {
    HPLUGIN_UNKNOWN = 0,
    HPLUGIN_STRATEGY,
    HPLUGIN_LAYER,

    HPLUGIN_MAX
} hplugin_type_t;

/*
 * Declare the external structure used to hold per-search instance
 * data.  This is the "struct hplugin_data" that should be defined by each
 * plug-in module.
 */
typedef struct hplugin_data hplugin_data_t;

/*
 * Plug-in instance data management function signatures.
 */
typedef hplugin_data_t* (*hook_alloc_t)(void);

/*
 * Strategy plug-in event function signatures.
 */
typedef int (*strategy_generate_t)(hplugin_data_t* data,
                                   hflow_t* flow, hpoint_t* point);
typedef int (*strategy_rejected_t)(hplugin_data_t* data,
                                   hflow_t* flow, hpoint_t* point);
typedef int (*strategy_analyze_t)(hplugin_data_t* data, htrial_t* trial);
typedef int (*strategy_best_t)(hplugin_data_t* data, hpoint_t* point);

/*
 * Layer plug-in event function signatures.
 */
typedef int (*layer_generate_t)(hplugin_data_t* data,
                                hflow_t* flow, htrial_t* trial);
typedef int (*layer_analyze_t)(hplugin_data_t* data,
                               hflow_t* flow, htrial_t* trial);

/*
 * Generic plug-in event function signatures.
 */
typedef int (*hook_init_t)(hplugin_data_t* data, hspace_t* space);
typedef int (*hook_join_t)(hplugin_data_t* data, const char* client);
typedef int (*hook_setcfg_t)(hplugin_data_t* data,
                             const char* key, const char* val);
typedef int (*hook_fini_t)(hplugin_data_t* data);

/*
 * Harmony structure to encapsulate the state of a Harmony plug-in.
 */
typedef struct hplugin {
    hplugin_type_t type;

    // Variable addresses, loaded from plugin by name.
    const char*        type_str;
    const char*        name;
    const hcfg_info_t* keyinfo;

    // Plug-in instance data initialization and storage.
    hook_alloc_t    alloc;
    hplugin_data_t* data;

    // Strategy event function addresses.
    struct hook_strategy {
        strategy_generate_t generate;
        strategy_rejected_t rejected;
        strategy_analyze_t  analyze;
        strategy_best_t     best;
    } strategy;

    // Layer event function addresses.
    struct hook_layer {
        layer_generate_t generate;
        layer_analyze_t  analyze;
    } layer;

    // Common event function addresses.
    hook_init_t   init;
    hook_join_t   join;
    hook_setcfg_t setcfg;
    hook_fini_t   fini;

    // Dynamic library handle.
    void* handle;
} hplugin_t;

#define HPLUGIN_INITIALIZER {HPLUGIN_UNKNOWN}
extern const hplugin_t hplugin_zero;

/*
 * Base structure management interface.
 */
int hplugin_open(hplugin_t* plugin, const char* filename, const char** errptr);
int hplugin_close(hplugin_t* plugin, const char** errptr);

/*
 * Event function calling interface.
 */
int hplugin_analyze(hplugin_t* plugin, hflow_t* flow, htrial_t* trial);
int hplugin_best(hplugin_t* plugin, hpoint_t* point);
int hplugin_generate(hplugin_t* plugin, hflow_t* flow, htrial_t* trial);
int hplugin_rejected(hplugin_t* plugin, hflow_t* flow, htrial_t* trial);

int hplugin_init(hplugin_t* plugin, hspace_t* space);
int hplugin_join(hplugin_t* plugin, const char* client);
int hplugin_setcfg(hplugin_t* plugin, const char* key, const char* val);
int hplugin_fini(hplugin_t* plugin);

#ifdef __cplusplus
}
#endif

#endif
