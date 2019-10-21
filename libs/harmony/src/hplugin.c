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

#include "hplugin.h"
#include "hutil.h"

#include <stdlib.h> // For NULL and free().
#include <string.h> // For memset().
#include <dlfcn.h>  // For dlopen(), dlsym(), dlerror(), and dlclose().

const hplugin_t hplugin_zero = HPLUGIN_INITIALIZER;

/*
 * Internal helper function prototypes.
 */
static int load_hooks(hplugin_t* plugin);
static int verify_hooks(hplugin_t* plugin, const char** errptr);
static int verify_type(hplugin_t* plugin, const char** errptr);

/*
 * Base structure management implementation.
 */

int hplugin_open(hplugin_t* plugin, const char* filename, const char** errptr)
{
    const char* errstr;

    // Tabula rasa.
    memset(plugin, 0, sizeof(*plugin));

    plugin->handle = dlopen(filename, RTLD_LAZY | RTLD_LOCAL);
    if (!plugin->handle) {
        errstr = dlerror();
        goto error;
    }

    // Find plug-in variables.
    plugin->type_str = dlsym(plugin->handle, "hplugin_type");
    plugin->name     = dlsym(plugin->handle, "hplugin_name");
    plugin->keyinfo  = dlsym(plugin->handle, "hplugin_keyinfo");

    if (verify_type(plugin, &errstr) != 0)
        return -1;

    if (load_hooks(plugin) != 0) {
        errstr = "Could not allocate memory for plug-in symbol";
        goto error;
    }

    if (verify_hooks(plugin, &errstr) != 0)
        return -1;

    if (plugin->alloc) {
        plugin->data = plugin->alloc();
        if (!plugin->data) {
            errstr = "Could not allocate private data for plug-in";
            goto error;
        }
    }

    return 0;

  error:
    if (plugin->handle)
        dlclose(plugin->handle);

    if (errptr)
        *errptr = errstr;
    return -1;
}

int hplugin_close(hplugin_t* plugin, const char** errptr)
{
    const char* errstr;

    if (dlclose(plugin->handle) != 0) {
        if (errptr) {
            errstr = dlerror();
            goto error;
        }
    }

    return 0;

  error:
    if (errptr)
        *errptr = errstr;
    return -1;
}

/*
 * Event function calling interface.
 */

int hplugin_analyze(hplugin_t* plugin, hflow_t* flow, htrial_t* trial)
{
    if (plugin->type == HPLUGIN_STRATEGY) {
        return plugin->strategy.analyze(plugin->data, trial);
    }
    else if (plugin->type == HPLUGIN_LAYER) {
        if (plugin->layer.analyze)
            return plugin->layer.analyze(plugin->data, flow, trial);
    }
    return 0;
}

int hplugin_best(hplugin_t* plugin, hpoint_t* point)
{
    if (plugin->type != HPLUGIN_STRATEGY)
        return -1;

    if (plugin->strategy.best)
        return plugin->strategy.best(plugin->data, point);
    else
        return 0;
}

int hplugin_generate(hplugin_t* plugin, hflow_t* flow, htrial_t* trial)
{
    if (plugin->type == HPLUGIN_STRATEGY) {
        hpoint_t* point = (hpoint_t*) &trial->point;
        return plugin->strategy.generate(plugin->data, flow, point);
    }
    else if (plugin->type == HPLUGIN_LAYER) {
        if (plugin->layer.generate)
            return plugin->layer.generate(plugin->data, flow, trial);
    }
    return 0;
}

int hplugin_rejected(hplugin_t* plugin, hflow_t* flow, htrial_t* trial)
{
    if (plugin->type != HPLUGIN_STRATEGY)
        return -1;

    hpoint_t* point = (hpoint_t*) &trial->point;
    return plugin->strategy.rejected(plugin->data, flow, point);
}

int hplugin_init(hplugin_t* plugin, hspace_t* space)
{
    if (plugin->init)
        return plugin->init(plugin->data, space);
    else
        return 0;
}

int hplugin_join(hplugin_t* plugin, const char* client)
{
    if (plugin->join)
        return plugin->join(plugin->data, client);
    else
        return 0;
}

int hplugin_setcfg(hplugin_t* plugin, const char* key, const char* val)
{
    if (plugin->setcfg)
        return plugin->setcfg(plugin->data, key, val);
    else
        return 0;
}

int hplugin_fini(hplugin_t* plugin)
{
    if (plugin->fini)
        return plugin->fini(plugin->data);
    else
        return 0;
}

/*
 * Internal helper function implementation.
 */

/*
 * ISO C forbids conversion of object pointers to function pointers,
 * making it difficult to use dlsym() for functions.  We get around
 * this by first casting to a word-length integer.  (ILP32/LP64
 * compilers assumed).
 */
#define dlfptr(x, y) ((void*) (void (*)(void))(long)(dlsym((x), (y))))

int load_hooks(hplugin_t* plugin)
{
    const char* prefix = plugin->name;
    void*       handle = plugin->handle;

    char* buf = NULL;
    int   len = 0;
    int   retval = 0;

    if (plugin->type == HPLUGIN_STRATEGY) {
        // Find strategy specific event functions.
        plugin->strategy.generate =
            (strategy_generate_t) dlfptr(handle, "strategy_generate");
        plugin->strategy.rejected =
            (strategy_rejected_t) dlfptr(handle, "strategy_rejected");
        plugin->strategy.analyze =
            (strategy_analyze_t) dlfptr(handle, "strategy_analyze");
        plugin->strategy.best =
            (strategy_best_t) dlfptr(handle, "strategy_best");
        prefix = "strategy";
    }
    else if (plugin->type == HPLUGIN_LAYER) {
        // Find layer specific event functions.
        if (snprintf_grow(&buf, &len, "%s_generate", prefix) < 0) goto error;
        plugin->layer.generate = (layer_generate_t) dlfptr(handle, buf);

        if (snprintf_grow(&buf, &len, "%s_analyze", prefix) < 0) goto error;
        plugin->layer.analyze = (layer_analyze_t) dlfptr(handle, buf);
    }

    // Load optional plug-in hooks.
    if (snprintf_grow(&buf, &len, "%s_alloc", prefix) < 0) goto error;
    plugin->alloc = (hook_alloc_t) dlfptr(handle, buf);

    if (snprintf_grow(&buf, &len, "%s_init", prefix) < 0) goto error;
    plugin->init = (hook_init_t) dlfptr(handle, buf);

    if (snprintf_grow(&buf, &len, "%s_join", prefix) < 0) goto error;
    plugin->join = (hook_join_t) dlfptr(handle, buf);

    if (snprintf_grow(&buf, &len, "%s_setcfg", prefix) < 0) goto error;
    plugin->setcfg = (hook_setcfg_t) dlfptr(handle, buf);

    if (snprintf_grow(&buf, &len, "%s_fini", prefix) < 0) goto error;
    plugin->fini = (hook_fini_t) dlfptr(handle, buf);

    goto cleanup;

  error:
    retval = -1;

  cleanup:
    free(buf);
    return retval;
}

int verify_hooks(hplugin_t* plugin, const char** errptr)
{
    if (plugin->type == HPLUGIN_STRATEGY) {
        if (!plugin->strategy.generate) {
            *errptr = "Strategy plug-in missing symbol: strategy_generate";
            return -1;
        }

        if (!plugin->strategy.analyze) {
            *errptr = "Strategy plug-in missing symbol: strategy_analyze";
            return -1;
        }

        if (!plugin->strategy.rejected) {
            *errptr = "Strategy plug-in missing symbol: strategy_rejected";
            return -1;
        }

        if (!plugin->strategy.best) {
            *errptr = "Strategy plug-in missing symbol: strategy_best";
            return -1;
        }
    }
    else if (plugin->type == HPLUGIN_LAYER) {
        if (!plugin->layer.generate &&
            !plugin->layer.analyze &&
            !plugin->init &&
            !plugin->join &&
            !plugin->setcfg)
        {
            *errptr = "Layer plug-in does not define any event functions";
            return -1;
        }
    }
    else {
        *errptr = "Improper Active Harmony plug-in type";
        return -1;
    }

    return 0;
}

int verify_type(hplugin_t* plugin, const char** errptr)
{
    if (!plugin->type_str) {
        *errptr = "Library is not an Active Harmony plug-in";
        return -1;
    }
    else if (strcmp(plugin->type_str, "strategy") == 0) {
        plugin->type = HPLUGIN_STRATEGY;
    }
    else if (strcmp(plugin->type_str, "layer") == 0) {
        if (!plugin->name) {
            *errptr = "Layer plug-in does not define its name";
            return -1;
        }
        plugin->type = HPLUGIN_LAYER;
    }
    else {
        *errptr = "Unknown Active Harmony plug-in type";
        return -1;
    }

    return 0;
}
