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
 * \page omega Omega Constraint (constraint.so)
 *
 * Active Harmony allows for basic bounds on tuning variables during
 * session specification, where each tuning variable is bounded
 * individually.  However, some target applications require tuning
 * variables that are dependent upon one another, reducing the number
 * of valid parameters from the strict Cartesian product.
 *
 * This processing layer allows for the specification of such variable
 * dependencies through algebraic statements.  For example, consider a
 * target application with tuning variables `x` and `y`.  If `x` may
 * never be greater than `y`, one could use the following statement:
 *
 *     x < y
 *
 * Also, if the sum of `x` and `y` must remain under a certain
 * constant, one could use the following statement:
 *
 *     x + y = 10
 *
 * If multiple constraint statements are specified, the logical
 * conjunction of the set is applied to the
 * [Search Space](\ref intro_space).
 *
 * \note This processing layer requires the Omega Calculator, which is
 * available at:<br>
 * <https://github.com/davewathaverford/the-omega-project/>.
 *
 * \note Some search strategies provide a `REJECT_METHOD`
 * configuration variable that can be used to specify how to deal with
 * rejected points.  This can have great affect on the productivity of
 * a tuning session.
 */

#include "hlayer.h"
#include "session-core.h"
#include "hspace.h"
#include "hpoint.h"
#include "hutil.h"
#include "hsockutil.h"
#include "hcfg.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/wait.h>
#include <signal.h>

/*
 * Name used to identify this plugin layer.
 * All Harmony plugin layers must define this variable.
 */
const char hplugin_name[] = "constraint";

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_OC_BIN, "oc",
      "Location of the Omega Calculator binary. The PATH environment "
      "variable will be searched if not found initially." },
    { CFGKEY_OC_CONSTRAINTS, NULL,
      "Constraint statements to be used during this session. "
      "This variable has precedence over " CFGKEY_OC_FILE "." },
    { CFGKEY_OC_FILE, NULL,
      "If the " CFGKEY_OC_CONSTRAINTS " variable is not specified, "
      "constraints will be loaded from this file." },
    { CFGKEY_OC_QUIET, "False",
      "Bounds suggestion and rejection messages can be suppressed "
      "by setting this variable to true." },
    { NULL }
};

#define MAX_CMD_LEN  4096
#define MAX_TEXT_LEN 1024

/*
 * Structure to hold all data needed by an individual search instance.
 *
 * To support multiple parallel search instances, no global variables
 * should be defined or used in this plug-in layer.  They should
 * instead be defined as a part of this structure.
 */
struct hplugin_data {
    hspace_t local_space;

    const char* omega_bin;
    char constraints[MAX_TEXT_LEN];

    char vars_text[MAX_TEXT_LEN];
    char bounds_text[MAX_TEXT_LEN];
    char user_text[MAX_TEXT_LEN];
    char point_text[MAX_TEXT_LEN];

    int quiet;
};

/*
 * Internal helper function prototypes.
 */
static int   strategy_cfg(hplugin_data_t* data, hspace_t* space);
static int   build_vars_text(hplugin_data_t* data, hspace_t* space);
static int   build_bounds_text(hplugin_data_t* data, hspace_t* space);
static int   build_user_text(hplugin_data_t* data);
static int   build_point_text(hplugin_data_t* data, hpoint_t* point);
static int   update_bounds(hplugin_data_t* data, hspace_t* space);
static int   check_validity(hplugin_data_t* data, hpoint_t* point);
static char* call_omega_calc(hplugin_data_t* data, const char* cmd);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* constraint_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search task.
 */
int constraint_init(hplugin_data_t* data, hspace_t* space)
{
    // Make a copy of the search space.
    hspace_copy(&data->local_space, space);

    // Initialize layer variables.
    if (strategy_cfg(data, space) != 0)
        return -1;

    if (build_vars_text(data, space) != 0)
        return -1;

    if (build_bounds_text(data, space) != 0)
        return -1;

    if (build_user_text(data) != 0)
        return -1;

    // Calculate the range for each tuning variable, given the constraints.
    if (update_bounds(data, space) != 0)
        return -1;

    return 0;
}

int constraint_generate(hplugin_data_t* data, hflow_t* flow, hpoint_t* point)
{
    flow->status = HFLOW_ACCEPT;
    if (!check_validity(data, point)) {
        flow->status = HFLOW_REJECT;
        flow->point.id = -1;

        if (!data->quiet) {
            fprintf(stderr, "Rejecting point: {");
            for (int i = 0; i < point->len; ++i) {
                const hval_t* val = &point->term[i];

                switch (data->local_space.dim[i].type) {
                case HVAL_INT:  fprintf(stderr, "%ld", val->value.i); break;
                case HVAL_REAL: fprintf(stderr, "%g", val->value.r); break;
                case HVAL_STR:  fprintf(stderr, "%s", val->value.s); break;
                default:        fprintf(stderr, "<INV>");
                }
                if (i < point->len - 1)
                    fprintf(stderr, ", ");
            }
            fprintf(stderr, "}\n");
        }
    }
    return 0;
}

/*
 * Free memory associated with this search task.
 */
int constraint_fini(hplugin_data_t* data)
{
    hspace_fini(&data->local_space);
    free(data);
    return 0;
}

/*
 * Internal helper function implementation.
 */

int strategy_cfg(hplugin_data_t* data, hspace_t* space)
{
    const char* cfgval;

    data->omega_bin = hcfg_get(search_cfg, CFGKEY_OC_BIN);
    if (!file_exists(data->omega_bin)) {
        data->omega_bin = search_path(data->omega_bin);
        if (!data->omega_bin) {
            search_error("Could not find Omega Calculator executable. "
                         "Use " CFGKEY_OC_BIN " to specify its location");
            return -1;
        }
    }

    data->quiet = hcfg_bool(search_cfg, CFGKEY_OC_QUIET);

    cfgval = hcfg_get(search_cfg, CFGKEY_OC_CONSTRAINTS);
    if (cfgval) {
        if (strlen(cfgval) >= sizeof(data->constraints)) {
            search_error("Constraint string too long");
            return -1;
        }
        strncpy(data->constraints, cfgval, sizeof(data->constraints));
    }
    else {
        size_t retval;
        FILE* fp;

        cfgval = hcfg_get(search_cfg, CFGKEY_OC_FILE);
        if (!cfgval) {
            search_error("No constraints specified.  Either "
                         CFGKEY_OC_CONSTRAINTS " or " CFGKEY_OC_FILE
                         " must be defined");
            return -1;
        }

        fp = fopen(cfgval, "r");
        if (!fp) {
            search_error("Could not open constraint file");
            return -1;
        }

        retval = fread(data->constraints, sizeof(char),
                       sizeof(data->constraints), fp);
        data->constraints[sizeof(data->constraints) - 1] = '\0';

        if (retval >= sizeof(data->constraints)) {
            search_error("Constraint file too large");
            return -1;
        }

        if (fclose(fp) != 0) {
            search_error("Could not close constraint file");
            return -1;
        }
    }
    return 0;
}

int build_vars_text(hplugin_data_t* data, hspace_t* space)
{
    data->vars_text[0] = '\0';
    for (int i = 0; i < space->len; ++i) {
        strcat(data->vars_text, space->dim[i].name);
        if (i < space->len - 1)
            strcat(data->vars_text, ", ");
    }
    return 0;
}

int build_bounds_text(hplugin_data_t* data, hspace_t* space)
{
    char* ptr = data->bounds_text;
    char* end = data->bounds_text + sizeof(data->bounds_text);

    data->bounds_text[0] = '\0';
    for (int i = 0; i < space->len; ++i) {
        hrange_t* range = &space->dim[i];

        // Fetch min and max according to variable type.
        switch (range->type) {
        case HVAL_INT:
            ptr += snprintf(ptr, end - ptr, "%ld <= %s <= %ld",
                            range->bounds.i.min, range->name,
                            range->bounds.i.max);
            break;

        case HVAL_REAL:
            ptr += snprintf(ptr, end - ptr, "%lf <= %s <= %lf",
                            range->bounds.r.min, range->name,
                            range->bounds.r.max);
            break;

        case HVAL_STR:
        default:
            return -1;
        }

        if (i < space->len - 1)
            ptr += snprintf(ptr, end - ptr, " && ");
    }
    return 0;
}

int build_user_text(hplugin_data_t* data)
{
    const char* ptr = data->constraints;
    const char* end;
    int len = 0;

    while (*ptr) {
        while (isspace(*ptr)) ++ptr;
        end = ptr + strcspn(ptr, "\n");

        len += end - ptr;
        if (len < sizeof(data->user_text))
            strncat(data->user_text, ptr, end - ptr);

        while (isspace(*end)) ++end;
        if (*end) {
            len += 4;
            if (len < sizeof(data->user_text))
                strcat(data->user_text, " && ");
        }

        if (len >= sizeof(data->user_text)) {
            search_error("User constraint string overflow");
            return -1;
        }
        ptr = end;
    }
    return 0;
}

int build_point_text(hplugin_data_t* data, hpoint_t* point)
{
    char* ptr = data->point_text;
    char* end = data->point_text + sizeof(data->point_text);

    data->point_text[0] = '\0';
    for (int i = 0; i < point->len; ++i) {
        // Fetch min and max according to variable type.
        switch (data->local_space.dim[i].type) {
        case HVAL_INT:
            ptr += snprintf(ptr, end - ptr, "%s = %ld",
                            data->local_space.dim[i].name,
                            point->term[i].value.i);
            break;

        case HVAL_REAL:
            ptr += snprintf(ptr, end - ptr, "%s = %f",
                            data->local_space.dim[i].name,
                            point->term[i].value.r);
            break;

        case HVAL_STR:
        default:
            return -1;
        }

        if (i < point->len - 1)
            ptr += snprintf(ptr, end - ptr, " && ");
    }
    return 0;
}

/*
 * XXX - We don't actually update the session search space just yet,
 * resulting in correct, but less optimal point generation.
 */
int update_bounds(hplugin_data_t* data, hspace_t* space)
{
    char cmd[MAX_CMD_LEN];
    char* retstr;
    int retval = 0;

    for (int i = 0; i < data->local_space.len; ++i) {
        hrange_t* range = &data->local_space.dim[i];

        // Write the domain text with variable constraints to the file.
        snprintf(cmd, sizeof(cmd),
                 "symbolic %s;\n"
                 "D:={[%s]: %s && %s};\n"
                 "range D;",
                 range->name, data->vars_text,
                 data->bounds_text, data->user_text);

        // Call omega calculator.
        retstr = call_omega_calc(data, cmd);
        if (!retstr)
            return -1;

        // Parse the result.
        while (*retstr) {
            if (strncmp(">>>", retstr, 3) == 0)
                goto nextline; // Skip echo lines.

            switch (range->type) {
            case HVAL_INT:
                retval = sscanf(retstr, "{ %ld <= %*s <= %ld }\n",
                                &range->bounds.i.min, &range->bounds.i.max);
                break;

            case HVAL_REAL:
                retval = sscanf(retstr, "{ %lf <= %*s <= %lf }\n",
                                &range->bounds.r.min, &range->bounds.r.max);
                break;

            case HVAL_STR:
            default:
                search_error("Constraint layer cannot handle string types");
                return -1;
            }

            if (retval != 2)
                fprintf(stderr, "Unexpected Omega Calculator output: %.*s\n",
                        (int) strcspn(retstr, "\n"), retstr);
            break;

          nextline:
            retstr += strcspn(retstr, "\n");
            if (*retstr)
                ++retstr;
        }
        if (retval != 2) {
            search_error("Error parsing Omega Calculator output");
            return -1;
        }
    }

    if (!data->quiet) {
        if (!hspace_equal(space, &data->local_space)) {
            fprintf(stderr, "For the given constraints, we suggest re-running"
                    " the session with these bounds:\n");

            for (int i = 0; i < data->local_space.len; ++i) {
                hrange_t* range = &data->local_space.dim[i];

                switch (range->type) {
                case HVAL_INT:
                    fprintf(stderr, "INT %s = {%ld, %ld, %ld}\n",
                            range->name,
                            range->bounds.i.min,
                            range->bounds.i.max,
                            range->bounds.i.step);
                    break;
                case HVAL_REAL:
                    fprintf(stderr, "REAL %s = {%lf, %lf, %lf}\n",
                            range->name,
                            range->bounds.r.min,
                            range->bounds.r.max,
                            range->bounds.r.step);
                    break;
                case HVAL_STR:
                default:
                    return -1;
                }
            }
        }
    }
    return 0;
}

int check_validity(hplugin_data_t* data, hpoint_t* point)
{
    char cmd[MAX_CMD_LEN];
    char* retstr;
    int retval = 0;

    if (build_point_text(data, point) != 0)
        return -1;

    snprintf(cmd, sizeof(cmd),
             "D:={[%s]: %s && %s && %s};\n"
             "D;",
             data->vars_text, data->bounds_text,
             data->user_text, data->point_text);

    // Call omega test.
    retstr = call_omega_calc(data, cmd);
    if (!retstr)
        return -1;

    // Parse result.
    while (*retstr) {
        char c;
        if (strncmp(">>>", retstr, 3) == 0)
            goto nextline; // Skip echo lines.

        retval = sscanf(retstr, " { %*s : FALSE %c ", &c);

        if (retval == 1)
            break;

      nextline:
        retstr += strcspn(retstr, "\n");
        if (*retstr)
            ++retstr;
    }

    // If sscanf ever succeeded, retval will be 1 and the point is invalid.
    return (retval != 1);
}

char* call_omega_calc(hplugin_data_t* data, const char* cmd)
{
    static char* buf = NULL;
    static int buf_cap = 4096;
    char* child_argv[2];
    char* ptr;
    pid_t oc_pid;
    int fd, count;

    if (buf == NULL) {
        buf = malloc(sizeof(char) * buf_cap);
        if (!buf)
            return NULL;
    }

    child_argv[0] = (char*) data->omega_bin;
    child_argv[1] = NULL;
    fd = socket_launch(data->omega_bin, child_argv, &oc_pid);
    if (!fd) {
        search_error("Could not launch Omega Calculator");
        return NULL;
    }

    if (socket_write(fd, cmd, strlen(cmd)) < 0) {
        search_error("Could not send command to Omega Calculator");
        return NULL;
    }

    if (shutdown(fd, SHUT_WR) != 0) {
        search_error("Internal error: Could not shutdown write to Omega");
        return NULL;
    }

    ptr = buf;
    while ( (count = socket_read(fd, ptr, buf_cap - (buf - ptr))) > 0) {
        ptr += count;
        if (ptr == buf + buf_cap) {
            if (array_grow(&buf, &buf_cap, sizeof(char)) != 0) {
                search_error("Internal error: Could not grow buffer for"
                             " Omega Calculator input");
                return NULL;
            }
            ptr = buf + strlen(buf);
        }
        *ptr = '\0';
    }
    if (count < 0) {
        search_error("Could not read output from Omega Calculator");
        return NULL;
    }

    if (close(fd) != 0) {
        search_error("Internal error: Could not close Omega socket");
        return NULL;
    }

    if (waitpid(oc_pid, NULL, 0) != oc_pid) {
        search_error("Internal error: Could not reap Omega process");
        return NULL;
    }

    return buf;
}
