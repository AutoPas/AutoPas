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
#define _XOPEN_SOURCE 600 // Needed for gethostname()

#include "hclient.h"
#include "hcfg.h"
#include "hspace.h"
#include "hpoint.h"
#include "hperf.h"
#include "hval.h"
#include "hmesg.h"
#include "hutil.h"
#include "hsockutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

typedef enum harmony_state_t {
    HARMONY_STATE_UNKNOWN,
    HARMONY_STATE_DETACHED,
    HARMONY_STATE_ATTACHED,
    HARMONY_STATE_READY,
    HARMONY_STATE_TESTING,

    HARMONY_STATE_MAX
} harmony_state_t;

/*
 * Harmony server session descriptor.
 *
 * A single session may host for multiple searches.
 */
struct hdesc {
    int     socket;
    hmesg_t mesg;

    char*  id;
    hcfg_t cfg; // Overriding configuration directives, such as those
                // taken from the environment and the command line.

    htask_t** tlist;
    int       tlist_len;
    int       tlist_cap;
    char* buf;
    int   buflen;
};

/*
 * Harmony search definition descriptor.
 *
 * Stores information related to how a search is defined, such as
 * input space parameters and output objectives.
 */
struct hdef {
    hspace_t space;
    hcfg_t   cfg;

    void** varloc;
    int    varloc_cap;
};

/*
 * Harmony search task descriptor.
 *
 * Stores information related to a search that is currently running,
 * such the current point to be tested.
 */
struct htask {
    harmony_state_t state;

    hspace_t space;
    hcfg_t   cfg;

    hdesc_t* hdesc;
    hmesg_t* mesg;
    int      dest;

    void** varloc; // Location where input parameters should be stored.
    int    varloc_cap;

    hpoint_t test; // Current point to test retrieved from the search task.
    hpoint_t best; // Best performing point across all search task clients.
    hperf_t  perf; // Performance of the current testing point.

    hpoint_t* curr; // This will point at either "test" or "best."
};

/*
 * Internal helper function prototypes.
 */
static htask_t* alloc_task(hdesc_t* hdesc);
static int      init_task(htask_t* htask, const hspace_t* space,
                          const hcfg_t* cfg);
static void     free_task(htask_t* htask);

static int   extend_perf(htask_t* htask);
static int   extend_varloc(void*** varloc, int* varloc_cap, hspace_t* space);
static int   find_var(htask_t* htask, const char* name);
static char* generate_id(hdesc_t* hdesc, int suffix);
static int   send_request(htask_t* htask, hmesg_type msg_type);
static int   set_varloc(htask_t* htask, const char* name, void* ptr);
static int   write_values(htask_t* htask);

/*
 * Global static variables.
 */
static const char* ah_errstr;
static char*       ah_errbuf;
static int         ah_errbuflen;

/*
 * Harmony descriptor implementation.
 */

/**
 * \defgroup api_desc Harmony Descriptor Management Functions
 *
 * A Harmony descriptor is an opaque structure that stores state
 * associated with a client's connection to a tuning session.  A
 * tuning session hosts one or more tuning search tasks, and is the
 * only intermediary by which clients may communicate with search
 * tasks.  The functions within this section allow the user to create
 * and manage the resources controlled by the descriptor.
 *
 * @{
 */

/**
 * \brief Allocate and initialize a new Harmony descriptor.
 *
 * All client API functions require a valid Harmony descriptor to
 * function correctly.  The descriptor is designed to be used as an
 * opaque type and no guarantees are made about the contents of its
 * structure.
 *
 * Heap memory is allocated for the descriptor, so be sure to call
 * ah_free() when it is no longer needed.
 *
 * \return Returns a new Harmony descriptor upon success, and `NULL`
 *         otherwise.
 */
hdesc_t* ah_alloc(void)
{
    hdesc_t* hdesc = calloc(1, sizeof(*hdesc));
    if (!hdesc) {
        ah_errstr = "Could not allocate memory for Harmony descriptor";
        return NULL;
    }

    hdesc->id = stralloc( generate_id(hdesc, 0) );
    if (!hdesc->id) {
        ah_errstr = "Could not allocate memory for client ID string";
        goto error;
    }

    if (hcfg_loadenv(&hdesc->cfg) != 0) {
        ah_errstr = "Could not copy environment into descriptor configuration";
        goto error;
    }

    hdesc->socket = -1;
    return hdesc;

  error:
    hcfg_fini(&hdesc->cfg);
    free(hdesc->id);
    free(hdesc);
    return NULL;
}

/**
 * \brief Find configuration directives in the command-line arguments.
 *
 * This function iterates through the command-line arguments (as
 * passed to the main() function) to search for Harmony configuration
 * directives in the form NAME=VAL.  Any arguments matching this
 * pattern are added to the configuration and removed from `argv`.
 *
 * Iteration through the command-line arguments stops prematurely if a
 * double-dash ("--") token is found.  All non-directive arguments are
 * shifted to the front of `argv`.  Finally, `argc` is updated to
 * reflect the remaining number of program arguments.
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 * \param argc  Address of the `argc` variable.
 * \param argv  Address of the command-line argument array
 *              (i.e., the value of `argv`).
 *
 * \return Returns the number of directives successfully taken from
 *         `argv`, or -1 on error.
 */
int ah_args(hdesc_t* hdesc, int* argc, char** argv)
{
    int skip = 0;
    int tail = 0;

    for (int i = 0; i < *argc; ++i) {
        if (strcmp(argv[i], "--") == 0)
            skip = 1;

        if (skip || hcfg_parse(&hdesc->cfg, argv[i], NULL) != 1)
            argv[tail++] = argv[i];
    }

    int removed = *argc - tail;
    *argc = tail;
    return removed;
}

/**
 * \brief Assign an identifying string to this client.
 *
 * Set the client identification string.  All messages sent to the
 * tuning session will be tagged with this string, allowing the
 * framework to distinguish clients from one another.  As such, care
 * should be taken to ensure this string is unique among all clients
 * participating in a tuning session.
 *
 * If this function is not called, the client will be assigned an ID
 * based on its hostname and process ID.
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 * \param id    Unique identification string.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_id(hdesc_t* hdesc, const char* id)
{
    if (!id) {
        ah_errstr = "Invalid ID string";
        return -1;
    }

    char* newid = stralloc(id);
    if (!newid) {
        ah_errstr = strerror(errno);
        return -1;
    }

    free(hdesc->id);
    hdesc->id = newid;

    return 0;
}

/**
 * \brief Establish a connection with a Harmony tuning session.
 *
 * This function enables communication with a tuning session, Active
 * Harmony's search management process.  It must succeed prior to
 * calling any search interface functions (e.g., ah_join(),
 * ah_launch(), etc.).
 *
 * If *host* is `NULL`, its value will be taken from the environment
 * variable `HARMONY_HOST`.  If `HARMONY_HOST` is not defined, the
 * environment variable `HARMONY_HOME` will be used to initiate a
 * private tuning session, which will be available only to the local
 * process.
 *
 * If *port* is 0, its value will be taken from the environment
 * variable `HARMONY_PORT`, if defined.  Otherwise, its value will be
 * taken from the src/defaults.h file.
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 * \param host  Host of the Harmony server (or `NULL`).
 * \param port  Port of the Harmony server.
 *
 * \return Returns 0 on success.  Otherwise, -1 is returned and a
 *         string indicating the nature of the failure may be
 *         retrieved from ah_error().
 */
int ah_connect(hdesc_t* hdesc, const char* host, int port)
{
    hcfg_t connect_cfg;
    int retval = 0;

    // Initialize a configuration with defaults.
    if (hcfg_init(&connect_cfg) != 0) {
        ah_errstr = "Could not seed temporary hcfg_t with defaults";
        goto error;
    }

    // Merge in overriding configuration directives.
    if (hcfg_merge(&connect_cfg, &hdesc->cfg) != 0) {
        ah_errstr = "Could not add overriding directives to temporary hcfg_t";
        goto error;
    }

    // Sanity check input.
    if (hdesc->socket != -1) {
        ah_errstr = "Descriptor already attached to another session";
        goto error;
    }

    if (!host)
        host = hcfg_get(&connect_cfg, CFGKEY_HARMONY_HOST);

    if (port == 0)
        port = hcfg_int(&connect_cfg, CFGKEY_HARMONY_PORT);

    if (!host) {
        char* path;
        const char* home;

        // Find the Active Harmony installation.
        home = hcfg_get(&connect_cfg, CFGKEY_HARMONY_HOME);
        if (!home) {
            ah_errstr = "No host or " CFGKEY_HARMONY_HOME " specified";
            goto error;
        }

        // Fork a local tuning session.
        path = sprintf_alloc("%s/libexec/" SESSION_CORE_EXECFILE, home);
        if (!path) {
            ah_errstr = "Could not allocate memory for session filename";
            goto error;
        }

        char* const child_argv[] = {path,
                                    (char*)home,
                                    NULL};
        hdesc->socket = socket_launch(path, child_argv, NULL);
        free(path);
    }
    else {
        hdesc->socket = tcp_connect(host, port);
    }

    if (hdesc->socket < 0) {
        ah_errstr = strerror(errno);
        goto error;
    }
    goto cleanup;

  error:
    retval = -1;

  cleanup:
    hcfg_fini(&connect_cfg);
    return retval;
}

/**
 * \brief Close the connection to a Harmony tuning session.
 *
 * Terminates communication with a Harmony session.  Task descriptors
 * associated with this session may continue to be used in a limited
 * capacity.  For example, best point retrieval is available
 * (ah_best()), but new testing values are not (ah_fetch()).
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_close(hdesc_t* hdesc)
{
    if (hdesc->socket < 0) {
        ah_errstr = "Descriptor already closed";
        return -1;
    }

    if (close(hdesc->socket) != 0) {
        snprintf_grow(&ah_errbuf, &ah_errbuflen, "Could not close socket: %s",
                      strerror(errno));
        ah_errstr = ah_errbuf;
    }

    // Inform all search task descriptors of their detached state.
    for (int i = 0; i < hdesc->tlist_len; ++i)
        hdesc->tlist[i]->state = HARMONY_STATE_DETACHED;

    // Reset the descriptor socket to prepare for reuse.
    hdesc->socket = -1;

    return 0;
}

/**
 * \brief Release resources associated with a Harmony client descriptor.
 *
 * If the descriptor is connected to a search session, this function
 * will call ah_leave() on any active search tasks this client is
 * currently participating in.
 *
 * \note This function will not free memory for search definition
 *       structures built using this descriptor.  They must be
 *       released separately using ah_def_fini().
 *
 * \warning This function will free memory for search tasks associated
 *          with this session descriptor.  This will invalidate
 *          pointers that were returned from functions like ah_start()
 *          and ah_join().  Be sure these are no longer in use before
 *          calling this function.
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 */
void ah_free(hdesc_t* hdesc)
{
    if (hdesc) {
        if (hdesc->socket > 0)
            close(hdesc->socket);

        for (int i = 0; i < hdesc->tlist_len; ++i)
            free_task(hdesc->tlist[i]);
        free(hdesc->tlist);

        hmesg_fini(&hdesc->mesg);
        hcfg_fini(&hdesc->cfg);
        free(hdesc->id);
        free(hdesc->buf);
        free(hdesc);
    }
    ah_error_clear();
}

/*
 * Search definition implementation.
 */

/**
 * @}
 *
 * \defgroup api_def Tuning Search Definition Functions
 *
 * These functions are used to define a Harmony tuning search.  They
 * allow users to articulate the details of each tuning variable, such
 * as value domain type (e.g., integers vs. real numbers) and valid
 * value ranges.  To be valid, a search definition must contain at
 * least one tuning variable.
 *
 * @{
 */

/**
 * \brief Generate an empty search definition.
 *
 * This structure is used with functions like ah_def_int() and
 * ah_def_enum() to programmatically define the parameters and
 * configuration environment for a new search.
 *
 * \return Returns a new Harmony search definition descriptor on
 *         success, and `NULL` otherwise.
 */
hdef_t* ah_def_alloc(void)
{
    hdef_t* hdef = calloc(1, sizeof(*hdef));
    if (!hdef) {
        ah_errstr = "Could not allocate memory for new hdef_t";
        return NULL;
    }
    return hdef;
}

/**
 * \brief Generate an empty search definition.
 *
 * This structure is used with functions like ah_def_int() and
 * ah_def_enum() to programmatically define the parameters and
 * configuration environment for a new search.
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 * \param name Name to associate with this search.  If NULL, then the
 *             library will attempt to generate a unique identifier
 *             during ah_start().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_name(hdef_t* hdef, const char* name)
{
    if (hspace_name(&hdef->space, name) != 0) {
        ah_errstr = "Could not allocate memory for session name";
        return -1;
    }
    return 0;
}

/**
 * \brief Generate a tuning session description from a file.
 *
 * Opens and parses a file for configuration and tuning variable
 * declarations.  The session name is copied from the filename.  If
 * `filename` is the string "-", `stdin` will be used instead.  To
 * load a file named "-", include the path (i.e., "./-").
 *
 * The file is parsed line by line.  Each line may contain either a
 * tuning variable or configuration variable definitions.  The line is
 * considered a tuning variable declaration if it begins with a valid
 * tuning variable type.
 *
 * Here are some examples of valid integer-domain variables:
 *
 *     int first_ivar = min:-10 max:10 step:2
 *     int second_ivar = min:1 max:100 # Integers have a default step of 1.
 *
 * Here are some examples of valid real-domain variables:
 *
 *     real first_rvar = min:0 max:1 step:1e-4
 *     real second_rvar = min:-0.5 max:0.5 step:0.001
 *
 * Here are some examples of valid enumerated-domain variables:
 *
 *     enum sort_algo   = bubble, insertion, quick, merge
 *     enum search_algo = forward backward binary # Commas are optional.
 *     enum quotes      = "Cogito ergo sum" \
 *                        "To be, or not to be." \
 *                        "What's up, doc?"
 *
 * Otherwise, the line is considered to be a configuration variable
 * definition, like the following examples:
 *
 *     STRATEGY=pro.so
 *     INIT_POINT = 1.41421, \
 *                  2.71828, \
 *                  3.14159
 *     LAYERS=agg.so:log.so
 *
 *     AGG_TIMES=5
 *     AGG_FUNC=median
 *
 *     LOG_FILE=/tmp/search.log
 *     LOG_MODE=w
 *
 * Variable names must conform to C identifier rules.  Comments are
 * preceded by the hash (#) character, and long lines may be joined
 * with a backslash (\) character.
 *
 * Parsing stops immediately upon error and ah_error() may be
 * used to determine the nature of the error.
 *
 * \param filename Name to associate with this variable.
 *
 * \return Returns a new hdef_t pointer on success, and NULL
 *         otherwise.
 */
hdef_t* ah_def_load(const char* filename)
{
    FILE* fp;
    hdef_t* hdef = ah_def_alloc();
    if (!hdef) return NULL;

    if (ah_def_name(hdef, filename) != 0)
        return NULL;

    if (strcmp(filename, "-") == 0) {
        fp = stdin;
    }
    else {
        fp = fopen(filename, "r");
        if (!fp) {
            snprintf_grow(&ah_errbuf, &ah_errbuflen,
                          "Could not open '%s' for reading: %s",
                          filename, strerror(errno));
            ah_errstr = ah_errbuf;
            return NULL;
        }
    }

    // Begin with an empty buffer string.
    char* buf = NULL;
    int   buflen = 0;
    char* end = NULL;
    int   linenum = 1;
    while (1) {
        char* line;
        int   count = file_read_line(fp, &buf, &buflen,
                                     &line, &end, &ah_errstr);
        if (count <  0) goto error;
        if (count == 0) break; // Loop exit.

        if ((strncmp(line, "int",  3) == 0 && isspace(line[3])) ||
            (strncmp(line, "real", 4) == 0 && isspace(line[4])) ||
            (strncmp(line, "enum", 4) == 0 && isspace(line[4])))
        {
            if (hspace_parse(&hdef->space, line, &ah_errstr) == -1)
                goto error;
        }
        else {
            if (hcfg_parse(&hdef->cfg, line, &ah_errstr) == -1)
                goto error;
        }
        linenum += count;
    }

    goto cleanup;

  error:
    ah_def_free(hdef);
    hdef = NULL;
    snprintf_grow(&ah_errbuf, &ah_errbuflen, "Parse error at %s:%d: %s",
                  filename, linenum, ah_errstr);
    ah_errstr = ah_errbuf;

  cleanup:
    free(buf);
    fclose(fp);
    return hdef;
}

/**
 * \brief Release the resources used by this definition descriptor.
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 */
void ah_def_free(hdef_t* hdef)
{
    if (hdef) {
        hspace_fini(&hdef->space);
        hcfg_fini(&hdef->cfg);
        free(hdef->varloc);
        free(hdef);
    }
}

/**
 * \brief Add an integer-domain variable to the search definition.
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 * \param name Name to associate with this variable.
 * \param min  Minimum range value (inclusive).
 * \param max  Maximum range value (inclusive).
 * \param step Minimum search increment.
 * \param ptr  Address inside the client where `sizeof(long)` bytes may
 *             be written when new values are retrieved via ah_fetch()
 *             or ah_best().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_int(hdef_t* hdef, const char* name,
               long min, long max, long step, long* ptr)
{
    int index = hdef->space.len;

    if (hspace_int(&hdef->space, name, min, max, step, &ah_errstr) != 0)
        return -1;

    if (extend_varloc(&hdef->varloc, &hdef->varloc_cap, &hdef->space) != 0)
        return -1;

    hdef->varloc[index] = ptr;
    return 0;
}

/**
 * \brief Add a real-domain variable to the search definition.
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 * \param name Name to associate with this variable.
 * \param min  Minimum range value (inclusive).
 * \param max  Maximum range value (inclusive).
 * \param step Minimum search increment.
 * \param ptr  Address inside the client where `sizeof(double)` bytes
 *             may be written when new values are retrieved via
 *             ah_fetch() or ah_best().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_real(hdef_t* hdef, const char* name,
                double min, double max, double step, double* ptr)
{
    int index = hdef->space.len;

    if (hspace_real(&hdef->space, name, min, max, step, &ah_errstr) != 0)
        return -1;

    if (extend_varloc(&hdef->varloc, &hdef->varloc_cap, &hdef->space) != 0)
        return -1;

    hdef->varloc[index] = ptr;
    return 0;
}

/**
 * \brief Add an enumerated-domain variable to the search definition.
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 * \param name Name to associate with this variable.
 * \param ptr  Address inside the client where `sizeof(char*)` bytes
 *             may be written when new values are retrieved via
 *             ah_fetch() or ah_best().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_enum(hdef_t* hdef, const char* name, const char** ptr)
{
    int index = hdef->space.len;

    if (hspace_enum(&hdef->space, name, NULL, &ah_errstr) != 0)
        return -1;

    if (extend_varloc(&hdef->varloc, &hdef->varloc_cap, &hdef->space) != 0)
        return -1;

    hdef->varloc[index] = ptr;
    return 0;
}

/**
 * \brief Append a value to an enumerated-domain variable.
 *
 * If the variable does not exist in the search definition, it will be
 * created.
 *
 * \param hdef  Definition descriptor returned from ah_def_alloc() or
 *              ah_def_load().
 * \param name  Name to associate with this variable.
 * \param value String that belongs in this enumeration.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_enum_value(hdef_t* hdef, const char* name, const char* value)
{
    return hspace_enum(&hdef->space, name, value, &ah_errstr);
}

/**
 * \brief Specify the strategy to use for this search.
 *
 * \param hdef     Definition descriptor returned from ah_def_alloc() or
 *                 ah_def_load().
 * \param strategy Filename of the strategy plug-in to use in this session.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_strategy(hdef_t* hdef, const char* strategy)
{
    if (hcfg_set(&hdef->cfg, CFGKEY_STRATEGY, strategy) != 0) {
        ah_errstr = "Could not modify " CFGKEY_STRATEGY " configuration key";
        return -1;
    }
    return 0;
}

/**
 * \brief Specify the list of plug-ins to use for this search.
 *
 * Plug-in layers are specified via a single string of filenames
 * separated by colon (`:`), semicolon (`;`), comma (`,`), or
 * whitespace characters.  The layers are loaded in list order, with
 * each successive layer placed further from the search strategy in
 * the center.
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 * \param list List of plug-ins to load with this session.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_layers(hdef_t* hdef, const char* list)
{
    if (hcfg_set(&hdef->cfg, CFGKEY_LAYERS, list) != 0) {
        ah_errstr = "Could not modify " CFGKEY_LAYERS " configuration key";
        return -1;
    }
    return 0;
}

/**
 * \brief Modify the initial configuration of a new Harmony search.
 *
 * This function allows the user to specify key/value pairs that will
 * exist in the search prior to calling ah_start().
 *
 * \param hdef Definition descriptor returned from ah_def_alloc() or
 *             ah_def_load().
 * \param key  List of plug-ins to load with this session.
 * \param val  List of plug-ins to load with this session.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_def_cfg(hdef_t* hdef, const char* key, const char* val)
{
    if (hcfg_set(&hdef->cfg, key, val) != 0) {
        if (!valid_id(key, strlen(key))) {
            snprintf_grow(&ah_errbuf, &ah_errbuflen,
                          "Invalid configuration key: %s", key);
            ah_errstr = ah_errbuf;
        }
        else {
            ah_errstr = "Could not modify configuration";
        }
        return -1;
    }
    return 0;
}

/*
 * Search task control implementation.
 */

/**
 * @}
 *
 * \defgroup api_task Tuning Search Task Control Functions
 *
 * These functions are used to control Harmony tuning search tasks.
 * They allow client programs to participate in a search by either
 * starting a new searc or joining an existing task.  Once
 * participating, the client may control the search by pausing or
 * resetting it.
 *
 * @{
 */

/**
 * \brief Start a new Harmony tuning search task.
 *
 * After connecting to a session using ah_connect(), this function
 * initiates a new tuning search task as defined by the given search
 * definition.
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 * \param hdef  Definition descriptor returned from ah_def_alloc() or
 *              ah_def_load().
 *
 * \return Returns a new Harmony search task descriptor on success,
 *         and `NULL` otherwise.
 */
htask_t* ah_start(hdesc_t* hdesc, hdef_t* hdef)
{
    // Sanity check input.
    if (hdesc->socket < 0) {
        ah_errstr = "Descriptor not yet connected to a session";
        return NULL;
    }

    if (hdef->space.len < 1) {
        ah_errstr = "No tuning variables defined";
        return NULL;
    }

    for (int i = 0; i < hdef->space.len; ++i) {
        if (hdef->space.dim[i].type == HVAL_STR &&
            hdef->space.dim[i].bounds.e.len == 0)
        {
            ah_errstr = "Enumerated variable (dimension) has no values";
            return NULL;
        }
    }

    htask_t* htask = alloc_task(hdesc);
    if (!htask) {
        ah_errstr = "Could not allocate memory for new htask_t";
        return NULL;
    }

    if (init_task(htask, &hdef->space, &hdef->cfg) != 0)
        goto error;

    // Copy variable location pointers from the search definition.
    if (htask->varloc_cap < hdef->varloc_cap) {
        void** newbuf = realloc(&htask->varloc,
                                hdef->varloc_cap * sizeof(*newbuf));
        if (!newbuf) {
            ah_errstr = "Could not allocate memory for variable locations";
            goto error;
        }
        htask->varloc = newbuf;
        htask->varloc_cap = hdef->varloc_cap;
    }
    memcpy(htask->varloc, hdef->varloc,
           hdef->varloc_cap * sizeof(hdef->varloc));

    // Prepare a Harmony message.
    htask->mesg->state.space = &htask->space;
    htask->mesg->data.cfg = &htask->cfg;

    if (send_request(htask, HMESG_SESSION) != 0)
        goto error;

    // All subsequent messages using this task should be sent to the
    // search who replied to our HMESG_SESSION request.
    //
    htask->dest = hdesc->mesg.src;
    htask->state = HARMONY_STATE_READY;
    return htask;

  error:
    free_task(htask);
    --hdesc->tlist_len;
    return NULL;
}

/**
 * \brief Join an established Harmony tuning search task.
 *
 * After connecting to a session using ah_connect(), this function
 * initiates a new tuning search task as defined by the given search
 * definition.
 *
 * \param hdesc Harmony descriptor returned from ah_alloc().
 * \param name  Name of an existing tuning search in the session.
 *
 * \return Returns a new htask_t pointer on success.  Otherwise, NULL
 *         is returned and ah_error() may be used to determine
 *         the nature of the error.
 */
htask_t* ah_join(hdesc_t* hdesc, const char* name)
{
    // Sanity check input.
    if (hdesc->socket < 0) {
        ah_errstr = "Descriptor not yet connected to a session";
        return NULL;
    }

    if (!name) {
        ah_errstr = "Invalid session name";
        return NULL;
    }

    // Allocate a new task structure and begin initialization.
    htask_t* htask = alloc_task(hdesc);
    if (!htask) {
        ah_errstr = "Could not allocate memory for new htask_t";
        return NULL;
    }

    // Prepare a Harmony message.
    hdesc->mesg.data.string = name;

    // Send the client registration message.
    if (send_request(htask, HMESG_JOIN) != 0)
        goto error;

    if (hdesc->mesg.status != HMESG_STATUS_OK)
        goto error;

    if (init_task(htask, hdesc->mesg.state.space, NULL) != 0)
        goto error;

    // All subsequent messages using this task should be sent to the
    // search who replied to our HMESG_SESSION request.
    //
    htask->dest = hdesc->mesg.src;
    htask->state = HARMONY_STATE_READY;
    return htask;

  error:
    free_task(htask);
    --hdesc->tlist_len;
    return NULL;
}

/**
 * \brief Pause a tuning search task.
 *
 * This instructs the session to halt the production of new trial
 * points for this search task, and return only the best point until
 * the search is resumed.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_pause(htask_t* htask)
{
    if (ah_set_cfg(htask, CFGKEY_PAUSED, "0") == NULL)
        return -1;
    else
        return 0;
}

/**
 * \brief Resume a tuning search task.
 *
 * This instructs the session to continue the production of new trial
 * points for this search task.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_resume(htask_t* htask)
{
    if (ah_set_cfg(htask, CFGKEY_PAUSED, "0") == NULL)
        return -1;
    else
        return 0;
}

/**
 * \brief Resume a tuning search task.
 *
 * This instructs the session to continue the production of new trial
 * points for this search task.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_restart(htask_t* htask)
{
    if (htask->state < HARMONY_STATE_ATTACHED) {
        ah_errstr = "Cannot restart a detached task descriptor";
        return -1;
    }

    // Prepare a Harmony message.
    htask->mesg->data.string = "restart";
    return send_request(htask, HMESG_COMMAND);
}

/**
 * \brief Leave a tuning search task.
 *
 * End participation with the search task.  The task descriptor will
 * remain available for limited set of functions (e.g., ah_best(),
 * ah_get_int(), etc.), but functions that participation in a task
 * (e.g., ah_fetch()) will result in error.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_leave(htask_t* htask)
{
    if (htask->state < HARMONY_STATE_ATTACHED) {
        ah_errstr = "Already detached from search task";
        return -1;
    }

    // Prepare a Harmony message.
    htask->mesg->data.string = "leave";
    return send_request(htask, HMESG_COMMAND);
}

/**
 * \brief Kill a tuning search task.
 *
 * If the task descriptor is currently attached, this function will
 * send a kill request to the remote session.  This instructs the
 * session to release the resources it holds for the task.  Any future
 * messages received by the session destine for this task will result
 * in error.
 *
 * \note Errors encountered while sending the kill request will be
 *       ignored, in case the search has already been killed by a
 *       competing client.
 *
 * If the task descriptor is currently detached (via ah_leave()), the
 * kill request is skipped.
 *
 * Either way, this function will release client-local resources
 * reserved for this task descriptor.  Using the now defunct task
 * descriptor is an error and will result in undefined behavior.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_kill(htask_t* htask)
{
    if (htask->state >= HARMONY_STATE_ATTACHED) {
        // Prepare a Harmony message.
        htask->mesg->data.string = "kill";

        if (send_request(htask, HMESG_COMMAND) != 0)
            return -1;
    }

    // Remove this task from the session descriptor task list.
    int newlen = --htask->hdesc->tlist_len;
    for (int i = 0; i < newlen; ++i) {
        if (htask->hdesc->tlist[i] == htask) {
            htask->hdesc->tlist[i] = htask->hdesc->tlist[newlen];
            break;
        }
    }

    // Free client-local task resources.
    free_task(htask);
    return 0;
}

/*
 * Search task interaction implementation.
 */

/**
 * @}
 *
 * \defgroup api_search Search Task Interaction Functions
 *
 * These functions are used by the client after it is participating in
 * a running search.  This is how clients perform activities such as
 * retrieving a point to test and reporting its associated
 * performance.
 *
 * @{
 */

/**
 * \brief Bind a client address to an integer-domain search variable
 *        (dimension).
 *
 * This function associates an address inside the client with a search
 * variable defined using ah_def_int().  Upon ah_fetch(), the value
 * retrieved for this tuning variable will be stored in `sizeof(long)`
 * bytes starting at address `ptr`.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param name  Search variable defined using ah_def_int().
 * \param ptr   Address inside the client where `sizeof(long)` bytes may
 *              be written to store the current testing value.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_bind_int(htask_t* htask, const char* name, long* ptr)
{
    return set_varloc(htask, name, ptr);
}

/**
 * \brief Bind a client address to an real-domain search variable
 *        (dimension).
 *
 * This function associates an address inside the client with a search
 * variable defined using ah_def_real().  Upon ah_fetch(), the value
 * retrieved for this tuning variable will be stored in
 * `sizeof(double)` bytes starting at address `ptr`.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param name  Search variable defined using ah_def_real().
 * \param ptr   Address inside the client where `sizeof(double)` bytes
 *              may be written to store the current testing value.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_bind_real(htask_t* htask, const char* name, double* ptr)
{
    return set_varloc(htask, name, ptr);
}

/**
 * \brief Bind a client address to an enumerated string-based search
 *        variable (dimension).
 *
 * This function associates an address inside the client with a search
 * variable defined using ah_def_enum().  Upon ah_fetch(), the value
 * retrieved for this tuning variable will be stored in
 * `sizeof(char*)` bytes starting at address `ptr`.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param name  Search variable defined using ah_def_int().
 * \param ptr   Address inside the client where `sizeof(char*)` bytes
 *              may be written to store the current testing value.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_bind_enum(htask_t* htask, const char* name, const char** ptr)
{
    return set_varloc(htask, name, ptr);
}

/**
 * \brief Return the current value of an integer-domain search
 *        variable.
 *
 * Finds an integer-domain tuning variable given its name and returns
 * its current value.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param name  Name of a tuning variable declared with ah_def_int().
 *
 * \return Returns the current value of a tuning variable, if the
 *         given name matches an integer-domain search variable.
 *         Otherwise, LONG_MIN is returned and ah_error() may be used
 *         to determine the nature of the error.
 */
long ah_get_int(htask_t* htask, const char* name)
{
    int idx = find_var(htask, name);

    if (idx < 0) {
        ah_errstr = "Could not find variable by name";
        return LONG_MIN;
    }

    if (htask->space.dim[idx].type != HVAL_INT) {
        ah_errstr = "Variable type mismatch";
        return LONG_MIN;
    }

    if (htask->varloc[idx])
        return *(long*)htask->varloc[idx];
    else
        return htask->curr->term[idx].value.i;
}

/**
 * \brief Return the current value of a real-domain search variable.
 *
 * Finds a real-domain tuning variable given its name and returns its
 * current value.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param name  Name of a tuning variable declared with ah_def_real().
 *
 * \return Returns the current value of a tuning variable, if the
 *         given name matches a real-domain search variable.
 *         Otherwise, NAN is returned and and ah_error() may be used
 *         to determine the nature of the error.
 */
double ah_get_real(htask_t* htask, const char* name)
{
    int idx = find_var(htask, name);

    if (idx < 0) {
        ah_errstr = "Could not find variable by name";
        return NAN;
    }

    if (htask->space.dim[idx].type != HVAL_REAL) {
        ah_errstr = "Variable type mismatch";
        return NAN;
    }

    if (htask->varloc[idx])
        return *(double*)htask->varloc[idx];
    else
        return htask->curr->term[idx].value.r;
}

/**
 * \brief Return the current value of an enumerated-domain search
 *        variable.
 *
 * Finds an enumerated-domain tuning variable given its name and
 * returns its current value.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param name  Name of a tuning variable declared with ah_def_enum().
 *
 * \return Returns the current value of a tuning variable, if the
 *         given name matches an enumerated-domain search variable.
 *         Otherwise, NULL is returned and ah_error() may be used to
 *         determine the nature of the error.
 */
const char* ah_get_enum(htask_t* htask, const char* name)
{
    int idx = find_var(htask, name);

    if (idx < 0) {
        ah_errstr = "Could not find variable by name";
        return NULL;
    }

    if (htask->space.dim[idx].type != HVAL_STR) {
        ah_errstr = "Variable type mismatch";
        return NULL;
    }

    if (htask->varloc[idx])
        return *(const char**)htask->varloc[idx];
    else
        return htask->curr->term[idx].value.s;
}

/**
 * \brief Get a key value from the search's configuration.
 *
 * Retrieve the string value associated with the given key in the
 * search task's configuration system.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param key   Config key to search for on the session.
 *
 * \return Returns a c-style string on success, and `NULL` otherwise.
 */
const char* ah_get_cfg(htask_t* htask, const char* key)
{
    hmesg_t* mesg = htask->mesg;

    if (!valid_id(key, strlen(key))) {
        ah_errstr = "Invalid key string";
        return NULL;
    }

    if (htask->hdesc->socket >= 0) {
        // Prepare a Harmony message.
        mesg->data.string = key;

        if (send_request(htask, HMESG_GETCFG) != 0)
            return NULL;

        if (mesg->status != HMESG_STATUS_OK)
            return NULL;

        char* val = strchr(mesg->data.string, '=');
        if (!val) {
            ah_errstr = "Malformed message received from server";
            return NULL;
        }
        *(val++) = '\0';

        // Store a copy in the task descriptors configuration cache.
        if (hcfg_set(&htask->cfg, mesg->data.string, val) != 0) {
            ah_errstr = "Could not modify task-local configuration cache";
            return NULL;
        }
    }
    return hcfg_get(&htask->cfg, mesg->data.string);
}

/**
 * \brief Set a new key/value pair in the search's configuration.
 *
 * Writes the new key/value pair into the search's run-time
 * configuration database.  If the key exists in the database, its
 * value is overwritten.  If val is `NULL`, the key will be erased
 * from the database.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param key   Config key to modify in the search.
 * \param val   Config value to associate with the key.
 *
 * \return Returns the original key value on success.  If the key did not
 *         exist prior to this call, an empty string ("") is returned.
 *         Otherwise, `NULL` is returned on error.
 *
 * \note The buffer which stores the previous configuration value may
 *       be used by other Harmony Client API functions.  If you need
 *       the value to persist beyond the next API call, you must make
 *       a copy.
 */
const char* ah_set_cfg(htask_t* htask, const char* key, const char* val)
{
    hmesg_t* mesg = htask->mesg;

    if (htask->state < HARMONY_STATE_ATTACHED) {
        ah_errstr = "Cannot modify the configuration of a detached task";
        return NULL;
    }

    if (!valid_id(key, strlen(key))) {
        ah_errstr = "Invalid key string";
        return NULL;
    }

    if (!val)
        val = "";

    // Prepare a Harmony message.
    snprintf_grow(&htask->hdesc->buf, &htask->hdesc->buflen,
                  "%s=%s", key, val);
    mesg->data.string = htask->hdesc->buf;

    if (send_request(htask, HMESG_SETCFG) != 0)
        return NULL;

    if (mesg->status != HMESG_STATUS_OK)
        return NULL;

    return mesg->data.string;
}

/**
 * \brief Fetch a new configuration from the Harmony search.
 *
 * If a new configuration is available, this function will update the
 * values of all registered variables.  Otherwise, it will configure
 * the system to run with the best known configuration thus far.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 0 if no registered variables were modified, 1 if
 *         any registered variables were modified, and -1 otherwise.
 */
int ah_fetch(htask_t* htask)
{
    if (htask->state == HARMONY_STATE_READY) {
        // Prepare a Harmony message.
        if (send_request(htask, HMESG_FETCH) != 0)
            return -1;

        if (htask->mesg->status == HMESG_STATUS_BUSY) {
            if (!htask->best.id) {
                // No best point is available yet.  Do not set values.
                return 0;
            }

            // Set current point to best point.
            htask->curr = &htask->best;
        }
        else if (htask->mesg->status == HMESG_STATUS_OK) {
            if (hpoint_copy(&htask->test, htask->mesg->data.point) != 0) {
                ah_errstr = "Error copying test point data";
                return -1;
            }

            if (hpoint_align(&htask->test, &htask->space) != 0) {
                ah_errstr = "Error aligning test point data";
                return -1;
            }
            htask->curr = &htask->test;
            htask->state = HARMONY_STATE_TESTING;
        }
        else {
            ah_errstr = "Invalid message received from server";
            return -1;
        }
    }

    // Update the variables from the content of the message.
    if (write_values(htask) != 0)
        return -1;

    // Initialize our internal performance array.
    for (int i = 0; i < htask->perf.len; ++i)
        htask->perf.obj[i] = NAN;

    // Client variables were changed.  Inform the user by returning 1.
    return 1;
}

/**
 * \brief Report the performance of a configuration to the Harmony
 *        search.
 *
 * Sends a complete performance report regarding the current
 * configuration to the Harmony session.  The `perf` parameter should
 * point to a floating-point double.
 *
 * For multi-objective search problems, the `perf` parameter should
 * point to an array of floating-point doubles containing all
 * performance values.  If all performance values have already been
 * reported via ah_report_one(), then `NULL` may be passed as the
 * performance pointer.  Unreported performance values will result in
 * error.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param perf  Performance vector for the current configuration.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_report(htask_t* htask, double* perf)
{
    if (htask->state < HARMONY_STATE_TESTING)
        return 0;

    if (perf) {
        memcpy(htask->perf.obj, perf,
               sizeof(*htask->perf.obj) * htask->perf.len);
    }

    for (int i = 0; i < htask->perf.len; ++i) {
        if (isnan(htask->perf.obj[i])) {
            ah_errstr = "Insufficient performance values to report";
            return -1;
        }
    }

    // Prepare a Harmony message.
    htask->mesg->data.point = &htask->test;
    htask->mesg->data.perf = &htask->perf;

    if (send_request(htask, HMESG_REPORT) != 0)
        return -1;

    if (htask->mesg->status != HMESG_STATUS_OK)
        return -1;

    htask->state = HARMONY_STATE_READY;
    return 0;
}

/**
 * \brief Report a single performance value for the current
 *        configuration.
 *
 * Allows performance values to be reported one at a time for
 * multi-objective search problems.
 *
 * Note that this function only caches values to send to the Harmony
 * search.  Once all values have been reported, ah_report() must
 * be called (passing `NULL` as the performance argument).
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 * \param index Objective index of the value to report.
 * \param value Performance measured for the current configuration.
 *
 * \return Returns 0 on success, and -1 otherwise.
 */
int ah_report_one(htask_t* htask, int index, double value)
{
    if (index != 0 || index >= htask->perf.len) {
        ah_errstr = "Invalid performance index";
        return -1;
    }

    if (htask->state == HARMONY_STATE_TESTING)
        htask->perf.obj[index] = value;

    return 0;
}

/**
 * \brief Set values under Harmony's control to the best known
 *        configuration.
 *
 * Bound memory (from using functions like ah_bind_int()) or values
 * retrieved (from functions like ah_get_real()) will be the best
 * known input values from the search thus far.
 *
 * If no best configuration exists (i.e., before any configurations
 * have been evaluated), this function will return an error.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 1 on success, and -1 otherwise.
 */
int ah_best(htask_t* htask)
{
    int retval = 0;

    // Make sure our best known point is valid.
    if (!htask->best.id) {
        ah_errstr = "Best point currently unavailable";
        return -1;
    }

    htask->curr = &htask->best;
    if (write_values(htask) != 0)
        return -1;

    return retval;
}

/**
 * \brief Retrieve the convergence state of the current search.
 *
 * \param htask Task descriptor returned from ah_start() or ah_join().
 *
 * \return Returns 1 if the search has converged, 0 if it has not,
 *         and -1 on error.
 */
int ah_converged(htask_t* htask)
{
    const char* retval = ah_get_cfg(htask, CFGKEY_CONVERGED);
    return (retval && retval[0] == '1');
}

/*
 * Harmony error reporting implementation.
 */

/**
 * @}
 *
 * \defgroup api_error Harmony Error Reporting Functions
 *
 * These functions may be used to retrieve information about the most
 * recent failure within Harmony API functions.
 *
 * @{
 */

/**
 * \brief Retrieve the most recent Harmony API error string.
 *
 * \return Returns a pointer to a string that describes the latest
 *         Harmony error, or `NULL` if no error has occurred since
 *         the last call to ah_error_clear().
 *
 * \note The buffer which stores the most recent error string may be
 *       overwritten by other Harmony Client API functions.  If you
 *       need the value to persist beyond the next API call, you must
 *       make a copy.
 */
const char* ah_error(void)
{
    return ah_errstr;
}

/**
 * \brief Clear the current task-level error string.
 *
 */
void ah_error_clear(void)
{
    free(ah_errbuf);
    ah_errstr = NULL;
    ah_errbuf = NULL;
    ah_errbuflen = 0;
}

/**
 * @}
 */

/*
 * Deprecated API function implementation.
 */
typedef struct hdesc_compat {
    hdesc_t* hdesc;
    hdef_t*  hdef;
    htask_t* htask;
} hdesc_compat_t;

static hdesc_compat_t* harmony_compat_find(hdesc_t* hdesc);

static hdesc_compat_t* harmony_compat_list;
static int harmony_compat_list_len;
static int harmony_compat_list_cap;

hdesc_t* harmony_init(int* argc, char*** argv)
{
    if (harmony_compat_list_len == harmony_compat_list_cap) {
        if (array_grow(&harmony_compat_list, &harmony_compat_list_cap,
                       sizeof(*harmony_compat_list)) != 0)
        {
            return NULL;
        }
    }

    hdesc_compat_t* compat = &harmony_compat_list[ harmony_compat_list_len ];
    compat->hdesc = ah_alloc();
    compat->hdef  = ah_def_alloc();
    compat->htask = NULL;
    ++harmony_compat_list_len;

    if (argc && argv)
        ah_args(compat->hdesc, argc, *argv);

    return compat->hdesc;
}

int harmony_parse_args(hdesc_t* hdesc, int argc, char** argv)
{
    return ah_args(hdesc, &argc, argv);
}

void harmony_fini(hdesc_t* hdesc)
{
    int index;
    for (index = 0; index < harmony_compat_list_len; ++index) {
        if (harmony_compat_list[index].hdesc == hdesc)
            break;
    }
    if (index == harmony_compat_list_len)
        return;

    ah_def_free(harmony_compat_list[index].hdef);
    ah_free(hdesc);

    int newlen = --harmony_compat_list_len;
    if (newlen && newlen != index)
        harmony_compat_list[index] = harmony_compat_list[newlen];
}

int harmony_int(hdesc_t* hdesc, const char* name,
                long min, long max, long step)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_def_int(compat->hdef, name, min, max, step, NULL);
}

int harmony_real(hdesc_t* hdesc, const char* name,
                 double min, double max, double step)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_def_real(compat->hdef, name, min, max, step, NULL);
}

int harmony_enum(hdesc_t* hdesc, const char* name, const char* value)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_def_enum_value(compat->hdef, name, value);
}

int harmony_session_name(hdesc_t* hdesc, const char* name)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return hspace_name(&compat->hdef->space, name);
}

int harmony_strategy(hdesc_t* hdesc, const char* strategy)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_def_strategy(compat->hdef, strategy);
}

int harmony_layers(hdesc_t* hdesc, const char* list)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_def_layers(compat->hdef, list);
}

int harmony_launch(hdesc_t* hdesc, const char* host, int port)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    if (ah_connect(hdesc, host, port) != 0)
        return -1;

    compat->htask = ah_start(hdesc, compat->hdef);
    if (!compat->htask)
        return -1;

    return 0;
}

int harmony_id(hdesc_t* hdesc, const char* id)
{
    return ah_id(hdesc, id);
}

int harmony_bind_int(hdesc_t* hdesc, const char* name, long* ptr)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);

    if (!compat->htask) {
        for (int i = 0; i < compat->hdef->space.len; ++i) {
            if (strcmp(compat->hdef->space.dim[i].name, name) == 0) {
                compat->hdef->varloc[i] = ptr;
                return 0;
            }
        }
        return -1;
    }
    return ah_bind_int(compat->htask, name, ptr);
}

int harmony_bind_real(hdesc_t* hdesc, const char* name, double* ptr)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);

    if (!compat->htask) {
        for (int i = 0; i < compat->hdef->space.len; ++i) {
            if (strcmp(compat->hdef->space.dim[i].name, name) == 0) {
                compat->hdef->varloc[i] = ptr;
                return 0;
            }
        }
        return -1;
    }
    return ah_bind_real(compat->htask, name, ptr);
}

int harmony_bind_enum(hdesc_t* hdesc, const char* name, const char** ptr)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);

    if (!compat->htask) {
        for (int i = 0; i < compat->hdef->space.len; ++i) {
            if (strcmp(compat->hdef->space.dim[i].name, name) == 0) {
                compat->hdef->varloc[i] = ptr;
                return 0;
            }
        }
        return -1;
    }
    return ah_bind_enum(compat->htask, name, ptr);
}

int harmony_join(hdesc_t* hdesc, const char* host, int port, const char* name)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);

    if (!compat->htask) {
        if (ah_connect(hdesc, host, port) != 0)
            return -1;

        compat->htask = ah_join(hdesc, name);
        if (!compat->htask)
            return -1;
    }

    return 0;
}

int harmony_leave(hdesc_t* hdesc)
{
    return ah_close(hdesc);
}

long harmony_get_int(hdesc_t* hdesc, const char* name)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_get_int(compat->htask, name);
}

double harmony_get_real(hdesc_t* hdesc, const char* name)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_get_real(compat->htask, name);
}

const char* harmony_get_enum(hdesc_t* hdesc, const char* name)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_get_enum(compat->htask, name);
}

char* harmony_getcfg(hdesc_t* hdesc, const char* key)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    const char* retval;

    if (!compat->htask)
        retval = hcfg_get(&compat->hdef->cfg, key);
    else
        retval = ah_get_cfg(compat->htask, key);

    return stralloc( retval );
}

char* harmony_setcfg(hdesc_t* hdesc, const char* key, const char* val)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    const char* retval;

    if (!compat->htask) {
        retval = hcfg_get(&compat->hdef->cfg, key);
        hcfg_set(&compat->hdef->cfg, key, val);
    }
    else {
        retval = ah_set_cfg(compat->htask, key, val);
    }
    return stralloc( retval );
}

int harmony_fetch(hdesc_t* hdesc)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_fetch(compat->htask);
}

int harmony_report(hdesc_t* hdesc, double perf)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_report(compat->htask, &perf);
}

int harmony_report_one(hdesc_t* hdesc, int index, double value)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_report_one(compat->htask, index, value);
}

int harmony_best(hdesc_t* hdesc)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_best(compat->htask);
}

int harmony_converged(hdesc_t* hdesc)
{
    hdesc_compat_t* compat = harmony_compat_find(hdesc);
    return ah_converged(compat->htask);
}

const char* harmony_error_string(hdesc_t* hdesc)
{
    return ah_error();
}

void harmony_error_clear(hdesc_t* hdesc)
{
    ah_error_clear();
}

hdesc_compat_t* harmony_compat_find(hdesc_t* hdesc)
{
    for (int i = 0; i < harmony_compat_list_len; ++i)
        if (harmony_compat_list[i].hdesc == hdesc)
            return &harmony_compat_list[i];
    return NULL;
}

/*
 * Internal helper function implementation.
 */
htask_t* alloc_task(hdesc_t* hdesc)
{
    // Allocate a new task structure and begin initialization.
    htask_t* htask = calloc(1, sizeof(*htask));
    if (!htask) {
        ah_errstr = "Could not allocate memory for new htask_t";
        return NULL;
    }

    // Extend the Harmony descriptor task list, if necessary.
    if (hdesc->tlist_len == hdesc->tlist_cap) {
        if (array_grow(&hdesc->tlist, &hdesc->tlist_cap,
                       sizeof(*hdesc->tlist)) != 0)
        {
            ah_errstr = "Could not grow htask_t list";
            return NULL;
        }
    }

    htask->hdesc = hdesc;
    htask->mesg = &hdesc->mesg;

    // Add the new task to the back of the task list.
    hdesc->tlist[ hdesc->tlist_len ] = htask;
    ++hdesc->tlist_len;

    return htask;
}

int init_task(htask_t* htask, const hspace_t* space, const hcfg_t* cfg)
{
    // Copy the input space from the search definition.
    if (hspace_copy(&htask->space, space) != 0) {
        ah_errstr = "Could not copy hspace_t from definition";
        return -1;
    }

    // Provide a default name, if necessary.
    if (!htask->space.name) {
        if (hspace_name(&htask->space, generate_id(htask->hdesc, 1)) != 0) {
            ah_errstr = "Could not provide a default search name";
            return -1;
        }
    }

    // Begin with global default configuration directives.
    if (hcfg_init(&htask->cfg) != 0) {
        ah_errstr = "Could not initialize task configuration";
        return -1;
    }

    // Merge search definition configuration directives, if necessary.
    if (cfg) {
        if (hcfg_merge(&htask->cfg, cfg) != 0) {
            ah_errstr = "Could not merge search configuration directives";
            return -1;
        }
    }

    // Merge overriding directives taken from the environment or command line.
    if (hcfg_merge(&htask->cfg, &htask->hdesc->cfg) != 0) {
        ah_errstr = "Could not merge overriding configuration directives";
        return -1;
    }

    if (extend_varloc(&htask->varloc, &htask->varloc_cap, &htask->space) != 0)
        return -1;

    if (extend_perf(htask) != 0)
        return -1;

    return 0;
}

void free_task(htask_t* htask)
{
    if (htask) {
        hperf_fini(&htask->perf);
        hpoint_fini(&htask->best);
        hpoint_fini(&htask->test);
        free(htask->varloc);
        hcfg_fini(&htask->cfg);
        hspace_fini(&htask->space);
        free(htask);
    }
}

int find_var(htask_t* htask, const char* name)
{
    for (int idx = 0; idx < htask->space.len; ++idx) {
        if (strcmp(htask->space.dim[idx].name, name) == 0)
            return idx;
    }
    return -1;
}

int extend_perf(htask_t* htask)
{
    int perf_len;

    perf_len = hcfg_int(&htask->cfg, CFGKEY_PERF_COUNT);
    if (perf_len < 1) {
        ah_errstr = "Invalid value for " CFGKEY_PERF_COUNT;
        return -1;
    }

    if (htask->perf.cap < perf_len) {
        if (hperf_init(&htask->perf, perf_len) != 0) {
            ah_errstr = "Could not allocate performance array";
            return -1;
        }
    }

    htask->perf.len = perf_len;
    return 0;
}

int extend_varloc(void*** varloc, int* varloc_cap, hspace_t* space)
{
    while (*varloc_cap < space->len) {
        if (array_grow(varloc, varloc_cap, sizeof(**varloc)) != 0) {
            ah_errstr = "Could not extend bound variable address array";
            return -1;
        }
    }
    return 0;
}

char* generate_id(hdesc_t* hdesc, int suffix)
{
    static int id_count = 0;
    int maxlen = sysconf(_SC_HOST_NAME_MAX) + 32; // Add room for PID and

    while (hdesc->buflen < maxlen) {
        if (array_grow(&hdesc->buf, &hdesc->buflen,
                       sizeof(*hdesc->buf)) != 0)
        {
            ah_errstr = "Could not extend temporary string buffer";
            return NULL;
        }
    }

    gethostname(hdesc->buf, maxlen);
    hdesc->buf[maxlen] = '\0'; // Portability safety.

    if (suffix)
        sprintf(hdesc->buf + strlen(hdesc->buf),
                "_%d_%x", getpid(), ++id_count);
    else
        sprintf(hdesc->buf + strlen(hdesc->buf),
                "_%d", getpid());

    return hdesc->buf;
}

int send_request(htask_t* htask, hmesg_type msg_type)
{
    hmesg_t* mesg = htask->mesg;

    mesg->src = 0;
    mesg->dest = htask->dest;
    mesg->type = msg_type;
    mesg->status = HMESG_STATUS_REQ;

    mesg->state.space = &htask->space;
    mesg->state.best = &htask->best;
    mesg->state.client = htask->hdesc->id;

    if (mesg_send(htask->hdesc->socket, mesg) < 1) {
        ah_errstr = "Error sending Harmony message to server";
        return -1;
    }

    if (mesg_recv(htask->hdesc->socket, mesg) < 1) {
        ah_errstr = "Error retrieving Harmony message from server";
        return -1;
    }

    // Special message indicating the search task was killed unexpectedly.
    if (mesg->type == HMESG_UNKNOWN && mesg->status == HMESG_STATUS_FAIL) {
        ah_errstr = "Search no longer exists";
        return -1;
    }

    if (mesg->type != msg_type) {
        ah_errstr = "Server response message mismatch";
        return -1;
    }

    if (mesg->status == HMESG_STATUS_FAIL) {
        ah_errstr = htask->mesg->data.string;
        return -1;
    }

    // Update local state from server.
    if (htask->space.id < mesg->state.space->id) {
        if (hspace_copy(&htask->space, mesg->state.space) != 0) {
            ah_errstr = "Could not update session search space";
            return -1;
        }
    }
    if (htask->best.id < mesg->state.best->id) {
        if (hpoint_copy(&htask->best, htask->mesg->state.best) != 0) {
            ah_errstr = "Could not update best known point";
            return -1;
        }
        if (hpoint_align(&htask->best, &htask->space) != 0) {
            ah_errstr = "Could not align best point to search space";
            return -1;
        }
    }
    return 0;
}

int set_varloc(htask_t* htask, const char* name, void* ptr)
{
    int idx = find_var(htask, name);
    if (idx < 0) {
        ah_errstr = "Tuning variable not found";
        return -1;
    }

    if (extend_varloc(&htask->varloc, &htask->varloc_cap, &htask->space) != 0)
        return -1;

    htask->varloc[idx] = ptr;
    return 0;
}

int write_values(htask_t* htask)
{
    if (htask->space.len != htask->curr->len) {
        ah_errstr = "Invalid internal point structure";
        return -1;
    }

    for (int i = 0; i < htask->space.len; ++i) {
        const hval_t* val = &htask->curr->term[i];

        if (!htask->varloc[i])
            continue;

        switch (htask->space.dim[i].type) {
        case HVAL_INT:  *(long*)       htask->varloc[i] = val->value.i; break;
        case HVAL_REAL: *(double*)     htask->varloc[i] = val->value.r; break;
        case HVAL_STR:  *(const char**)htask->varloc[i] = val->value.s; break;
        default:
            ah_errstr = "Invalid space dimension value type";
            return -1;
        }
    }
    return 0;
}
