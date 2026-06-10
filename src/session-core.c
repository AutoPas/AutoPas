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
#define _XOPEN_SOURCE 600 // Needed for srand48() and S_ISSOCK.

#include "session-core.h"
#include "hspace.h"
#include "hcfg.h"
#include "hmesg.h"
#include "hplugin.h"
#include "hutil.h"
#include "hsockutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/select.h>

/*
 * Structure to encapsulate a Harmony plug-in and its run-time state.
 */
typedef struct pstack {
    hplugin_t plugin;

    // Index list of trials waiting in the generation side of this layer.
    int* wait_generate;
    int wait_generate_len;
    int wait_generate_cap;

    // Index list of trials waiting in the analyze side of this layer.
    int* wait_analyze;
    int wait_analyze_len;
    int wait_analyze_cap;
} pstack_t;

/*
 * Structure to encapsulate the state of a single search instance.
 */
typedef struct hsearch {
    int open;

    hspace_t space;       // Search space definition.
    hcfg_t   cfg;         // Configuration environment.
    hpoint_t best;        // Best known point.

    pstack_t*  pstack;     // Plug-in stack.
    int        pstack_len;
    int        pstack_cap;

    // Search specific pseudo-random number generator state.
    unsigned short xsubi[3];

    // Workflow state.
    hflow_t flow;
    int     curr_layer;
    int     paused_id;

    // List of all points generated, but not yet returned to the strategy.
    htrial_t* pending;
    int       pending_cap;
    int       pending_len;

    // These variables control the capacity of the pending list.
    int clients;
    int per_client;

    // List of all trials (point/performance pairs) waiting for client fetch.
    int* ready;
    int  ready_head;
    int  ready_tail;
    int  ready_cap;

    char* buf;
    int   buf_len;
    const char* errmsg;
} hsearch_t;

/*
 * Callback registration system.
 */
typedef struct callback {
    int        fd;        // Listen on this file descriptor for incoming data.
    int        layer_idx; // Return to this layer index when data is ready.
    cb_func_t  func;      // Call this function to process incoming data.
    void*      data;      // Instance-specific data pointer.
    hsearch_t* search;    // Search this callback is associated with.
} callback_t;

callback_t* cbs; // List of callbacks.
int         cbs_len;
int         cbs_cap;

/*
 * Variables used for select().
 */
struct timeval  polltime;
struct timeval* pollstate;
fd_set fds;
int maxfd;

/*
 * Other global variables.
 */
static const char* home_dir;
static hsearch_t** slist;
static int slist_cap;

/*
 * Pointer to the current search instance for use from within the
 * functions exported for pluggable modules.
 */
static hsearch_t* current_search;
const  hcfg_t*    search_cfg;

/*
 * Base search structure management prototypes.
 */
static int  open_search(hsearch_t* search, hmesg_t* mesg);
static void close_search(hsearch_t* search);
static void free_search(hsearch_t* search);
static int  load_strategy(hsearch_t* search, const char* file);
static int  load_layers(hsearch_t* search, const char* list);
static int  load_plugin(hsearch_t* search, const char* filename,
                        hplugin_type_t type);

/*
 * Internal helper function prototypes.
 */
static int        generate_trial(hsearch_t* search);
static int        plugin_workflow(hsearch_t* search, int trial_idx);
static int        workflow_transition(hsearch_t* search, int trial_idx);
static hsearch_t* find_search(hmesg_t* mesg);
static int        handle_callback(callback_t* cb);
static int        handle_session(hsearch_t* search, hmesg_t* mesg);
static int        handle_join(hsearch_t* search, hmesg_t* mesg);
static int        handle_getcfg(hsearch_t* search, hmesg_t* mesg);
static int        handle_setcfg(hsearch_t* search, hmesg_t* mesg);
static int        handle_best(hmesg_t* mesg);
static int        handle_fetch(hsearch_t* search, hmesg_t* mesg);
static int        handle_report(hsearch_t* search, hmesg_t* mesg);
static int        handle_reject(hsearch_t* search, int trial_idx);
static int        handle_command(hsearch_t* search, hmesg_t* mesg);
static int        handle_wait(hsearch_t* search, int trial_idx);
static int        extend_lists(hsearch_t* search, int target_cap);
static void       reverse_array(void* ptr, int head, int tail);
static int        update_state(hmesg_t* mesg, hsearch_t* search);
static void       set_current(hsearch_t* search);

/*
 * Core session routines begin here.
 */
int main(int argc, char* argv[])
{
    struct stat sb;
    int retval;
    fd_set ready_fds;
    hmesg_t mesg = HMESG_INITIALIZER;

    if (argc < 2) {
        fprintf(stderr, "%s should not be launched manually.\n", argv[0]);
        return -1;
    }
    home_dir = argv[1];

    // Check that we have been launched correctly by checking that
    // STDIN_FILENO is a socket descriptor.
    //
    // Print an error message to stderr if this is not the case.  This
    // should be the only time an error message is printed to stdout
    // or stderr.
    //
    if (fstat(STDIN_FILENO, &sb) < 0) {
        perror("Could not determine the status of STDIN");
        return -1;
    }

    if (!S_ISSOCK(sb.st_mode)) {
        fprintf(stderr, "%s should not be launched manually.\n", argv[0]);
        return -1;
    }

    // Ignore SIGINT.  Graceful exit comes from the parent process
    // closing their end of our STDIN file descriptor.
    //
    if (signal(SIGINT, SIG_IGN) == SIG_ERR) {
        perror("Error in requesting to ignore SIGINT");
        return -1;
    }

    // Initialize global data structures.
    pollstate = &polltime;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    maxfd = STDIN_FILENO;

    while (1) {
        ready_fds = fds;
        retval = select(maxfd + 1, &ready_fds, NULL, NULL, pollstate);
        if (retval < 0) {
            perror("Error during main select loop of session-core");
            break;
        }

        // Launch callbacks, if needed.
        for (int i = 0; i < cbs_len; ++i) {
            if (FD_ISSET(cbs[i].fd, &ready_fds))
                handle_callback(&cbs[i]);
        }

        // Handle hmesg_t, if needed.
        if (FD_ISSET(STDIN_FILENO, &ready_fds)) {
            retval = mesg_recv(STDIN_FILENO, &mesg);
            if (retval == 0) break;
            if (retval <  0) {
                perror("Error receiving message in session-core");
                break;
            }

            hsearch_t* search = find_search(&mesg);
            if (!search)
                goto error;

            set_current(search);
            search->flow.status = HFLOW_ACCEPT;

            if (search->open) {
                hcfg_set(&search->cfg, CFGKEY_CURRENT_CLIENT,
                         mesg.state.client);
            }

            switch (mesg.type) {
            case HMESG_SESSION: retval = handle_session(search, &mesg); break;
            case HMESG_JOIN:    retval = handle_join(search, &mesg); break;
            case HMESG_GETCFG:  retval = handle_getcfg(search, &mesg); break;
            case HMESG_SETCFG:  retval = handle_setcfg(search, &mesg); break;
            case HMESG_BEST:    retval = handle_best(&mesg); break;
            case HMESG_FETCH:   retval = handle_fetch(search, &mesg); break;
            case HMESG_REPORT:  retval = handle_report(search, &mesg); break;
            case HMESG_COMMAND: retval = handle_command(search, &mesg); break;

            default:
                search->errmsg = "Invalid message type";
                goto error;
            }
            if (retval != 0)
                goto error;

            if (update_state(&mesg, search) != 0) {
                search->errmsg = "Could not update message session state";
                goto error;
            }

            hcfg_set(&search->cfg, CFGKEY_CURRENT_CLIENT, NULL);
            set_current(NULL);
            goto reply;

          error:
            mesg.status = HMESG_STATUS_FAIL;
            if (search)
                mesg.data.string = search->errmsg;
            else
                mesg.data.string = "No matching search in session-core";

          reply:
            // Swap the source and destination fields to reply.
            mesg.src  ^= mesg.dest;
            mesg.dest ^= mesg.src;
            mesg.src  ^= mesg.dest;

            if (mesg_send(STDIN_FILENO, &mesg) < 1)
                fprintf(stderr, "%s: Error sending reply: %s\n",
                        argv[0], mesg.data.string);

            search_cfg = NULL;
        }

        // Generate more points to test.
        do {
            pollstate = NULL;
            for (int i = 0; i < slist_cap; ++i) {
                hsearch_t* search = slist[i];

                if (!search || !search->open)
                    continue;

                // Generate a single trial for this search.
                set_current(search);
                generate_trial(search);
                set_current(NULL);

                if (!pollstate && search->flow.status != HFLOW_WAIT &&
                    search->pending_len < search->pending_cap)
                {
                    pollstate = &polltime;
                }
            }
        } while (pollstate);
    }

    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i]) {
            if (slist[i]->open)
                close_search(slist[i]);
            free_search(slist[i]);
        }
    }
    free(slist);
    free(cbs);
    hmesg_fini(&mesg);

    return retval;
}

/*
 * Base search structure management implementation.
 */

int open_search(hsearch_t* search, hmesg_t* mesg)
{
    const char* ptr;

    // Initialize the session.
    if (hspace_copy(&search->space, mesg->state.space) != 0) {
        search->errmsg = "Could not copy session data";
        return -1;
    }

    if (hcfg_copy(&search->cfg, mesg->data.cfg) != 0) {
        search->errmsg = "Could not copy configuration environment";
        return -1;
    }

    // Control the random number generator state.
    long seed = hcfg_int(&search->cfg, CFGKEY_RANDOM_SEED);
    if (seed == -1) {
        // Default seed is based on time and hsearch_t object address.
        seed = time(NULL);

        search->xsubi[0] = (((unsigned long) seed) >> 16) & 0xFFFF;
        search->xsubi[1] = (((unsigned long) seed)      ) & 0xFFFF;
        search->xsubi[2] = (((unsigned long) search)    ) & 0xFFFF;
    }
    else {
        search->xsubi[0] = (((unsigned long) seed) >> 32) & 0xFFFF;
        search->xsubi[1] = (((unsigned long) seed) >> 16) & 0xFFFF;
        search->xsubi[2] = (((unsigned long) seed)      ) & 0xFFFF;
    }

    // Overwrite CFGKEY_HARMONY_HOME.
    if (hcfg_set(&search->cfg, CFGKEY_HARMONY_HOME, home_dir) != 0) {
        search->errmsg = "Could not set " CFGKEY_HARMONY_HOME " for session";
        return -1;
    }

    if (hcfg_int(&search->cfg, CFGKEY_PERF_COUNT) < 1) {
        search->errmsg = "Invalid " CFGKEY_PERF_COUNT " configuration value";
        return -1;
    }

    search->per_client = hcfg_int(&search->cfg, CFGKEY_GEN_COUNT);
    if (search->per_client < 1) {
        search->errmsg = "Invalid " CFGKEY_GEN_COUNT " configuration value";
        return -1;
    }

    int expected = hcfg_int(&search->cfg, CFGKEY_CLIENT_COUNT);
    if (expected < 1) {
        search->errmsg = "Invalid " CFGKEY_CLIENT_COUNT " configuration value";
        return -1;
    }

    if (extend_lists(search, expected * search->per_client) != 0)
        return -1;

    // Any existing points in the pending list must be re-initialized.
    search->pending_len = 0;
    for (int i = 0; i < search->pending_cap; ++i)
        ((hpoint_t*) &search->pending[i].point)->id = 0;

    // Any index values in the ready list must be re-initialized.
    search->ready_head = 0;
    search->ready_tail = 0;
    for (int i = 0; i < search->ready_cap; ++i)
        search->ready[i] = -1;

    ptr = hcfg_get(&search->cfg, CFGKEY_STRATEGY);
    if (!ptr) {
        if (expected > 1)
            ptr = "pro.so"; // Default strategy for multi-client searches.
        else
            ptr = "nm.so";  // Default strategy for single-client searches.
    }

    // Begin region where plug-in hooks may be executed.  Make sure the
    // correct data is available in the global pointers.
    //
    search_cfg = &search->cfg;

    // Load and initialize requested strategy plug-in.
    if (load_strategy(search, ptr) != 0)
        return -1;

    // Load and initialize requested layer plug-ins.
    ptr = hcfg_get(&search->cfg, CFGKEY_LAYERS);
    if (ptr) {
        if (load_layers(search, ptr) != 0)
            return -1;
    }

    // End region where plug-in hooks may be executed.
    //
    search_cfg = NULL;

    search->best.id = 0;
    search->paused_id = 0;
    search->clients = 1;
    search->open = 1;

    set_current(search);
    for (int i = 0; i < search->pstack_len; ++i) {
        hplugin_t* plugin = &search->pstack[i].plugin;

        if (hplugin_init(plugin, &search->space) != 0)
            return -1;
    }
    set_current(NULL);

    return 0;
}

void close_search(hsearch_t* search)
{
    for (int i = search->pstack_len - 1; i >= 0; --i) {
        hplugin_t* plugin = &search->pstack[i].plugin;

        if (hplugin_fini(plugin) != 0)
            fprintf(stderr, "Error finalizing plug-in: %s\n", search->errmsg);

        if (hplugin_close(plugin, &search->errmsg) != 0)
            fprintf(stderr, "Error closing plug-in: %s\n", search->errmsg);
    }
    search->open = 0;
}

void free_search(hsearch_t* search)
{
    if (!search)
        return;

    for (int i = 0; i < search->pending_cap; ++i) {
        hpoint_fini((hpoint_t*) &search->pending[i].point);
        hperf_fini(&search->pending[i].perf);
    }

    hpoint_fini(&search->flow.point);
    hpoint_fini(&search->best);
    hcfg_fini(&search->cfg);
    hspace_fini(&search->space);

    free(search->buf);
    free(search->ready);
    free(search->pending);
    free(search->pstack);
    free(search);
}

/*
 * Loads strategy given name of library file.
 * Checks that strategy defines required functions,
 * and then calls the strategy's init function (if defined)
 */
int load_strategy(hsearch_t* search, const char* file)
{
    char** pathbuf     = &search->buf;
    int*   pathbuf_len = &search->buf_len;

    search->pstack_len = 0;

    if (snprintf_grow(pathbuf, pathbuf_len,
                      "%s/libexec/%s", home_dir, file) < 0)
    {
        search->errmsg = "Could not allocate full pathname for plug-in";
        return -1;
    }

    return load_plugin(search, *pathbuf, HPLUGIN_STRATEGY);
}

#define LAYER_DELIM_CHARS ":;, \t\n"
int load_layers(hsearch_t* search, const char* list)
{
    char** pathbuf = &search->buf;
    int*   pathlen = &search->buf_len;

    while (*list) {
        // Skip any delimiter characters.
        list += strspn(list, LAYER_DELIM_CHARS);

        // Find the end of the current filename.
        int span = strcspn(list, LAYER_DELIM_CHARS);

        if (snprintf_grow(pathbuf, pathlen,
                          "%s/libexec/%.*s", home_dir, span, list) < 0)
        {
            search->errmsg = "Could not allocate full pathname for plug-in";
            return -1;
        }

        if (load_plugin(search, *pathbuf, HPLUGIN_LAYER) != 0)
            return -1;

        list += span;
    }
    return 0;
}

int load_plugin(hsearch_t* search, const char* pathname, hplugin_type_t type)
{
    if (search->pstack_len == search->pstack_cap) {
        if (array_grow(&search->pstack, &search->pstack_cap,
                       sizeof(*search->pstack)) != 0)
        {
            search->errmsg = "Could not extend plug-in stack";
            return -1;
        }
    }

    hplugin_t* plugin = &search->pstack[ search->pstack_len ].plugin;
    if (hplugin_open(plugin, search->buf, &search->errmsg) != 0)
        return -1;

    // At this point, resources have been allocated.  Increment
    // the counter now so, even if this function encounters an
    // error, these resources may be properly freed during
    // close_search().
    //
    ++search->pstack_len;

    if (plugin->type != type) {
        search->errmsg = "Inappropriate Harmony plug-in type loaded";
        return -1;
    }

    if (plugin->keyinfo) {
        if (hcfg_reginfo(&search->cfg, plugin->keyinfo) != 0) {
            search->errmsg = "Could not register default configuration";
            return -1;
        }
    }

    return 0;
}

/*
 * Internal helper function implementation.
 */

int generate_trial(hsearch_t* search)
{
    int idx;
    htrial_t* trial = NULL;

    // Find a free point.
    for (idx = 0; idx < search->pending_cap; ++idx) {
        trial = &search->pending[idx];
        if (!trial->point.id)
            break;
    }
    if (idx == search->pending_cap) {
        search->errmsg = "Point generation overflow";
        return -1;
    }

    // Reset the performance for this trial.
    hperf_reset(&trial->perf);

    // Call strategy generation routine.
    hplugin_t* strategy = &search->pstack[0].plugin;
    if (hplugin_generate(strategy, &search->flow, trial) != 0)
        return -1;

    if (search->flow.status == HFLOW_ACCEPT) {
        ++search->pending_len;

        // Begin generation workflow for new point.
        search->curr_layer = 1;
        return plugin_workflow(search, idx);
    }
    return 0;
}

int plugin_workflow(hsearch_t* search, int trial_idx)
{
    htrial_t* trial = &search->pending[trial_idx];

    while (search->curr_layer != 0 && search->curr_layer < search->pstack_len)
    {
        int        stack_idx = abs(search->curr_layer);
        hplugin_t* plugin    = &search->pstack[stack_idx].plugin;

        search->flow.status = HFLOW_ACCEPT;
        if (search->curr_layer < 0) {
            // Analyze workflow.
            if (hplugin_analyze(plugin, &search->flow, trial) != 0)
                return -1;
        }
        else {
            // Generate workflow.
            if (hplugin_generate(plugin, &search->flow, trial) != 0)
                return -1;
        }

        int retval = workflow_transition(search, trial_idx);
        if (retval < 0) return -1;
        if (retval > 0) return  0;
    }

    if (search->curr_layer == 0) {
        // Completed analysis layers.  Send trial to strategy.
        hplugin_t* strategy = &search->pstack[0].plugin;
        if (hplugin_analyze(strategy, NULL, trial) != 0)
            return -1;

        // Remove point data from pending list.
        ((hpoint_t*) &trial->point)->id = 0;
        --search->pending_len;

        // Point generation attempts may begin again.
        pollstate = &polltime;
    }
    else if (search->curr_layer == search->pstack_len) {
        // Completed generation layers.  Enqueue trial in ready queue.
        if (search->ready[search->ready_tail] != -1) {
            search->errmsg = "Ready queue overflow";
            return -1;
        }

        search->ready[search->ready_tail] = trial_idx;
        search->ready_tail = (search->ready_tail + 1) % search->ready_cap;
    }
    else {
        search->errmsg = "Invalid current plug-in layer";
        return -1;
    }

    return 0;
}

// Layer state machine transitions.
int workflow_transition(hsearch_t* search, int trial_idx)
{
    switch (search->flow.status) {
    case HFLOW_ACCEPT:
        search->curr_layer += 1;
        break;

    case HFLOW_WAIT:
        if (handle_wait(search, trial_idx) != 0) {
            return -1;
        }
        return 1;

    case HFLOW_RETURN:
    case HFLOW_RETRY:
        search->curr_layer = -search->curr_layer;
        break;

    case HFLOW_REJECT:
        if (handle_reject(search, trial_idx) != 0)
            return -1;

        if (search->flow.status == HFLOW_WAIT)
            return 1;

        search->curr_layer = 1;
        break;

    default:
        search->errmsg = "Invalid workflow status";
        return -1;
    }
    return 0;
}

int handle_reject(hsearch_t* search, int trial_idx)
{
    if (search->curr_layer < 0) {
        search->errmsg = "REJECT invalid for analysis workflow";
        return -1;
    }

    // Regenerate this rejected point.
    htrial_t*  trial    = &search->pending[trial_idx];
    hplugin_t* strategy = &search->pstack[0].plugin;
    return hplugin_rejected(strategy, &search->flow, trial);
}

int handle_wait(hsearch_t* search, int trial_idx)
{
    int   idx = abs(search->curr_layer);
    int** list;
    int*  len;
    int*  cap;

    if (search->curr_layer < 0) {
        list = &search->pstack[idx].wait_analyze;
        len  = &search->pstack[idx].wait_analyze_len;
        cap  = &search->pstack[idx].wait_analyze_cap;
    }
    else {
        list = &search->pstack[idx].wait_generate;
        len  = &search->pstack[idx].wait_generate_len;
        cap  = &search->pstack[idx].wait_generate_cap;
    }

    if (*len == *cap) {
        int i;

        array_grow(list, cap, sizeof(int));
        for (i = *len; i < *cap; ++i)
            (*list)[i] = -1;
    }

    if ((*list)[*len] != -1) {
        search->errmsg = "Could not append to wait list";
        return -1;
    }

    (*list)[*len] = trial_idx;
    ++(*len);
    return 0;
}

hsearch_t* find_search(hmesg_t* mesg)
{
    int idx;

    if (mesg->type == HMESG_SESSION || mesg->type == HMESG_JOIN) {
        int open_slot = slist_cap;

        // Find an active search by name.
        for (idx = 0; idx < slist_cap; ++idx) {
            if (!slist[idx] || !slist[idx]->open) {
                // Save the first open search slot we find.
                if (open_slot == slist_cap)
                    open_slot = idx;
                continue;
            }

            const char* name;
            if (mesg->type == HMESG_SESSION)
                name = mesg->state.space->name;
            else
                name = mesg->data.string;

            if (strcmp(name, slist[idx]->space.name) == 0)
                return slist[idx];
        }
        if (open_slot == slist_cap) {
            if (array_grow(&slist, &slist_cap, sizeof(*slist)) != 0)
                return NULL;
        }

        if (!slist[open_slot])
            slist[open_slot] = calloc(1, sizeof(*slist[open_slot]));
        return slist[open_slot];
    }
    else {
        // Sanity check that the given index is valid and active.
        idx = mesg->dest;
        if (idx >= 0 && idx < slist_cap && slist[idx] && slist[idx]->open)
            return slist[idx];
    }
    return NULL;
}

int handle_callback(callback_t* cb)
{
    hsearch_t* search = cb->search;
    htrial_t** trial_list;
    int* list;
    int* len;
    int  trial_idx, retval;

    search->curr_layer = cb->layer_idx;
    int        idx    = abs(search->curr_layer);
    pstack_t*  pstack = &search->pstack[idx];

    // The idx variable represents layer plug-in index for now.
    if (search->curr_layer < 0) {
        list = pstack->wait_analyze;
        len = &pstack->wait_analyze_len;
    }
    else {
        list = pstack->wait_generate;
        len = &pstack->wait_generate_len;
    }

    if (*len < 1) {
        search->errmsg = "Callback on layer with empty waitlist";
        return -1;
    }

    // Prepare a list of htrial_t pointers.
    trial_list = malloc(*len * sizeof(htrial_t*));
    for (int i = 0; i < *len; ++i)
        trial_list[i] = &search->pending[ list[i] ];

    // Reusing idx to represent waitlist index.  (Shame on me.)
    idx = cb->func(cb->fd, cb->data, &search->flow, *len, trial_list);
    free(trial_list);

    trial_idx = list[idx];
    retval = workflow_transition(search, trial_idx);
    if (retval < 0) return -1;
    if (retval > 0) return  0;

    --(*len);
    list[ idx] = list[*len];
    list[*len] = -1;
    return plugin_workflow(search, trial_idx);
}

int handle_session(hsearch_t* search, hmesg_t* mesg)
{
    if (search->open) {
        search->errmsg = "Search name already exists";
        return -1;
    }

    // Initialize a new session slot.
    if (open_search(search, mesg) != 0)
        return -1;

    // Find the search index.
    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i] == search) {
            mesg->dest = i;
            break;
        }
    }

    mesg->status = HMESG_STATUS_OK;
    return 0;
}

int handle_join(hsearch_t* search, hmesg_t* mesg)
{
    if (!search->open) {
        search->errmsg = "Could not find existing search to join";
        return -1;
    }

    // Grow the pending and ready queues.
    ++search->clients;
    if (extend_lists(search, search->clients * search->per_client) != 0)
        return -1;

    // Launch all join hooks defined in the plug-in stack.
    for (int i = 0; i < search->pstack_len; ++i) {
        hplugin_t* plugin = &search->pstack[i].plugin;

        if (hplugin_join(plugin, mesg->state.client) != 0)
            return -1;
    }

    // Find the search index.
    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i] == search) {
            mesg->dest = i;
            break;
        }
    }

    mesg->state.space = &search->space;
    mesg->status = HMESG_STATUS_OK;
    return 0;
}

int handle_getcfg(hsearch_t* search, hmesg_t* mesg)
{
    // Prepare getcfg response message for client.
    const char* cfgval = hcfg_get(&search->cfg, mesg->data.string);
    if (cfgval == NULL)
        cfgval = "";

    snprintf_grow(&search->buf, &search->buf_len, "%s=%s",
                  mesg->data.string, cfgval);
    mesg->data.string = search->buf;
    mesg->status = HMESG_STATUS_OK;
    return 0;
}

int handle_setcfg(hsearch_t* search, hmesg_t* mesg)
{
    char* sep = (char*) strchr(mesg->data.string, '=');
    const char* oldval;

    if (!sep) {
        search->errmsg = "Malformed configuration request";
        return -1;
    }
    *sep = '\0';

    // Store the original value, possibly allocating memory for it.
    oldval = hcfg_get(&search->cfg, mesg->data.string);
    if (oldval) {
        snprintf_grow(&search->buf, &search->buf_len, "%s", oldval);
        oldval = search->buf;
    }

    if (search_setcfg(mesg->data.string, sep + 1) != 0) {
        search->errmsg = "Error setting session configuration variable";
        return -1;
    }

    // Prepare setcfg response message for client.
    mesg->data.string = oldval;
    mesg->status = HMESG_STATUS_OK;
    return 0;
}

int handle_best(hmesg_t* mesg)
{
    mesg->status = HMESG_STATUS_OK;
    return 0;
}

int handle_fetch(hsearch_t* search, hmesg_t* mesg)
{
    int idx = search->ready[search->ready_head], paused;

    // Check if the session is paused.
    paused = hcfg_bool(&search->cfg, CFGKEY_PAUSED);

    if (!paused && idx >= 0) {
        // Send the next point on the ready queue.
        mesg->data.point = &search->pending[idx].point;

        // Remove the first point from the ready queue.
        search->ready[search->ready_head] = -1;
        search->ready_head = (search->ready_head + 1) % search->ready_cap;

        mesg->status = HMESG_STATUS_OK;
    }
    else {
        // Ready queue is empty, or session is paused.
        // Send the best known point.
        hplugin_t* strategy = &search->pstack[0].plugin;
        if (hplugin_best(strategy, &search->best) != 0)
            return -1;

        mesg->state.best = &search->best;
        search->paused_id = mesg->data.point->id;
        mesg->status = HMESG_STATUS_BUSY;
    }
    return 0;
}

int handle_report(hsearch_t* search, hmesg_t* mesg)
{
    int idx;
    htrial_t* trial = NULL;

    // Find the associated trial in the pending list.
    for (idx = 0; idx < search->pending_cap; ++idx) {
        trial = &search->pending[idx];
        if (trial->point.id == mesg->data.point->id)
            break;
    }
    if (idx == search->pending_cap) {
        if (mesg->data.point->id == search->paused_id) {
            mesg->status = HMESG_STATUS_OK;
            return 0;
        }
        else {
            search->errmsg = "Rouge point support not yet implemented";
            return -1;
        }
    }
    search->paused_id = 0;

    // Update performance in our local records.
    hperf_copy(&trial->perf, mesg->data.perf);

    // Begin the workflow at the outermost analysis layer.
    search->curr_layer = -search->pstack_len + 1;
    if (plugin_workflow(search, idx) != 0)
        return -1;

    mesg->status = HMESG_STATUS_OK;
    return 0;
}

int handle_command(hsearch_t* search, hmesg_t* mesg)
{
    if (strcmp(mesg->data.string, "restart") == 0) {
        if (search_restart() != 0)
            goto error;
    }
    else if (strcmp(mesg->data.string, "leave") == 0) {
        --search->clients;
    }
    else if (strcmp(mesg->data.string, "kill") == 0) {
        close_search(search);
    }
    else {
        search->errmsg = "Received unknown command";
        goto error;
    }

    mesg->status = HMESG_STATUS_OK;
    return 0;

  error:
    mesg->status = HMESG_STATUS_FAIL;
    return -1;
}

int update_state(hmesg_t* mesg, hsearch_t* search)
{
    // Send updated search space to client if necessary.
    if (mesg->state.space->id < search->space.id) {
        mesg->state.space = &search->space;
    }
    else if (mesg->type != HMESG_JOIN) {
        mesg->state.space = &mesg->unpacked_space;
        mesg->unpacked_space.id = 0;
    }

    // Refresh the best known point from the strategy.
    if (search->open) {
        hplugin_t* strategy = &search->pstack[0].plugin;
        if (hplugin_best(strategy, &search->best) != 0)
            return -1;
    }

    // Send it to the client if necessary.
    if (!mesg->state.best || mesg->state.best->id < search->best.id) {
        mesg->state.best = &search->best;
    }
    else {
        mesg->state.best = &mesg->unpacked_best;
        mesg->unpacked_best.id = 0;
    }
    return 0;
}


int extend_lists(hsearch_t* search, int target_cap)
{
    int orig_cap = search->pending_cap;

    if (orig_cap >= target_cap && orig_cap != 0)
        return 0;

    if (   search->ready
        && search->ready_tail <= search->ready_head
        && search->ready[search->ready_head] != -1)
    {
        int i = search->ready_cap - search->ready_head;

        // Shift ready array to align head with array index 0.
        reverse_array(search->ready, 0, search->ready_cap);
        reverse_array(search->ready, 0, i);
        reverse_array(search->ready, i, search->ready_cap);

        search->ready_head = 0;
        search->ready_tail = (search->ready_cap
                              - search->ready_tail
                              + search->ready_head);
    }

    search->ready_cap = target_cap;
    if (array_grow(&search->ready, &search->ready_cap,
                   sizeof(*search->ready)) != 0)
    {
        search->errmsg = "Could not extend ready array";
        return -1;
    }

    search->pending_cap = target_cap;
    if (array_grow(&search->pending, &search->pending_cap,
                   sizeof(*search->pending)) != 0)
    {
        search->errmsg = "Could not extend pending array";
        return -1;
    }

    for (int i = orig_cap; i < search->pending_cap; ++i)
        search->ready[i] = -1;

    return 0;
}

void reverse_array(void* ptr, int head, int tail)
{
    unsigned long* arr = (unsigned long*) ptr;

    while (head < --tail) {
        // Swap head and tail entries.
        arr[head] ^= arr[tail];
        arr[tail] ^= arr[head];
        arr[head] ^= arr[tail];
        ++head;
    }
}

void set_current(hsearch_t* search)
{
    current_search = search;
    search_cfg     = &search->cfg;
}

/*
 * Exported functions for pluggable modules.
 *
 * Since these function will be called from within plug-in module
 * hooks, there will exist the notion of a "current search instance."
 * These functions rely on the current search pointer to be set
 * correctly, and access search data exclusively through it.
 */

/*
 * Retrieve the best known configuration, as determined by the strategy.
 */
int search_best(hpoint_t* pt)
{
    hplugin_t* strategy = &current_search->pstack[0].plugin;
    return hplugin_best(strategy, pt);
}

/*
 * Register an "analyze" callback to handle file descriptor data.
 */
int search_callback_analyze(int fd, void* data, cb_func_t func)
{
    if (cbs_len >= cbs_cap) {
        if (array_grow(&cbs, &cbs_cap, sizeof(*cbs)) != 0)
            return -1;
    }

    cbs[ cbs_len ].fd        = fd;
    cbs[ cbs_len ].search    = current_search;
    cbs[ cbs_len ].layer_idx = -current_search->curr_layer;
    cbs[ cbs_len ].data      = data;
    cbs[ cbs_len ].func      = func;
    ++cbs_len;

    FD_SET(fd, &fds);
    if (maxfd < fd)
        maxfd = fd;

    return 0;
}

/*
 * Register a "generate" callback to handle file descriptor data.
 */
int search_callback_generate(int fd, void* data, cb_func_t func)
{
    if (cbs_len >= cbs_cap) {
        if (array_grow(&cbs, &cbs_cap, sizeof(*cbs)) != 0)
            return -1;
    }

    cbs[ cbs_len ].fd        = fd;
    cbs[ cbs_len ].search    = current_search;
    cbs[ cbs_len ].layer_idx = current_search->curr_layer;
    cbs[ cbs_len ].data      = data;
    cbs[ cbs_len ].func      = func;
    ++cbs_len;

    FD_SET(fd, &fds);
    if (maxfd < fd)
        maxfd = fd;

    return 0;
}

/*
 * Set the error message for the current search instance.
 */
void search_error(const char* msg)
{
    current_search->errmsg = msg;
}

/*
 * Trigger a restart of the current search instance.
 */
int search_restart(void)
{
    // Re-initialize all plug-ins associated with this search.
    for (int i = 0; i < current_search->pstack_len; ++i) {
        hplugin_t* plugin = &current_search->pstack[i].plugin;

        if (hplugin_init(plugin, &current_search->space) != 0)
            return -1;
    }

    return 0;
}

/*
 * Central interface for shared configuration between pluggable modules.
 */
int search_setcfg(const char* key, const char* val)
{
    if (hcfg_set(&current_search->cfg, key, val) != 0)
        return -1;

    // Make sure setcfg callbacks are triggered after the
    // configuration has been set.  Otherwise, any subsequent calls to
    // setcfg further down the stack will be nullified when we return
    // to this frame.
    //
    for (int i = 0; i < current_search->pstack_len; ++i) {
        hplugin_t* plugin = &current_search->pstack[i].plugin;

        if (hplugin_setcfg(plugin, key, val) != 0)
            return -1;
    }

    return 0;
}

/*
 * Return a session-specific pseudo-random double value.
 */
double search_drand48(void)
{
    return erand48(current_search->xsubi);
}

/*
 * Return a session-specific pseudo-random long integer value.
 */
long int search_lrand48(void)
{
    return nrand48(current_search->xsubi);
}
