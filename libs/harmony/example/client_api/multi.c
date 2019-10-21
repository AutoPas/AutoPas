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
#define _XOPEN_SOURCE 600 // Needed for lrand48() and snprintf().

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

#include "hclient.h"
#include "defaults.h"

#define MAX_EVALS 256
#define MAX_RUNNING_TASKS 16
#define MAX_STARTED_TASKS 64

typedef enum task_state {
    TASK_STATE_RUNNING = 0,
    TASK_STATE_LEAVE,
    TASK_STATE_CONVERGED,
    TASK_STATE_ERROR
} task_state_t;

typedef struct sinfo {
    int          id;
    htask_t*     htask;
    int          count;
    long         ival;
    double       rval;
    const char*  sval;
    task_state_t state;
} sinfo_t;

const char* fruits[] = {"apples",
                        "bananas",
                        "cherries",
                        "figs",
                        "grapes",
                        "oranges",
                        "peaches",
                        "pineapple",
                        "pears",
                        "watermelon",
                        "strawberries",
                        NULL};

/*
 * A simple performance function is defined here for illustration purposes.
 */
double application(long ival, double rval, const char* string)
{
    int i;
    double sval = 0.0;
    for (i = 0; string[i]; ++i)
        sval += string[i];
    return sval * ival / rval;
}

hdef_t* define_search(void)
{
    hdef_t* hdef = ah_def_alloc();
    if (!hdef) {
        fprintf(stderr, "Error allocating a definition descriptor");
        goto error;
    }

    // Define a tuning variable that resides in the integer domain.
    if (ah_def_int(hdef, "i_var",  1, 1000, 1, NULL) != 0) {
        fprintf(stderr, "Error defining an integer tuning variable");
        goto error;
    }

    // Define a tuning variable that resides in the real domain.
    if (ah_def_real(hdef, "r_var", 0.0001, 1.0, 0.0001, NULL) != 0) {
        fprintf(stderr, "Error defining a real tuning variable");
        goto error;
    }

    for (int i = 0; fruits[i]; ++i) {
        if (ah_def_enum_value(hdef, "s_var", fruits[i]) != 0) {
            fprintf(stderr, "Error defining an enumerated tuning variable");
            goto error;
        }
    }
    return hdef;

  error:
    ah_def_free(hdef);
    return NULL;
}

static const char* strategy[] = {
    "exhaustive.so",
    "random.so",
    "nm.so",
    "pro.so"
};

int start_search(hdesc_t* hdesc, hdef_t* hdef, sinfo_t* sinfo, int id)
{
    char namebuf[32];

    // Tabula rasa.
    memset(sinfo, 0, sizeof(*sinfo));

    // Provide a unique name to the search.
    snprintf(namebuf, sizeof(namebuf), "multi_%03d", id);
    if (ah_def_name(hdef, namebuf) != 0) {
        fprintf(stderr, "Error setting search name to %s", namebuf);
        return -1;
    }

    // Provide a different strategy for each search.
    if (ah_def_strategy(hdef, strategy[id % 4]) != 0) {
        fprintf(stderr, "Error setting search strategy to %s",
                strategy[id % 4]);
        return -1;
    }

    // Start a new search in the session.
    sinfo->htask = ah_start(hdesc, hdef);
    if (!sinfo->htask) {
        fprintf(stderr, "Error starting search %d", id);
        return -1;
    }

    if (ah_bind_int(sinfo->htask, "i_var", &sinfo->ival) != 0) {
        fprintf(stderr, "Error binding 'i_var' for search %d", id);
        return -1;
    }

    if (ah_bind_real(sinfo->htask, "r_var", &sinfo->rval) != 0) {
        fprintf(stderr, "Error binding 'r_var' for search %d", id);
        return -1;
    }

    if (ah_bind_enum(sinfo->htask, "s_var", &sinfo->sval) != 0) {
        fprintf(stderr, "Error binding 's_var' for search %d", id);
        return -1;
    }

    sinfo->id = id;
    return 0;
}

int end_search(sinfo_t* sinfo)
{
    int retval = 0;

    if (!sinfo->htask)
        return 0;

    if (sinfo->state == TASK_STATE_ERROR) {
        printf("%4d: Error after %4d evals: %s.\n",
               sinfo->id, sinfo->count, ah_error());

        ah_kill(sinfo->htask);
        ah_error_clear();
    }
    else {
        if (sinfo->state == TASK_STATE_CONVERGED) {
            printf("%4d: Converged after %4d evals.", sinfo->id, sinfo->count);
        }
        else if (sinfo->state == TASK_STATE_LEAVE) {
            printf("%4d: Left after %4d evals.", sinfo->id, sinfo->count);
        }

        if (ah_best(sinfo->htask) != 0) {
            fprintf(stderr, "\nError retrieving best search");
            goto error;
        }

        printf(" Best: %4ld, %.4f, %s\n",
               sinfo->ival, sinfo->rval, sinfo->sval);

        if (ah_kill(sinfo->htask) != 0) {
            fprintf(stderr, "Error killing search #%d", sinfo->id);
            goto error;
        }
    }
    goto cleanup;

  error:
    retval = -1;

  cleanup:
    sinfo->htask = NULL;
    return retval;
}

void shuffle(int* order)
{
    for (int i = 0; i < MAX_RUNNING_TASKS; ++i) {
        int j = lrand48() % (i + 1);
        order[i] = order[j];
        order[j] = i;
    }
}

int main(int argc, char* argv[])
{
    hdef_t* hdef = NULL;
    int retval = 0;
    sinfo_t slist[MAX_RUNNING_TASKS] = {{0}};
    int order[MAX_RUNNING_TASKS];
    int started_tasks = 0;
    int done;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            fprintf(stderr, "Usage: %s"
                    " [KEY_1=VAL_1] ... [KEY_N=VAL_N]\n\n", argv[0]);
            return 0;
        }
    }

    // Initialize a Harmony client.
    hdesc_t* hdesc = ah_alloc();
    if (hdesc == NULL) {
        fprintf(stderr, "Error initializing a Harmony descriptor");
        goto cleanup;
    }
    ah_args(hdesc, &argc, argv);

    printf("Starting Harmony...\n");
    if (ah_connect(hdesc, NULL, 0) != 0) {
        fprintf(stderr, "Error connecting to Harmony session");
        goto error;
    }

    // Define a new tuning search.
    hdef = define_search();
    if (!hdef)
        goto error;

    // Begin a new tuning search.
    for (int i = 0; i < MAX_RUNNING_TASKS; ++i) {
        if (start_search(hdesc, hdef, &slist[i], ++started_tasks) != 0)
            goto error;
    }

    // Main tuning loop.
    done = 0;
    while (!done) {
        // Fetch new values from the search task in a random order.
        shuffle(order);
        for (int i = 0; i < MAX_RUNNING_TASKS; ++i) {
            int idx = order[i];
            if (!slist[idx].htask)
                continue;

            sinfo_t* sinfo = &slist[idx];
            int hresult = ah_fetch(sinfo->htask);
            if (hresult < 0) {
                sinfo->state = TASK_STATE_ERROR;
                if (end_search(sinfo) != 0)
                    goto error;
            }
            else if (hresult == 0) {
                fprintf(stderr, "Waiting for point on #%d.\n", sinfo->id);
                sleep(1);
                --i;
                continue;
            }
        }

        // Report new values to the search task in a different random order.
        shuffle(order);
        for (int i = 0; i < MAX_RUNNING_TASKS; ++i) {
            int idx = order[i];
            if (!slist[idx].htask)
                continue;

            sinfo_t* sinfo = &slist[idx];
            double perf = application(sinfo->ival, sinfo->rval, sinfo->sval);

            if (sinfo->state != TASK_STATE_ERROR) {
                if (ah_report(sinfo->htask, &perf) != 0) {
                    sinfo->state = TASK_STATE_ERROR;
                }
                else if (ah_converged(sinfo->htask)) {
                    sinfo->state = TASK_STATE_CONVERGED;
                }
                else if (++sinfo->count >= MAX_EVALS) {
                    sinfo->state = TASK_STATE_LEAVE;
                }
            }
        }

        for (int i = 0; i < MAX_RUNNING_TASKS; ++i) {
            sinfo_t* sinfo = &slist[i];

            if (sinfo->state != TASK_STATE_RUNNING) {
                if (end_search(sinfo) != 0)
                    goto error;

                if (started_tasks < MAX_STARTED_TASKS) {
                    if (start_search(hdesc, hdef, sinfo, ++started_tasks) != 0)
                        goto error;
                }
            }
        }

        done = 1;
        for (int i = 0; i < MAX_RUNNING_TASKS; ++i) {
            if (slist[i].state == TASK_STATE_RUNNING)
                done = 0;
        }
    }
    goto cleanup;

  error:
    fprintf(stderr, ": %s\n", ah_error());
    retval = -1;

  cleanup:
    // Close the connection to the tuning session.
    if (ah_close(hdesc) != 0) {
        fprintf(stderr, "Error disconnecting from Harmony session: %s.\n",
                ah_error());
    }
    ah_def_free(hdef);
    ah_free(hdesc);
    return retval;
}
