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

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <math.h>

#include "hclient.h"

/*
 * For illustration purposes, the performance here is defined by following
 * simple definition:
 *   perf = (p1 - 15)^2 + (p2 - 30)^2 + (p3 - 45)^2 +
 *          (p4 - 60)^2 + (p5 - 75)^2 + (p6 - 90)^2
 *
 * So the theoretical minimum can be found at point:
 *      (15, 30, 45, 60, 75, 90)
 *
 * And a reasonable search range for all parameters is [1-100].
 */
long application(long p1, long p2, long p3, long p4, long p5, long p6)
{
    long perf =
        (p1-15) * (p1-15) +
        (p2-30) * (p2-30) +
        (p3-45) * (p3-45) +
        (p4-60) * (p4-60) +
        (p5-75) * (p5-75) +
        (p6-90) * (p6-90);
    return perf;
}

htask_t* start_search(hdesc_t* hdesc, const char* name)
{
    hdef_t* hdef;
    htask_t* htask;

    hdef = ah_def_alloc();
    if (!hdef) {
        fprintf(stderr, "Error allocating search definition: %s\n",
                ah_error());
        return NULL;
    }

    if (ah_def_name(hdef, name)                      != 0 ||
        ah_def_layers(hdef, "xmlWriter.so")          != 0 ||
        ah_def_int(hdef, "param_1", 1, 100, 1, NULL) != 0 ||
        ah_def_int(hdef, "param_2", 1, 100, 1, NULL) != 0 ||
        ah_def_int(hdef, "param_3", 1, 100, 1, NULL) != 0 ||
        ah_def_int(hdef, "param_4", 1, 100, 1, NULL) != 0 ||
        ah_def_int(hdef, "param_5", 1, 100, 1, NULL) != 0 ||
        ah_def_int(hdef, "param_6", 1, 100, 1, NULL) != 0)
    {
        fprintf(stderr, "Error during session setup: %s\n", ah_error());
        goto error;
    }

    // Begin a new tuning session.
    htask = ah_start(hdesc, hdef);
    if (!htask) {
        fprintf(stderr, "Could not start tuning search task: %s\n",
                ah_error());
        goto error;
    }
    goto cleanup;

  error:
    htask = NULL;

  cleanup:
    ah_def_free(hdef);
    return htask;
}

int main(int argc, char* argv[])
{
    hdesc_t* hdesc;
    htask_t* htask;
    const char* name;
    int i, retval, loop = 200;
    double perf = -HUGE_VAL;

    // Variables to hold the application's runtime tunable parameters.
    // Once bound to a Harmony tuning session, these variables will be
    // modified upon ah_fetch() to a new testing configuration.
    //
    long param_1;
    long param_2;
    long param_3;
    long param_4;
    long param_5;
    long param_6;

    retval = 0;
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            fprintf(stderr, "Usage: %s [session_name]"
                    " [KEY_1=VAL_1] .. [KEY_N=VAL_N]\n\n", argv[0]);
            return 0;
        }
    }

    // Initialize a Harmony client.
    hdesc = ah_alloc();
    if (hdesc == NULL) {
        fprintf(stderr, "Error allocating Harmony descriptor: %s\n",
                ah_error());
        goto error;
    }
    ah_args(hdesc, &argc, argv);

    // Connect to Harmony session.
    if (ah_connect(hdesc, NULL, 0) != 0) {
        fprintf(stderr, "Error connecting to Harmony session: %s\n",
                ah_error());
        goto error;
    }

    // Process the program arguments.
    name = "XML_example";
    if (argc > 1)
        name = argv[1];

    htask = start_search(hdesc, name);
    if (!htask)
        goto error;

    // Bind the session variables to local variables.
    if (ah_bind_int(htask, "param_1", &param_1) != 0 ||
        ah_bind_int(htask, "param_2", &param_2) != 0 ||
        ah_bind_int(htask, "param_3", &param_3) != 0 ||
        ah_bind_int(htask, "param_4", &param_4) != 0 ||
        ah_bind_int(htask, "param_5", &param_5) != 0 ||
        ah_bind_int(htask, "param_6", &param_6) != 0)
    {
        fprintf(stderr, "Error binding tuning variable to local memory: %s\n",
                ah_error());
        goto error;
    }

    // Main loop.
    for (i = 0; !ah_converged(htask) && i < loop; i++) {
        int hresult = ah_fetch(htask);
        if (hresult < 0) {
            fprintf(stderr, "Error fetching values from search task: %s\n",
                    ah_error());
            goto error;
        }
        else if (hresult == 0) {
            // New values were not available at this time.
            // Bundles remain unchanged by Harmony system.
        }
        else if (hresult > 0) {
            // The Harmony system modified the variable values.
            // Do any systemic updates to deal with such a change.
        }

        // Run one full iteration of the application (or code variant).
        //
        // Here our application is rather simple. Definition of performance can
        // be user-defined. Depending on application, it can be MFlops/sec,
        // time to complete the entire run of the application, cache hits vs.
        // misses and so on.
        //
        // For searching the parameter space in a Transformation framework,
        // just run different parameterized code variants here. A simple
        // mapping between the parameters and the code-variants is needed to
        // call the appropriate code variant.
        //
        perf = application(param_1, param_2, param_3,
                           param_4, param_5, param_6);

        if (hresult > 0) {
            // Only print performance if new values were fetched.
            printf("%ld, %ld, %ld, %ld, %ld, %ld = %lf\n",
                   param_1, param_2, param_3,
                   param_4, param_5, param_6, perf);
        }

        // Report the performance we've just measured.
        if (ah_report(htask, &perf) != 0) {
            fprintf(stderr, "Error reporting performance to search task: %s\n",
                    ah_error());
            goto error;
        }
    }


    // Leave the session.
    if (ah_close(hdesc) != 0) {
        fprintf(stderr, "Error disconnecting from Harmony session: %s\n",
                ah_error());
        goto error;
    }
    goto cleanup;

  error:
    retval = -1;

  cleanup:
    ah_free(hdesc);
    return retval;
}
