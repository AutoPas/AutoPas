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
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>

#include "hclient.h"
#include "defaults.h"

#define MAX_LOOP 500

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

int main(int argc, char* argv[])
{
    hdesc_t* hdesc;
    hdef_t*  hdef = NULL;
    htask_t* htask;
    int i, retval;
    double perf;
    char* filename = "session.cfg";

    retval = 0;
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            fprintf(stderr, "Usage: %s [filename]"
                    " [KEY_1=VAL_1] ... [KEY_N=VAL_N]\n\n", argv[0]);
            return 0;
        }
    }

    // Initialize a Harmony client.
    hdesc = ah_alloc();
    if (hdesc == NULL) {
        fprintf(stderr, "Error initializing a Harmony session");
        goto cleanup;
    }
    ah_args(hdesc, &argc, argv);

    if (argc > 1)
        filename = argv[1];

    // Load a session definition file.
    hdef = ah_def_load(filename);
    if (!hdef) {
        fprintf(stderr, "Error loading session file");
        goto error;
    }

    printf("Starting Harmony...\n");
    if (ah_connect(hdesc, NULL, 0) != 0) {
        fprintf(stderr, "Error connecting to Harmony session");
        goto error;
    }

    // Begin a new tuning search.
    htask = ah_start(hdesc, hdef);
    if (!htask) {
        fprintf(stderr, "Error starting Harmony search");
        goto error;
    }

    // Main tuning loop.
    for (i = 0; !ah_converged(htask) && i < MAX_LOOP; ++i) {
        int hresult = ah_fetch(htask);
        if (hresult < 0) {
            fprintf(stderr, "Error fetching values from search task");
            goto error;
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
        perf = application(ah_get_int(htask, "i_var"),
                           ah_get_real(htask, "r_var"),
                           ah_get_enum(htask, "fruits"));

        printf("(%4ld, %.4lf, \"%s\") = %lf\n",
               ah_get_int(htask, "i_var"),
               ah_get_real(htask, "r_var"),
               ah_get_enum(htask, "fruits"),
               perf);

        // Report the performance we've just measured.
        if (ah_report(htask, &perf) != 0) {
            fprintf(stderr, "Error reporting performance to search task");
            goto error;
        }
    }

    if (!ah_converged(htask)) {
        printf("*\n");
        printf("* Leaving tuning search after %d iterations.\n", MAX_LOOP);
        printf("*\n");
    }

    if (ah_best(htask) < 0) {
        fprintf(stderr, "Error retrieving best tuning point");
        goto error;
    }
    perf = application(ah_get_int(htask, "i_var"),
                       ah_get_real(htask, "r_var"),
                       ah_get_enum(htask, "fruits"));

    printf("(%4ld, %.4lf, \"%s\") = %lf (* Best point found. *)\n",
           ah_get_int(htask, "i_var"),
           ah_get_real(htask, "r_var"),
           ah_get_enum(htask, "fruits"),
           perf);
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
