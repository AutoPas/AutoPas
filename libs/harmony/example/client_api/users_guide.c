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
 * This is the example code included in our User's Guide.
 *
 * It provides sample code for connecting to a session, establishing a
 * new search task, and interacting with that task.
 *
 * The code in this file does not check the return values of API
 * calls.  This allows us to favor clarity over functionality.  We do
 * not recommended this practice for production code.  The source file
 * exists merely to ensure that our example compiles and runs without
 * error.
 *
 * See our other example programs for proper usage of our API.
 */

/** [snippet_session_connect] */
#include <stdio.h>
#include "hclient.h"

hdesc_t* init_tuner(int* argc, char** argv)
{
    hdesc_t* hdesc;

    /* Initialize a Harmony tuning session handle. */
    hdesc = ah_alloc();

    /* Read any Harmony configuration directives passed on the
     * command line.  This step is optional (but often helpful!).
     */
    ah_args(hdesc, argc, argv);

    /* Connect to a Harmony session.  Since we're passing NULL as
     * our hostname, the HARMONY_HOST environment variable will be
     * used if defined.  Otherwise, a local session will be spawned.
     */
    ah_connect(hdesc, NULL, 0);

    return hdesc;
}
/** [snippet_session_connect] */

void advanced_config(hdef_t* hdef);

/** [snippet_task_def] */
#include <stdio.h>
#include "hclient.h"

htask_t* new_search(hdesc_t* hdesc)
{
    hdef_t*  hdef;
    htask_t* htask;

    /* Allocate a new Harmony search definition handle. */
    hdef = ah_def_alloc();

    /* Give the tuning session a unique name. */
    ah_def_name(hdef, "Example");

    /* Add an integer variable called "intVar1" to the session.
     * Its value may range between 1 and 100 (inclusive).
     */
    ah_def_int(hdef, "intVar1", 1, 100, 1, NULL);

    /* Add another integer variable called "intVar2" to the session.
     * Its value may range between 2 and 200 (inclusive) by
     * strides of 2.
     */
    ah_def_int(hdef, "intVar2", 2, 200, 2, NULL);

    /* Add a real-valued variable called "realVar" to the session.
     * Its value may range between 0.0 and 0.5 (inclusive), using the
     * full precision available by an IEEE double.
     */
    ah_def_real(hdef, "realVar", 0.0, 0.5, 0.0, NULL);

    /* Add a string-valued variable called "strVar" to the session. */
    ah_def_enum(hdef, "strVar", NULL);

    /* The variable "strVar" may be "apples", "oranges", "peaches",
     * or "pears".
     */
    ah_def_enum_value(hdef, "strVar", "apples");
    ah_def_enum_value(hdef, "strVar", "oranges");
    ah_def_enum_value(hdef, "strVar", "peaches");
    ah_def_enum_value(hdef, "strVar", "pears");

    /* Here we pass the definition handle to a subroutine to further
     * define the search task.  See the "Advanced Search Task
     * Configuration" example for more details.
     */
    advanced_config(hdef);

    /* Start a new tuning task in the Harmony session.  It will be
     * based on the search we've described above.
     */
    htask = ah_start(hdesc, hdef);

    /* Don't forget to free the definition handle! */
    ah_def_free(hdef);

    return htask;
}
/** [snippet_task_def] */

/** [snippet_task_def_advanced] */
#include <stdio.h>
#include "hclient.h"

void advanced_config(hdef_t* hdef)
{
    /* Use the Parallel Rank Order simplex-based as the strategy for this
     * session, instead of the default Nelder-Mead.
     */
    ah_def_strategy(hdef, "pro.so");

    /* Instruct the strategy to use an initial simplex roughly half
     * the size of the search space.
     */
    ah_def_cfg(hdef, "INIT_RADIUS", "0.50");

    /* This tuning session should surround the search strategy with
     * a logger layer first, and an aggregator layer second.
     */
    ah_def_layers(hdef, "log.so:agg.so");

    /* Instruct the logger to use /tmp/tuning.run as the logfile. */
    ah_def_cfg(hdef, "LOG_FILE", "/tmp/tuning.run");

    /* Instruct the aggregator to collect 10 performance values for
     * each point, and allow the median performance to continue through
     * the feedback loop.
     */
    ah_def_cfg(hdef, "AGG_TIMES", "10");
    ah_def_cfg(hdef, "AGG_FUNC", "median");
}
/** [snippet_task_def_advanced] */

#define work(a, b, c, d) (0.0)

/** [snippet_task_interact] */
#include <stdio.h>
#include "hclient.h"

int main(int argc, char* argv[])
{
    /* Harmony data types. */
    hdesc_t* hdesc;
    htask_t* htask;

    /* Variables for this application's run-time tunable parameters. */
    long        var1;
    long        var2;
    double      var3;
    const char* var4;

    /* These functions are defined in the examples above. */
    hdesc = init_tuner(&argc, argv);
    htask = new_search(hdesc);

    /* Bind local memory to values fetched from Harmony.
     * This allows ah_fetch() to modify local variables automatically.
     */
    ah_bind_int(htask,  "intVar1", &var1);
    ah_bind_int(htask,  "intVar2", &var2);
    ah_bind_real(htask, "realVar", &var3);
    ah_bind_enum(htask, "strVar",  &var4);

    /* Loop until the session has reached a converged state. */
    while (!ah_converged(htask))
    {
        /* Define a variable to hold the resulting performance value. */
        double perf;

        /* Retrieve new values from the Harmony Session. */
        ah_fetch(htask);

        /* The local variables var1, var2, var3, and var4 have now
         * been updated and are ready for use.
         *
         * This is where a normal application would do some work
         * using these variables and measure the performance.
         * Since this is a simple example, we'll pretend the
         * function "work()" will take the variables, and produce a
         * performance value.
         */
        perf = work(var1, var2, var3, var4);

        /* Report the performance back to the session. */
        ah_report(htask, &perf);
    }

    /* Leave the search task. */
    ah_leave(htask);

    /* Free the session descriptor.
     * All associated task handles are also destroyed by this call. */
    ah_free(hdesc);

    return 0;
}
/** [snippet_task_interact] */
