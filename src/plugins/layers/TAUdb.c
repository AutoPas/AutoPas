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
 * \page TAUdb TAUdb Interface (TAUdb.so)
 *
 * \warning This processing layer is still considered pre-beta.
 *
 * This processing layer uses the TAU Performance System's C API to
 * keep a log of point/performance pairs to disk as they flow through
 * the auto-tuning [feedback loop](\ref intro_feedback).
 *
 * The `LIBTAUDB` [build variable](\ref start_install) must be defined
 * at build time for this plug-in to be available, since it is
 * dependent on a library distributed with TAU.  The distribution of
 * TAU is available here:
 * - http://www.cs.uoregon.edu/research/tau/downloads.php
 *
 * And `LIBTAUDB` should point to the `tools/src/taudb_c_api`
 * directory within that distribution.
 */

#include "hlayer.h"
#include "session-core.h"
#include "hspace.h"
#include "hpoint.h"
#include "hutil.h"
#include "hcfg.h"

#include "taudb_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>

/*
 * Name used to identify this plugin layer.
 * All Harmony plugin layers must define this variable.
 */
const char hplugin_name[] = "TAUdb";

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_TAUDB_NAME, NULL,
      "Name of the PostgreSQL database." },
    { CFGKEY_TAUDB_STORE_METHOD, "one_time",
      "Determines when statistics are computed: real_time: With each "
      "database write. one_time: At session cleanup time." },
    { CFGKEY_TAUDB_STORE_NUM, "0",
      "Number of reports to cache before writing to database." },
    { NULL }
};

#define REG_STR_LEN 32

typedef struct cinfo {
    char* id;
    int timer;
} cinfo_t;

/*
 * Structure to hold all data needed by an individual search instance.
 *
 * To support multiple parallel search instances, no global variables
 * should be defined or used in this plug-in layer.  They should
 * instead be defined as a part of this structure.
 */
struct hplugin_data {
    hspace_t* space;

    TAUDB_METRIC*      metric;
    TAUDB_TRIAL*       taudb_trial;
    TAUDB_CONNECTION*  connection;
    TAUDB_THREAD*      thread;
    TAUDB_TIMER_GROUP* timer_group;

    cinfo_t* client;
    int param_max;
    int client_max;
    int save_counter;
    int taudb_store_type; // 1 for real time, 0 for one-time
    int total_record_num;
};

/*
 * Internal helper function prototypes.
 */
static int client_idx(hplugin_data_t* data);
static int save_timer_parameter(hplugin_data_t* data,
                                TAUDB_TIMER* timer, htrial_t* trial);
static int init_tau_trial(hplugin_data_t* data,
                          const char* appName, const char* trialname);
static TAUDB_THREAD* create_tau_thread(hplugin_data_t* data, int threadNum);
static void get_tau_metadata(hplugin_data_t* data, TAUDB_THREAD* thread,
                             char* opsys, char* machine,
                             char* release, char* hostname,
                             char* procnum, char* cpuvendor,
                             char* cpumodel, char* clkfreq,
                             char* cachesize);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* TAUdb_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search task.
 */
int TAUdb_init(hplugin_data_t* data, hspace_t* space)
{
    char* tmpstr;

    // Connecting to TAUdb.
    tmpstr = hcfg_get(search_cfg, CFGKEY_TAUDB_NAME);
    if (!tmpstr) {
        search_error("TAUdb connection failed: config file not found");
        return -1;
    }

    data->connection = taudb_connect_config(tmpstr);

    // Check if the connection has been established.
    taudb_check_connection(data->connection);

    // Initializing trial.
    if (init_tau_trial(data, space->name, space->name) != 0) {
        search_error("Failed to create TAUdb trial");
        return -1;
    }

    // Create a metric.
    data->metric = taudb_create_metrics(1);
    data->metric->name = taudb_strdup("TIME");
    taudb_add_metric_to_trial(data->taudb_trial, data->metric);

    // Initializing timer group.
    data->timer_group = taudb_create_timer_groups(1);
    data->timer_group->name = taudb_strdup("Harmony Perf");
    taudb_add_timer_group_to_trial(data->taudb_trial, data->timer_group);

    data->param_max = space->len;
    data->client_max = hcfg_int(search_cfg, CFGKEY_CLIENT_COUNT);
    data->taudb_trial->node_count = data->client_max;

    tmpstr = hcfg_get(search_cfg, CFGKEY_TAUDB_STORE_METHOD);
    if (strcmp(tmpstr, "real_time") == 0) {
        data->taudb_store_type = 0;
    }
    else if (strcmp(tmpstr, "one_time") == 0) {
        data->taudb_store_type = 1;
    }
    else {
        search_error("Invalid value for " CFGKEY_TAUDB_STORE_METHOD
                     " configuration key");
        return -1;
    }
    data->total_record_num = hcfg_int(search_cfg, CFGKEY_TAUDB_STORE_NUM);

    data->thread = create_tau_thread(data, data->client_max);

    // Client id map to thread id.
    free(data->client);
    data->client = malloc(data->client_max * sizeof(cinfo_t));
    if (!data->client) {
        search_error("Could not allocate client list");
        return -1;
    }
    memset(data->client, 0, data->client_max * sizeof(cinfo_t));

    data->space = space;
    return 0;
}

int TAUdb_join(hplugin_data_t* data, const char* client)
{
    int idx;
    char node_name[REG_STR_LEN];
    char sys_name[REG_STR_LEN];
    char release[REG_STR_LEN];
    char machine[REG_STR_LEN];
    char proc_num[REG_STR_LEN];
    char cpu_vendor[REG_STR_LEN];
    char cpu_model[REG_STR_LEN];
    char cpu_freq[REG_STR_LEN];
    char cache_size[REG_STR_LEN];

    idx = client_idx(data);
    if (idx < 0)
        return -1;

    sscanf(client, "%[^$$]$$%[^$$]$$%[^$$]$$%[^$$]$$"
           "%[0-9]$$%[^$$]$$%[^$$]$$%[^$$]$$%[^$$]\n",
           node_name, sys_name, release, machine,
           proc_num, cpu_vendor, cpu_model, cpu_freq, cache_size);

    get_tau_metadata(data, &data->thread[idx], sys_name, machine, release,
                     node_name, proc_num, cpu_vendor, cpu_model, cpu_freq,
                     cache_size);
    return 0;
}

/* Invoked upon client reports performance
 * This routine stores the reported performance to TAUdb
 */
int TAUdb_analyze(hplugin_data_t* data, hflow_t* flow, htrial_t* ah_trial)
{
    int idx;
    char timer_name[REG_STR_LEN];

    //taudb_save_metrics(connection, taudb_trial, 1);

    // Create a timer.
    TAUDB_TIMER* timer = taudb_create_timers(1);
    TAUDB_TIMER_VALUE* timer_value = taudb_create_timer_values(1);
    TAUDB_TIMER_CALLPATH* timer_callpath = taudb_create_timer_callpaths(1);
    TAUDB_TIMER_CALL_DATA* timer_call_data = taudb_create_timer_call_data(1);
    //TAUDB_TIMER_GROUP* timer_group = taudb_create_timer_groups(1);

    idx = client_idx(data);
    if (idx < 0)
        return -1;

    // Parse input string and get param information.
    //
    //timer_group->name = taudb_strdup("Harmony Perf");

    ++data->client[idx].timer;
    snprintf(timer_name, sizeof(timer_name), "Runtime_%d_%d",
             idx, data->client[idx].timer);

    timer->name = taudb_strdup(timer_name);

    taudb_add_timer_to_trial(data->taudb_trial, timer);

    taudb_add_timer_to_timer_group(data->timer_group, timer);

    timer_callpath->timer = timer;
    timer_callpath->parent = NULL;
    taudb_add_timer_callpath_to_trial(data->taudb_trial, timer_callpath);

    timer_call_data->key.timer_callpath = timer_callpath;
    timer_call_data->key.thread = &data->thread[idx];
    timer_call_data->calls = 1;
    timer_call_data->subroutines = 0;
    taudb_add_timer_call_data_to_trial(data->taudb_trial, timer_call_data);

    timer_value->metric = data->metric;
    timer_value->inclusive = hperf_unify(&ah_trial->perf);
    timer_value->exclusive = hperf_unify(&ah_trial->perf);
    timer_value->inclusive_percentage = 100.0;
    timer_value->exclusive_percentage = 100.0;
    timer_value->sum_exclusive_squared = 0.0;
    taudb_add_timer_value_to_timer_call_data(timer_call_data, timer_value);

    if (save_timer_parameter(data, timer, ah_trial) != 0)
        return -1;

    // Save the trial.
    if (data->taudb_store_type == 0) {
        if (data->save_counter < data->total_record_num) {
            ++data->save_counter;
        }
        else {
            taudb_compute_statistics(data->taudb_trial);
            taudb_save_trial(data->connection, data->taudb_trial, 1, 1);
            data->taudb_store_type = -1;
        }
    }
    else if (data->taudb_store_type == 1) {
        if (data->save_counter < data->total_record_num) {
            ++data->save_counter;
        }
        else {
            taudb_save_trial(data->connection, data->taudb_trial, 1, 1);
            data->save_counter = 0;
        }
    }
    return 0;
}

/*
 * Finalize a search task.
 */
int TAUdb_fini(hplugin_data_t* data)
{
    boolean update = 1;
    boolean cascade = 1;
    taudb_compute_statistics(data->taudb_trial);
    taudb_save_trial(data->connection, data->taudb_trial, update, cascade);

    // Disconnect from TAUdb.
    taudb_disconnect(data->connection);

    for (int i = 0; i < data->client_max; ++i)
        free( data->client[i].id );
    free(data->client);

    free(data);
    return 0;
}

/*
 * Internal helper function implementation.
 */

int client_idx(hplugin_data_t* data)
{
    int i;
    const char* curr_id;

    curr_id = hcfg_get(search_cfg, CFGKEY_CURRENT_CLIENT);
    if (!curr_id) {
        search_error("Request detected from invalid client");
        return -1;
    }

    for (i = 0; i < data->client_max && data->client[i].id; ++i) {
        if (strcmp(curr_id, data->client[i].id) == 0)
            return i;
    }

    if (i == data->client_max) {
        search_error("Too few clients estimated by " CFGKEY_CLIENT_COUNT);
        return -1;
    }

    data->client[i].id = stralloc(curr_id);
    if (!data->client[i].id) {
        search_error("Could not allocate client id memory");
        return -1;
    }
    return i;
}

int save_timer_parameter(hplugin_data_t* data, TAUDB_TIMER* timer,
                         htrial_t* trial)
{
    int i;
    char buf[32];
    TAUDB_TIMER_PARAMETER* param;

    param = taudb_create_timer_parameters(data->param_max);
    for (i = 0; i < data->param_max; i++) {
        const hval_t* val = &trial->point.term[i];

        // Save name.
        param[i].name = taudb_strdup(data->space->dim[i].name);

        // Get value.
        switch (val->type) {
        case HVAL_INT:
            snprintf(buf, sizeof(buf), "%ld", val->value.i);
            break;

        case HVAL_REAL:
            snprintf(buf, sizeof(buf), "%f", val->value.r);
            break;

        case HVAL_STR:
            snprintf(buf, sizeof(buf), "%s", val->value.s);
            break;

        default:
            search_error("Invalid point value detected");
            return -1;
        }
        param[i].value = taudb_strdup(buf);
        taudb_add_timer_parameter_to_timer(timer, &param[i]);
    }

    return 0;
}

/*
 * Initialize a trial.
 */
int init_tau_trial(hplugin_data_t* data, const char* appName,
                   const char* trialname)
{
    char startTime[32];
    struct tm* current;
    time_t now;

    // Check if the connection has been established.
    taudb_check_connection(data->connection);

    // Create a new trial.
    data->taudb_trial = taudb_create_trials(1);
    data->taudb_trial->name = taudb_strdup(trialname);

    // Set the data source to other.
    data->taudb_trial->data_source =
        taudb_get_data_source_by_id(
            taudb_query_data_sources(data->connection), 999);

    // Create metadata.
    //num_metadata = count_num_metadata();
    TAUDB_PRIMARY_METADATA* pm = taudb_create_primary_metadata(2);
    pm[0].name = taudb_strdup("Application");
    pm[0].value = taudb_strdup(appName);
    taudb_add_primary_metadata_to_trial(data->taudb_trial, &(pm[0]));

    // Get the start time of the task.
    now = time(0);
    current = localtime(&now);
    snprintf(startTime, 64, "%d%d%d", (int)current->tm_hour,
             (int)current->tm_min, (int)current->tm_sec);
    //pm = taudb_create_primary_metadata(5);
    pm[1].name = taudb_strdup("StartTime");
    pm[1].value = taudb_strdup(startTime);
    taudb_add_primary_metadata_to_trial(data->taudb_trial, &(pm[1]));

    boolean update = 0;
    boolean cascade = 1;
    taudb_save_trial(data->connection, data->taudb_trial, update, cascade);

    return 0;
}

TAUDB_THREAD* create_tau_thread(hplugin_data_t* data, int num)
{
    //int ctr = 0;
    data->thread = taudb_create_threads(num);
    for (int i = 0; i < num; i++) {
        data->thread[i].node_rank = i;
        data->thread[i].thread_rank = 1;
        data->thread[i].context_rank = 1;
        data->thread[i].index = 1;
        taudb_add_thread_to_trial(data->taudb_trial, &data->thread[i]);
    }
    return data->thread;
}

/*
 * Get per client metadata.
 */
void get_tau_metadata(hplugin_data_t* data, TAUDB_THREAD* thr, char* opsys,
                      char* machine, char* release, char* hostname,
                      char* procnum, char* cpuvendor, char* cpumodel,
                      char* clkfreq, char* cachesize)
{
    //TAUDB_SECONDARY_METADATA* cur = taudb_create_secondary_metadata(1);
    TAUDB_SECONDARY_METADATA* sm = taudb_create_secondary_metadata(1);

    // Loading os info.
    fprintf(stderr, "Loading OS information.\n");
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("OS");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(opsys);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "OS name loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("Machine");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(machine);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "Machine name loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("Release");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(release);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "Release name loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("HostName");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(hostname);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "Host name loaded.\n");

    // Loading CPU info.
    fprintf(stderr, "Loading CPU information.\n");
    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("ProcNum");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(procnum);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "Processor num loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("CPUVendor");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(cpuvendor);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "CPU vendor loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("CPUModel");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(cpumodel);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "CPU model loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("ClockFreq");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(clkfreq);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "Clock frequency loaded.\n");

    sm = taudb_create_secondary_metadata(1);
    sm->key.thread = thr;
    sm->key.name = taudb_strdup("CacheSize");
    sm->value = malloc(sizeof(char*));
    sm->value[0] = taudb_strdup(cachesize);
    taudb_add_secondary_metadata_to_trial(data->taudb_trial, sm);
    fprintf(stderr, "Cache size loaded.\n");

    //taudb_add_secondary_metadata_to_secondary_metadata(root, &cur);
    fprintf(stderr, "Saving secondary metadata...\n");
    boolean update = 1;
    boolean cascade = 1;
    //taudb_add_secondary_metadata_to_trial(data->taudb_trial, &(sm[0]));
    //taudb_save_secondary_metadata(data->connection, data->taudb_trial,
    //                              update);
    taudb_save_trial(data->connection, data->taudb_trial, update, cascade);
    fprintf(stderr, "Secondary metadata saving complete.\n");
}
