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
 * \page codesvr Code Server (codegen.so)
 *
 * This processing layer passes messages from the tuning session to a
 * running code generation server.  Details on how to configure and
 * run a code generation tuning session are provided in
 * `code-server/README` of the distribution tarball.
 *
 * ## Code Server URLs ##
 * The code server uses a proprietary set of URL's to determine the
 * destination for various communications.  They take the following
 * form:
 *
 *     dir://<path>
 *     ssh://[user@]<host>[:port]/<path>
 *     tcp://<host>[:port]
 *
 * All paths are considered relative.  Use an additional slash to
 * specify an absolute path.  For example:
 *
 *     dir:///tmp
 *     ssh://code.server.net:2222//tmp/codegen
 */

#include "hlayer.h"
#include "session-core.h"
#include "hspace.h"
#include "hcfg.h"
#include "hmesg.h"
#include "hpoint.h"
#include "hutil.h"
#include "hsockutil.h"

#include <stdlib.h>
#include <string.h>
#include <errno.h>

/*
 * Name used to identify this plugin layer.
 * All Harmony plugin layers must define this variable.
 */
const char hplugin_name[] = "codegen";

/*
 * Configuration variables used in this plugin.
 * These will automatically be registered by session-core upon load.
 */
const hcfg_info_t hplugin_keyinfo[] = {
    { CFGKEY_SERVER_URL, NULL,
      "Destination of messages from session to code server." },
    { CFGKEY_TARGET_URL, NULL,
      "Destination of binary files from code server to target application." },
    { CFGKEY_REPLY_URL, NULL,
      "Destination of reply messages from code server to session. "
      "A reasonable directory in " CFGKEY_TMPDIR " will be generated "
      "by default if left blank." },
    { CFGKEY_TMPDIR, "/tmp",
      "Session local path to store temporary files." },
    { CFGKEY_SLAVE_LIST, NULL,
      "Comma separated list of 'host n' pairs, where 'n' slaves will run "
      "on 'host'.  For example 'bunker.cs.umd.edu 4, nexcor.cs.umd.edu 4' "
      "will instruct the code server to spawn 8 total slaves between two "
      "machines." },
    { CFGKEY_SLAVE_PATH, NULL,
      "Path on slave host to store generated binaries before being sent "
      "to the target application." },
    { NULL }
};

typedef enum {
    CODEGEN_STATUS_UNKNOWN = 0x0,
    CODEGEN_STATUS_REQUESTED,
    CODEGEN_STATUS_COMPLETE,

    CODEGEN_STATUS_MAX
} codegen_status;

typedef struct codegen_log {
    codegen_status status;
    hpoint_t point;
} codegen_log_t;

/*
 * Structure to hold all data needed by an individual search instance.
 *
 * To support multiple parallel search instances, no global variables
 * should be defined or used in this plug-in layer.  They should
 * instead be defined as a part of this structure.
 */
struct hplugin_data {
    codegen_log_t* cglog;
    int cglog_len, cglog_cap;

    hmesg_t mesg;
    int sockfd;

    char* buf;
    int buflen;
};

/*
 * Long-running code generation callback prototype.
 */
static int codegen_callback(int fd, void* data,
                            hflow_t* flow, int n, htrial_t** trial);

/*
 * Internal helper function prototypes.
 */
static int cglog_find(hplugin_data_t* data, const hpoint_t* pt);
static int cglog_insert(hplugin_data_t* data, const hpoint_t* pt);
static int url_connect(hplugin_data_t* data, const char* url);

/*
 * Allocate memory for a new search task.
 */
hplugin_data_t* codegen_alloc(void)
{
    hplugin_data_t* retval = calloc(1, sizeof(*retval));
    if (!retval)
        return NULL;

    return retval;
}

/*
 * Initialize (or re-initialize) data for this search instance.
 */
int codegen_init(hplugin_data_t* data, hspace_t* space)
{
    const char* url;

    url = hcfg_get(search_cfg, CFGKEY_TARGET_URL);
    if (!url || url[0] == '\0') {
        search_error("Destination URL for"
                     " generated code objects not specified");
        return -1;
    }

    url = hcfg_get(search_cfg, CFGKEY_SERVER_URL);
    if (!url || url[0] == '\0') {
        search_error("Codegen server URL not specified");
        return -1;
    }

    data->sockfd = url_connect(data, url);
    if (data->sockfd == -1) {
        search_error("Invalid codegen server URL");
        return -1;
    }

    /* In the future, we should rewrite the code generator system to
     * avoid using hmesg_t types.  Until then, generate a fake
     * HMESG_SESSION message to maintain compatibility.
     */
    data->mesg.type = HMESG_SESSION;
    data->mesg.status = HMESG_STATUS_REQ;
    data->mesg.state.space = space;
    data->mesg.data.cfg = search_cfg;

    if (mesg_send(data->sockfd, &data->mesg) < 1)
        return -1;

    if (mesg_recv(data->sockfd, &data->mesg) < 1)
        return -1;

    // TODO: Need a way to unregister a callback for re-initialization.
    if (search_callback_generate(data->sockfd, data, codegen_callback) != 0) {
        search_error("Could not register callback for codegen plugin");
        return -1;
    }

    return 0;
}

/*
 * Called after search driver generates a candidate point, but before
 * that point is returned to the client API.
 *
 * This routine should fill the flow variable appropriately and return
 * 0 upon success.  Otherwise, it should call search_error() with a
 * human-readable error message and return -1.
 */
int codegen_generate(hplugin_data_t* data, hflow_t* flow, htrial_t* trial)
{
    int i;

    i = cglog_find(data, &trial->point);
    if (i >= 0) {
        if (data->cglog[i].status == CODEGEN_STATUS_COMPLETE)
            flow->status = HFLOW_ACCEPT;

        if (data->cglog[i].status == CODEGEN_STATUS_REQUESTED)
            flow->status = HFLOW_WAIT;

        return 0;
    }

    if (cglog_insert(data, &trial->point) != 0) {
        search_error("Could not grow codegen log");
        return -1;
    }

    data->mesg.type = HMESG_FETCH;
    data->mesg.status = HMESG_STATUS_OK;
    data->mesg.data.point = &trial->point;

    if (mesg_send(data->sockfd, &data->mesg) < 1) {
        search_error( strerror(errno) );
        return -1;
    }

    flow->status = HFLOW_WAIT;
    return 0;
}

/*
 * Long-running code generation callback implementation.
 */
int codegen_callback(int fd, void* data_ptr,
                     hflow_t* flow, int n, htrial_t** trial)
{
    hplugin_data_t* data = (hplugin_data_t*) data_ptr;
    int i;

    if (mesg_recv(fd, &data->mesg) < 1)
        return -1;

    i = cglog_find(data, data->mesg.data.point);
    if (i < 0) {
        search_error("Could not find point from code server in log");
        return -1;
    }
    data->cglog[i].status = CODEGEN_STATUS_COMPLETE;

    // Search waitlist for index of returned point.
    for (i = 0; i < n; ++i) {
        if (trial[i]->point.id == data->mesg.data.point->id) {
            flow->status = HFLOW_ACCEPT;
            return i;
        }
    }

    search_error("Could not find point from code server in waitlist");
    return -1;
}

/*
 * Free memory associated with this search task.
 */
int codegen_fini(hplugin_data_t* data)
{
    close(data->sockfd);
    free(data->buf);
    free(data);
    data = NULL;

    return 0;
}

/*
 * Internal helper function implementation.
 */

int cglog_find(hplugin_data_t* data, const hpoint_t* point)
{
    for (int i = 0; i < data->cglog_len; ++i) {
        if (hpoint_eq(point, &data->cglog[i].point))
            return i;
    }
    return -1;
}

int cglog_insert(hplugin_data_t* data, const hpoint_t* point)
{
    if (data->cglog_len == data->cglog_cap) {
        if (array_grow(&data->cglog, &data->cglog_cap,
                       sizeof(*data->cglog)) != 0)
            return -1;
    }

    hpoint_copy(&data->cglog[data->cglog_len].point, point);
    data->cglog[data->cglog_len].status = CODEGEN_STATUS_REQUESTED;
    ++data->cglog_len;

    return 0;
}

int url_connect(hplugin_data_t* data, const char* url)
{
    const char* ptr;
    char* helper_argv[2];
    int count;

    ptr = strstr(url, "//");
    if (!ptr)
        return -1;

    if (strncmp("dir:", url, ptr - url) == 0 ||
        strncmp("ssh:", url, ptr - url) == 0)
    {
        ptr = hcfg_get(search_cfg, CFGKEY_HARMONY_HOME);
        if (!ptr)
            return -1;

        count = snprintf_grow(&data->buf, &data->buflen,
                              "%s/libexec/codegen-helper", ptr);
        if (count < 0)
            return -1;

        helper_argv[0] = data->buf;
        helper_argv[1] = NULL;
        return socket_launch(data->buf, helper_argv, NULL);
    }
    else if (strncmp("tcp:", url, ptr - url) == 0) {
        // Not implemented yet.
    }
    return -1;
}
