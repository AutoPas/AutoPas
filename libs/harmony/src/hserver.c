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

#include "hserver.h"
#include "httpsvr.h"
#include "hcfg.h"
#include "hmesg.h"
#include "hperf.h"
#include "hsockutil.h"
#include "hutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <libgen.h>
#include <signal.h>
#include <math.h>
#include <getopt.h> // For getopt_long(). Requires _GNU_SOURCE.

#include <sys/types.h>
#include <sys/select.h>
#include <sys/socket.h>

#include <netinet/in.h>
#include <arpa/inet.h>
#include <netinet/tcp.h>

/*
 * Internal helper function prototypes.
 */
static int  launch_session(void);
static int  verbose(const char* fmt, ...);
static int  parse_opts(int argc, char* argv[]);
static int  vars_init(int argc, char* argv[]);
static int  network_init(void);
static int  handle_new_connection(int fd);
static int  handle_unknown_connection(int fd);
static int  handle_client_socket(int fd);
static int  handle_session_socket(void);
static void update_flags(sinfo_t* sinfo, const char* keyval);
static int  update_state(sinfo_t* sinfo);
static int  append_http_log(sinfo_t* sinfo, const hpoint_t* pt, double perf);
static void close_client(int fd);
static void sigint_handler(int signum);

/*
 * Search-related helper function prototypes.
 */
static int      find_search_by_id(int id);
static int      find_search_by_name(const char* name);
static sinfo_t* open_search(int fd);
static sinfo_t* join_search(const char* name, int fd);
static void     close_search(sinfo_t* sinfo);
static void     fini_search(sinfo_t* sinfo);

/*
 * Integer list internal helper function prototypes.
 */
static int add_value(ilist_t* list, int value);
static int find_value(ilist_t* list, int value);
static int remove_index(ilist_t* list, int slot);
static int remove_value(ilist_t* list, int slot);

/*
 * File-local variables.
 */
static int    listen_port = DEFAULT_PORT;
static int    listen_socket;
static int    session_fd;
static fd_set listen_set;
static int    highest_socket;

static ilist_t unknown_fds;
static ilist_t client_fds;
static ilist_t http_fds;

static char* harmony_dir;
static char* session_bin;
static hmesg_t mesg;

static int verbose_flag;
static int done;

/*
 * Exported variables.
 */
sinfo_t* slist;
int slist_cap;

void usage(const char* prog)
{
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "OPTIONS:\n"
"  -p, --port=PORT   Port to listen to on the local host. (Default: %d)\n"
"  -v, --verbose     Print additional information during operation.\n\n",
            listen_port);
}

int main(int argc, char* argv[])
{
    int i, fd_count, retval;
    fd_set ready_set;

    // Parse user options.
    if (parse_opts(argc, argv) != 0)
        return -1;

    // Initialize global variable state.
    if (vars_init(argc, argv) < 0)
        return -1;

    // Initialize the HTTP user interface.
    if (http_init(harmony_dir) < 0)
        return -1;

    // Initialize socket networking service.
    if (network_init() < 0)
        return -1;

    // Launch the underlying search session.
    if (launch_session() != 0)
        return -1;

    while (!done) {
        ready_set = listen_set;
        fd_count = select(highest_socket + 1, &ready_set, NULL, NULL, NULL);
        if (fd_count == -1) {
            if (errno == EINTR)
                continue;

            perror("Error selecting active sockets");
            break;
        }

        if (fd_count > 0) {
            // Before all else, handle input from session process.
            if (FD_ISSET(session_fd, &ready_set)) {
                if (handle_session_socket() != 0)
                    goto shutdown;
            }

            // Handle new connections.
            if (FD_ISSET(listen_socket, &ready_set)) {
                retval = handle_new_connection(listen_socket);
                if (retval > 0) {
                    FD_SET(retval, &listen_set);
                    if (highest_socket < retval)
                        highest_socket = retval;
                }
            }

            // Handle unknown connections (Unneeded if we switch to UDP).
            for (i = 0; i < unknown_fds.len; ++i) {
                if (FD_ISSET(unknown_fds.slot[i], &ready_set)) {
                    if (handle_unknown_connection(unknown_fds.slot[i]) != 0)
                        FD_CLR(unknown_fds.slot[i], &listen_set);
                    remove_index(&unknown_fds, i--);
                }
            }

            // Handle Harmony messages.
            for (i = 0; i < client_fds.len; ++i) {
                if (FD_ISSET(client_fds.slot[i], &ready_set)) {
                    retval = handle_client_socket(client_fds.slot[i]);
                    if (retval > 0) goto shutdown;
                    if (retval < 0) {
                        close_client(client_fds.slot[i]);
                        remove_index(&client_fds, i--);
                    }
                }
            }

            // Handle http requests.
            for (i = 0; i < http_fds.len; ++i) {
                if (FD_ISSET(http_fds.slot[i], &ready_set)) {
                    if (handle_http_socket(http_fds.slot[i]) != 0) {
                        FD_CLR(http_fds.slot[i], &listen_set);
                        remove_index(&http_fds, i--);
                    }
                }
            }
        }
    }

  shutdown:
    for (i = 0; i < slist_cap; ++i)
        fini_search(&slist[i]);
    free(slist);

    free(unknown_fds.slot);
    free(client_fds.slot);
    free(http_fds.slot);

    free(harmony_dir);
    free(session_bin);
    hmesg_fini(&mesg);

    return 0;
}

int verbose(const char* fmt, ...)
{
    int retval;
    va_list ap;

    if (!verbose_flag)
        return 0;

    va_start(ap, fmt);
    retval = vfprintf(stderr, fmt, ap);
    va_end(ap);

    return retval;
}

int parse_opts(int argc, char* argv[])
{
    int c;
    static struct option long_options[] = {
        {"port",    required_argument, NULL, 'p'},
        {"verbose", no_argument,       NULL, 'v'},
        {NULL, 0, NULL, 0}
    };

    while (1) {
        c = getopt_long(argc, argv, "p:v", long_options, NULL);
        if (c == -1)
            break;

        switch(c) {
        case 'p': listen_port = atoi(optarg); break;
        case 'v': verbose_flag = 1; break;

        case ':':
            usage(argv[0]);
            fprintf(stderr, "\nOption ('%c') requires an argument.\n", optopt);
            break;

        case '?':
        default:
            usage(argv[0]);
            fprintf(stderr, "\nInvalid argument ('%c').\n", optopt);
            return -1;
        }
    }

    return 0;
}

int vars_init(int argc, char* argv[])
{
    char* tmppath;
    char* binfile;

    // Ignore signal for writes to broken pipes/sockets.
    if (signal(SIGPIPE, SIG_IGN) == SIG_ERR) {
        perror("Error ignoring SIGPIPE");
        return -1;
    }

    // Use SIGINT as a method for graceful server exit.
    if (signal(SIGINT, sigint_handler) == SIG_ERR) {
        perror("Error installing handler for SIGINT");
        return -1;
    }

    //Determine directory where this binary is located
    tmppath = stralloc(argv[0]);
    binfile = stralloc(basename(tmppath));
    free(tmppath);

    if ( (tmppath = getenv(CFGKEY_HARMONY_HOME))) {
        harmony_dir = stralloc(tmppath);
        verbose(CFGKEY_HARMONY_HOME " is %s\n", harmony_dir);
    }
    else {
        if (strchr(argv[0], '/'))
            tmppath = stralloc(argv[0]);
        else
            tmppath = stralloc(search_path(binfile));

        harmony_dir = dirname(tmppath);
        if (strcmp(harmony_dir, ".") == 0)
            harmony_dir = stralloc("..");
        else
            harmony_dir = stralloc(dirname(harmony_dir));
        free(tmppath);

        verbose("Detected %s/ as HARMONY_HOME\n", harmony_dir);
    }
    free(binfile);

    // Find supporting binaries and shared objects.
    session_bin = sprintf_alloc("%s/libexec/" SESSION_CORE_EXECFILE,
                                harmony_dir);
    if (!file_exists(session_bin)) {
        fprintf(stderr, "Could not find support files in "
                CFGKEY_HARMONY_HOME "\n");
        return -1;
    }

    // Prepare the search information lists.
    if (array_grow(&slist, &slist_cap, sizeof(*slist)) != 0) {
        perror("Could not allocate memory for session file descriptor list");
        return -1;
    }
    for (int i = 0; i < slist_cap; ++i)
        slist[i].id = -1; // Initialize slist with non-zero defaults.

    return 0;
}

int network_init(void)
{
    int optval;
    struct sockaddr_in addr;

    // Create a listening socket.
    verbose("Listening on TCP port: %d\n", listen_port);
    listen_socket = socket(AF_INET, SOCK_STREAM, 0);
    if (listen_socket < 0) {
        perror("Could not create listening socket");
        return -1;
    }

    // Set socket options.
    //
    // Try to make the port reusable and have it close as fast as
    // possible without waiting in unnecessary wait states on close.
    //
    optval = 1;
    if (setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR,
                   &optval, sizeof(optval)) < 0)
    {
        perror("Could not set options on listening socket");
        return -1;
    }

    // Initialize the socket address.
    addr.sin_family      = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port        = htons((unsigned short)listen_port);

    // Bind the socket to the desired port.
    if (bind(listen_socket, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("Could not bind socket to listening address");
        return -1;
    }

    // Set the file descriptor set.
    FD_ZERO(&listen_set);
    FD_SET(listen_socket, &listen_set);
    highest_socket = listen_socket;

    if (listen(listen_socket, SOMAXCONN) < 0) {
        perror("Could not listen on listening socket");
        return -1;
    }

    return 0;
}

int handle_new_connection(int fd)
{
    struct sockaddr_in addr;
    socklen_t addrlen = sizeof(addr);
    int newfd;

    newfd = accept(fd, (struct sockaddr*)&addr, &addrlen);
    if (newfd < 0) {
        perror("Error accepting new connection");
        return -1;
    }

    // The peer may not have sent data yet.  To prevent blocking, lets
    // stash the new socket in the unknown_list until we know it has data.
    //
    if (add_value(&unknown_fds, newfd) < 0)
        return -1;

    verbose("Accepted connection from %s as socket %d\n",
            inet_ntoa(addr.sin_addr), newfd);
    return newfd;
}

int handle_unknown_connection(int fd)
{
    unsigned int header;

    int readlen = recv(fd, &header, sizeof(header), MSG_PEEK | MSG_DONTWAIT);
    if (readlen < 0 || (unsigned int)readlen < sizeof(header)) {
        // Error on recv, or insufficient data.  Close the connection.
        printf("Can't determine type for socket %d.  Closing.\n", fd);
        if (close(fd) < 0)
            perror("Error closing connection");
        return -1;
    }
    else if (ntohl(header) == HMESG_MAGIC) {
        // This is a communication socket from a Harmony client.
        if (add_value(&client_fds, fd) < 0) {
            if (close(fd) != 0)
                perror("Error closing Harmony connection");
        }
    }
    else {
        // Consider this an HTTP communication socket.
        if (http_fds.len >= http_connection_limit) {
            printf("Hit HTTP connection limit on socket %d.  Closing.\n", fd);
            http_send_error(fd, 503, NULL);
            if (close(fd) != 0)
                perror("Error closing HTTP connection");
            return -1;
        }

        if (add_value(&http_fds, fd) < 0) {
            if (close(fd) != 0)
                perror("Error closing HTTP connection");
        }
    }
    return 0;
}

int handle_client_socket(int fd)
{
    int idx;
    double perf;
    sinfo_t* sinfo = NULL;

    int retval = mesg_recv(fd, &mesg);
    if (retval <  0) goto error;
    if (retval == 0) return -1;

    // Sanity check input.
    if (mesg.type != HMESG_SESSION && mesg.type != HMESG_JOIN) {
        if (mesg.dest < 0 || mesg.dest >= slist_cap ||
            slist[ mesg.dest ].id == -1)
        {
            mesg.data.string = "Invalid message destination";
            goto error;
        }
        sinfo = &slist[ mesg.dest ];
    }

    switch (mesg.type) {
    case HMESG_SESSION:
        sinfo = open_search(fd);
        if (!sinfo)
            goto error;
        break;

    case HMESG_JOIN:
        sinfo = join_search(mesg.data.string, fd);
        if (!sinfo)
            goto error;
        break;

    case HMESG_SETCFG:
        update_flags(sinfo, mesg.data.string);
        break;

    case HMESG_REPORT:
        perf = hperf_unify(mesg.data.perf);

        for (idx = 0; idx < sinfo->fetched_len; ++idx) {
            if (sinfo->fetched[idx].id == mesg.data.point->id)
                break;
        }
        if (idx < sinfo->fetched_len) {
            const hpoint_t* pt = &sinfo->fetched[idx];

            // Copy point from fetched list to HTTP log.
            if (append_http_log(sinfo, pt, perf) != 0) {
                mesg.data.string = "Could not append to HTTP log";
                goto error;
            }

            // Remove point from fetched list.
            --sinfo->fetched_len;
            if (idx < sinfo->fetched_len) {
                if (hpoint_copy(&sinfo->fetched[idx],
                                &sinfo->fetched[sinfo->fetched_len]) != 0)
                {
                    mesg.data.string = "Could not remove fetch list point";
                    goto error;
                }
            }
        }
        else {
            // Copy point from fetched list to HTTP log.
            if (append_http_log(sinfo, &sinfo->best, perf) != 0) {
                mesg.data.string = "Could not copy best point to HTTP log";
                goto error;
            }
        }
        break;

    case HMESG_GETCFG:
    case HMESG_BEST:
    case HMESG_FETCH:
    case HMESG_COMMAND:
        break;

    default:
        goto error;
    }

    // Overwrite the message source and destination.
    mesg.src = fd;
    if (mesg.type == HMESG_SESSION || mesg.type == HMESG_JOIN)
        mesg.dest = -1;
    else
        mesg.dest = sinfo->id;

    retval = mesg_forward(session_fd, &mesg);
    if (retval == 0) goto shutdown;
    if (retval <  0) {
        mesg.data.string = "Could not forward message to session";
        goto error;
    }

    // Record that the client has sent a request.
    if (add_value(&sinfo->request, fd) < 0) {
        mesg.data.string = "Server error: Could not add to request list";
        goto error;
    }

    return 0;

  error:
    mesg.dest   = mesg.src;
    mesg.src    = -1;
    mesg.status = HMESG_STATUS_FAIL;
    mesg_send(fd, &mesg);
    return 0; // Do not close client socket on failure.

  shutdown:
    fprintf(stderr, "Session socket closed. Shutting server down.\n");
    return 1;
}

int handle_session_socket(void)
{
    int close_flag = 0;

    int retval = mesg_recv(session_fd, &mesg);
    if (retval < 1) {
        if (retval == 0) fprintf(stderr, "Session socket closed.");
        if (retval <  0) fprintf(stderr, "Malformed message from session.");
        fprintf(stderr, " Shutting server down.\n");

        return -1;
    }

    int slist_idx;
    if (mesg.type != HMESG_SESSION)
        slist_idx = find_search_by_id(mesg.src); // Real search ID.
    else
        slist_idx = find_search_by_id(-mesg.dest); // Temporary bootstrap ID.

    if (slist_idx == -1) {
        fprintf(stderr, "No info found for session message.  Ignoring.\n");
        return 0;
    }
    sinfo_t* sinfo = &slist[slist_idx];

    switch (mesg.type) {
    case HMESG_SESSION:
        if (mesg.status == HMESG_STATUS_OK) {
            // Update search information with real ID
            sinfo->id = mesg.src;

            // Add the client originating client as a search participant.
            if (add_value(&sinfo->client, mesg.dest) < 0) {
                mesg.data.string = "Server error: Could not grow client list";
                goto error;
            }
        }
        else {
            // Search was not established.  Close it.
            close_flag = 1;
        }
        break;

    case HMESG_JOIN:
        if (mesg.status == HMESG_STATUS_OK) {
            if (add_value(&sinfo->client, mesg.dest) < 0) {
                mesg.data.string = "Server error: Could not grow client list";
                goto error;
            }
        }
        break;

    case HMESG_COMMAND:
        if (mesg.status == HMESG_STATUS_OK) {
            if (strcmp(mesg.data.string, "leave") == 0)
                remove_value(&sinfo->client, mesg.dest);

            if (strcmp(mesg.data.string, "kill") == 0)
                close_flag = 1;
        }
        break;

    case HMESG_GETCFG:
        update_flags(sinfo, mesg.data.string);
        break;

    case HMESG_FETCH:
        if (mesg.status == HMESG_STATUS_OK) {
            // Log this point before we forward it to the client.
            if (sinfo->fetched_len == sinfo->fetched_cap) {
                if (array_grow(&sinfo->fetched, &sinfo->fetched_cap,
                               sizeof(*sinfo->fetched)) != 0)
                {
                    mesg.data.string = "Server error: Couldn't grow fetch log";
                    goto error;
                }
            }

            if (hpoint_copy(&sinfo->fetched[sinfo->fetched_len],
                            mesg.data.point) != 0)
            {
                mesg.data.string = "Server error: Couldn't add to HTTP log";
                goto error;
            }

            if (hpoint_align(&sinfo->fetched[sinfo->fetched_len],
                             &sinfo->space) != 0)
            {
                mesg.data.string = "Server error: Couldn't align point";
                goto error;
            }
            ++sinfo->fetched_len;
        }
        break;

    case HMESG_REPORT:
        if (mesg.status == HMESG_STATUS_OK)
            ++sinfo->reported;
        break;

    case HMESG_SETCFG:
    case HMESG_BEST:
        break;

    default:
        mesg.data.string = "Server error: Invalid message type from session";
        goto error;
    }

    if (sinfo && update_state(sinfo) != 0) {
        mesg.data.string = "Server error: Could not update search state";
        goto error;
    }

    // Messages with a negative destination field came within hserver itself.
    // No need to forward those messages.
    //
    if (mesg.dest >= 0) {
        // Overwrite the message source.
        mesg.src = slist_idx;

        if (mesg_forward(mesg.dest, &mesg) < 1) {
            perror("Error forwarding message to client");
            close_client(mesg.dest);
        }

        // Record that this client has seen a response to its request.
        remove_value(&sinfo->request, mesg.dest);
    }

    if (close_flag)
        close_search(sinfo);

    return 0;

  error:
    fprintf(stderr, "%s. Shutting server down.\n", mesg.data.string);
    return -1;
}

int launch_session(void)
{
    // Fork and exec a session handler.
    char* const child_argv[] = {session_bin,
                                harmony_dir,
                                NULL};
    session_fd = socket_launch(session_bin, child_argv, NULL);
    if (session_fd < 0) {
        perror("Could not launch session process");
        return -1;
    }

    FD_SET(session_fd, &listen_set);
    if (highest_socket < session_fd)
        highest_socket = session_fd;

    return 0;
}

void update_flags(sinfo_t* sinfo, const char* keyval)
{
    int val = 0;

    sscanf(keyval, CFGKEY_CONVERGED "=%n", &val);
    if (val) {
        if (hcfg_parse_bool(&keyval[val]))
            sinfo->flags |= FLAG_CONVERGED;
        else
            sinfo->flags &= ~FLAG_CONVERGED;
        return;
    }

    sscanf(keyval, CFGKEY_PAUSED "=%n", &val);
    if (val) {
        if (hcfg_parse_bool(&keyval[val]))
            sinfo->flags |= FLAG_PAUSED;
        else
            sinfo->flags &= ~FLAG_PAUSED;
        return;
    }
}

int update_state(sinfo_t* sinfo)
{
    if (mesg.type == HMESG_SESSION ||
        mesg.type == HMESG_JOIN)
    {
        // Session state doesn't need to be updated yet.
        return 0;
    }

    if (mesg.state.space->id > sinfo->space.id) {
        if (hspace_copy(&sinfo->space, mesg.state.space) != 0) {
            perror("Could not copy search space");
            return -1;
        }
    }

    if (mesg.state.best->id > sinfo->best.id) {
        int i;
        for (i = sinfo->log_len - 1; i >= 0; --i) {
            if (sinfo->log[i].pt.id == mesg.state.best->id) {
                sinfo->best_perf = sinfo->log[i].perf;
                break;
            }
        }
        if (i < 0)
            sinfo->best_perf = NAN;

        if (hpoint_copy(&sinfo->best, mesg.state.best) != 0) {
            perror("Internal error copying hpoint to best");
            return -1;
        }

        if (hpoint_align(&sinfo->best, &sinfo->space) != 0) {
            perror("Could not align best point to search space");
            return -1;
        }
    }
    return 0;
}


int append_http_log(sinfo_t* sinfo, const hpoint_t* pt, double perf)
{
    http_log_t* entry;

    // Extend HTTP log if necessary.
    if (sinfo->log_len == sinfo->log_cap) {
        if (array_grow(&sinfo->log, &sinfo->log_cap,
                       sizeof(*sinfo->log)) != 0)
        {
            perror("Could not grow HTTP log");
            return -1;
        }
    }
    entry = &sinfo->log[sinfo->log_len];

    if (hpoint_copy(&entry->pt, pt) != 0) {
        perror("Internal error copying point into HTTP log");
        return -1;
    }

    entry->perf = perf;
    if (gettimeofday(&entry->stamp, NULL) != 0)
        return -1;

    ++sinfo->log_len;
    return 0;
}

int request_command(sinfo_t* sinfo, const char* command)
{
    int retval = 0;

    mesg.dest = sinfo->id;
    mesg.src = -1;
    mesg.type = HMESG_COMMAND;
    mesg.status = HMESG_STATUS_REQ;
    mesg.state.space = &sinfo->space;
    mesg.state.best = &sinfo->best;
    mesg.state.client = "<hserver>";
    mesg.data.string = command;

    if (mesg_send(session_fd, &mesg) < 1)
        retval = -1;

    return retval;
}

int request_refresh(sinfo_t* sinfo)
{
    mesg.dest = sinfo->id;
    mesg.src = -1;
    mesg.type = HMESG_GETCFG;
    mesg.status = HMESG_STATUS_REQ;
    mesg.state.space = &sinfo->space;
    mesg.state.best = &sinfo->best;
    mesg.state.client = "<hserver>";

    mesg.data.string = CFGKEY_CONVERGED;
    if (mesg_send(session_fd, &mesg) < 1)
        return -1;

    mesg.data.string = CFGKEY_PAUSED;
    if (mesg_send(session_fd, &mesg) < 1)
        return -1;

    return 0;
}

int request_setcfg(sinfo_t* sinfo, const char* key, const char* val)
{
    char* buf = sprintf_alloc("%s=%s", key, val ? val : "");
    int retval = 0;

    mesg.dest = sinfo->id;
    mesg.src = -1;
    mesg.type = HMESG_SETCFG;
    mesg.status = HMESG_STATUS_REQ;
    mesg.data.string = buf;
    mesg.state.space = &sinfo->space;
    mesg.state.best = &sinfo->best;
    mesg.state.client = "<hserver>";

    if (mesg_send(session_fd, &mesg) < 1)
        retval = -1;

    free(buf);
    return retval;
}

/*
 * Handles an unexpected client departure.  The session must be
 * informed if the client was actively participating in any searches.
 */
void close_client(int fd)
{
    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i].id == -1)
            continue;

        if (remove_value(&slist[i].client, fd) == 0)
            request_command(&slist[i], "leave");
    }

    if (shutdown(fd, SHUT_RDWR) != 0 || close(fd) != 0)
        perror("Error closing client connection post error");

    FD_CLR(fd, &listen_set);
}

/*
 * Search-related helper function implementation.
 */

int find_search_by_id(int id)
{
    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i].id == id)
            return i;
    }
    return -1;
}

int find_search_by_name(const char* name)
{
    int idx = -1;

    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i].id == -1) {
            if (idx < 0)
                idx = i; // Mark the first available index as we pass it.
        }
        else if (strcmp(slist[i].space.name, name) == 0) {
            return i;
        }
    }

    // Extend slist, if necessary.
    if (idx < 0) {
        idx = slist_cap;
        if (array_grow(&slist, &slist_cap, sizeof(*slist)) != 0) {
            mesg.data.string = "Server error: Could not extend search list";
            return -1;
        }

        // Initialize any sinfo_t objects that were created.
        for (int i = idx; i < slist_cap; ++i)
            slist[i].id = -1;
    }
    return idx;
}

sinfo_t* open_search(int fd)
{
    int idx = find_search_by_name(mesg.state.space->name);
    if (idx < 0)
        return NULL;

    sinfo_t* sinfo = &slist[idx];
    if (sinfo->id != -1) {
        mesg.data.string = "Search name already exists";
        return NULL;
    }

    // Initialize the available sinfo slot.
    if (hspace_copy(&sinfo->space, mesg.state.space) != 0) {
        mesg.data.string = "Server error: Could not copy initial space info";
        return NULL;
    }


    // Initialize HTTP server fields.
    if (gettimeofday(&sinfo->start, NULL) != 0) {
        mesg.data.string = "Server error: Could not set search start time";
        return NULL;
    }

    const char* cfgstr = hcfg_get(mesg.data.cfg, CFGKEY_STRATEGY);
    if (!cfgstr) {
        int clients = hcfg_int(mesg.data.cfg, CFGKEY_CLIENT_COUNT);
        if (clients < 1)
            cfgstr = "???";
        else if (clients == 1)
            cfgstr = "nm.so";
        else
            cfgstr = "pro.so";
    }
    free(sinfo->strategy);
    sinfo->strategy = stralloc(cfgstr);

    sinfo->client.len = 0;
    sinfo->request.len = 0;
    sinfo->best.id = 0;
    sinfo->best_perf = HUGE_VAL;
    sinfo->flags = 0x0;
    sinfo->log_len = 0;
    sinfo->fetched_len = 0;
    sinfo->reported = 0;

    // The true ID of the search won't be known until the session
    // responds positively to the HMESG_SESSION request.  Until that
    // time, base the temporary ID on the requesting client's socket
    // number.  To avoid conflicting with real IDs, the socket number
    // is negated.
    //
    sinfo->id = -fd;

    return sinfo;
}

/*
 * Sanity checks when a client requests to join a search task.
 */
sinfo_t* join_search(const char* name, int fd)
{
    int s_idx = find_search_by_name(name);
    if (s_idx < 0)
        return NULL;

    sinfo_t* sinfo = &slist[s_idx];
    if (sinfo->id == -1) {
        mesg.data.string = "No search currently exists with this name";
        return NULL;
    }

    // A single client may not join a search more than once concurrently.
    if (find_value(&sinfo->client, fd) < sinfo->client.len) {
        mesg.data.string = "Client already participating in this search";
        return NULL;
    }
    return sinfo;
}

/*
 * Update server bookkeeping after the search session reports it has
 * successfully killed the search task.
 */
void close_search(sinfo_t* sinfo)
{
    sinfo->id = -1;

    // Prepare a "dead search task" message.
    mesg.dest = -1;
    mesg.src = -1;
    mesg.type = HMESG_UNKNOWN;
    mesg.status = HMESG_STATUS_FAIL;

    // Inform clients awaiting a response from the dead search.
    for (int i = 0; i < sinfo->request.len; ++i)
        mesg_send(sinfo->request.slot[i], &mesg);
}

/*
 * Update server bookkeeping after the search session reports it has
 * successfully killed the search task.
 */
void fini_search(sinfo_t* sinfo)
{
    for (int i = 0; i < sinfo->fetched_cap; ++i)
        hpoint_fini(&sinfo->fetched[i]);
    free(sinfo->fetched);

    for (int i = 0; i < sinfo->log_cap; ++i)
        hpoint_fini(&sinfo->log[i].pt);
    free(sinfo->log);

    free(sinfo->strategy);
    hspace_fini(&sinfo->space);

    free(sinfo->client.slot);
    free(sinfo->request.slot);
    hpoint_fini(&sinfo->best);
}

/*
 * Integer list internal helper function implementation.
 */

int add_value(ilist_t* list, int value)
{
    int idx = find_value(list, value);
    if (idx < list->len)
        return 0;

    if (idx == list->cap) {
        if (array_grow(&list->slot, &list->cap, sizeof(*list->slot)) != 0)
            return -1;
    }

    list->slot[idx] = value;
    ++list->len;
    return 1;
}

int find_value(ilist_t* list, int value)
{
    int idx;
    for (idx = 0; idx < list->len; ++idx) {
        if (list->slot[idx] == value)
            break;
    }
    return idx;
}

int remove_index(ilist_t* list, int idx)
{
    if (idx < 0 || idx >= list->len)
        return -1;

    if (idx < --list->len)
        list->slot[idx] = list->slot[ list->len ];
    return 0;
}

int remove_value(ilist_t* list, int value)
{
    int idx = find_value(list, value);
    return remove_index(list, idx);
}

void sigint_handler(int signum)
{
    fprintf(stderr, "\nCaught signal %d. Shutting down the server.\n", signum);
    done = 1;
}
