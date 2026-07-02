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
#define _XOPEN_SOURCE 600 // Needed for S_ISSOCK, sigaction(), and
                          // gethostname().

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <netdb.h>
#include <time.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/select.h>
#include <sys/socket.h>

#include "hmesg.h"
#include "hspace.h"
#include "hcfg.h"
#include "hutil.h"
#include "hsockutil.h"

#define POLL_TIME 250000

/*
 * Internal helper function prototypes.
 */
void graceful_exit(int signum);
int  init_signals(void);
int  init_comm(void);
int  url_parse(const char* url);
int  reply_url_build(void);
int  mesg_write(int id);
int  mesg_read(int id);
int  read_loop(int fd, char* readbuf, int readlen);
int  write_loop(int fd, char* writebuf, int writelen);

/*
 * Global variables.
 */
hmesg_t mesg = HMESG_INITIALIZER;
char*   reply_dir;
int     reply_dir_created, signal_caught;
char*   scp_cmd;
char*   scp_dst;
char*   tmp_dir;
char*   sess_name;

char* buf;
int   buflen;

int main(int argc, char* argv[])
{
    struct stat sb;
    struct dirent* dent;
    DIR* dirfd;
    int next_id, curr_id, count, retval;
    fd_set fds, ready_fds;
    struct timeval polltime;

    // Check that we have been launched correctly by checking that
    // STDIN_FILENO is a socket descriptor.
    //
    // Print an error message to stderr if this is not the case.  This
    // should be the only time an error message is printed to stdout
    // or stderr.
    //
    if (fstat(STDIN_FILENO, &sb) != 0) {
        perror("Could not determine the status of STDIN");
        return -1;
    }

    if (!S_ISSOCK(sb.st_mode)) {
        fprintf(stderr, "%s should not be launched manually.\n", argv[0]);
        return -1;
    }

    retval = 0;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);

    buflen = 4096;
    buf = malloc(buflen);
    if (!buf) {
        mesg.data.string = "Could not allocate temporary memory buffer";
        goto error;
    }

    if (init_signals() < 0)
        goto error;

    if (init_comm() < 0)
        goto error;

    curr_id = -1;
    next_id = 0;

    while (!signal_caught) {
        ready_fds = fds;
        polltime.tv_sec = 0;
        polltime.tv_usec = POLL_TIME;
        count = select(STDIN_FILENO + 1, &ready_fds, NULL, NULL, &polltime);
        if (count == -1) {
            if (errno == EINTR)
                continue;

            mesg.data.string = strerror(errno);
            goto error;
        }

        // Forward incoming message from Harmony server to code server.
        if (FD_ISSET(STDIN_FILENO, &ready_fds)) {
            if (mesg_recv(STDIN_FILENO, &mesg) < 1)
                goto error;

            if (mesg_write(next_id) < 0)
                goto error;

            ++next_id;
        }

        // Forward incoming message from code server to Harmony server.
        do {
            count = mesg_read(curr_id);
            if (count < 0) goto error;
            if (count > 0) ++curr_id;
        } while (count > 0);
    }
    goto cleanup;

  error:
    mesg.status = HMESG_STATUS_FAIL;
    mesg_send(STDIN_FILENO, &mesg);
    retval = -1;

  cleanup:
    dirfd = opendir(reply_dir);
    while ( (dent = readdir(dirfd))) {
        if (strcmp(dent->d_name, ".") == 0 || strcmp(dent->d_name, "..") == 0)
            continue;

        count = snprintf_grow(&buf, &buflen, "%s/%s", reply_dir, dent->d_name);
        if (count >= 0)
            remove(buf);
    }
    closedir(dirfd);

    if (reply_dir_created)
        remove(reply_dir);
    free(reply_dir);

    return retval;
}

/*
 * Internal helper function implementation.
 */
void graceful_exit(int signum)
{
    fprintf(stderr, "\nCaught signal %d.  Exiting.\n", signum);
    signal_caught = 1;
}

int init_signals(void)
{
    struct sigaction new_action;

    new_action.sa_handler = graceful_exit;
    sigemptyset(&new_action.sa_mask);
    new_action.sa_flags = 0;

    if (sigaction(SIGINT, &new_action, NULL) < 0)
        return -1;

    if (sigaction(SIGTERM, &new_action, NULL) < 0)
        return -1;

    return 0;
}

int init_comm(void)
{
    const char* cfgval;
    hcfg_t cfg = hcfg_zero;

    // Read a session message from the codegen plugin.
    if (mesg_recv(STDIN_FILENO, &mesg) < 1) {
        mesg.data.string = "Socket or unpacking error";
        return -1;
    }

    if (mesg.type != HMESG_SESSION || mesg.status != HMESG_STATUS_REQ) {
        mesg.data.string = "Invalid initial message";
        return -1;
    }

    if (hcfg_copy(&cfg, mesg.data.cfg) != 0) {
        mesg.data.string = "Could not copy config environment during init";
        return -1;
    }

    cfgval = hcfg_get(&cfg, CFGKEY_TMPDIR);
    if (!cfgval)
        cfgval = "";

    tmp_dir = stralloc(cfgval);
    if (!tmp_dir) {
        mesg.data.string = "Could not copy " CFGKEY_TMPDIR " config var";
        goto error;
    }

    sess_name = stralloc(mesg.state.space->name);
    if (!sess_name) {
        mesg.data.string = "Could not copy session name";
        goto error;
    }

    // Get the URL for the code generation server, and open a socket to it.
    cfgval = hcfg_get(&cfg, CFGKEY_SERVER_URL);
    if (!cfgval) {
        mesg.data.string = "Codegen server URL not specified";
        goto error;
    }

    if (url_parse(cfgval) < 0)
        goto error;

    if (scp_cmd && !hcfg_get(&cfg, CFGKEY_REPLY_URL)) {
        // User did not supply a reply URL.  Build one ourselves.
        if (reply_url_build() < 0)
            goto error;

        if (hcfg_set(&cfg, CFGKEY_REPLY_URL, buf) != 0) {
            mesg.data.string = "Could not set reply URL config val";
            goto error;
        }
    }
    else if (!scp_cmd) {
        if (hcfg_set(&cfg, CFGKEY_REPLY_URL, reply_dir) != 0) {
            mesg.data.string = "Could not set reply URL config val";
            goto error;
        }
    }

    if (hmesg_pack(&mesg) < 0) {
        mesg.data.string = "Could not allocate memory for initial message";
        goto error;
    }

    // Write the initial message for the code server.
    if (mesg_write(-1) < 0)
        goto error;

    hcfg_fini(&cfg);
    return 0;

  error:
    if (sess_name)
        free(sess_name);
    hcfg_fini(&cfg);
    return -1;
}

int url_parse(const char* url)
{
    struct stat sb;
    const char* ptr;
    const char* ssh_host;
    const char* ssh_port;
    const char* ssh_path;
    int host_len, port_len = 0;

    ptr = strstr(url, "//");
    if (!ptr)
        return -1;
    ptr += 2;

    if (strncmp("dir://", url, ptr - url) == 0) {
        // Codegen path is local.  Check that the path exists as a
        // directory, and use it for incoming and outgoing messages.
        //
        if (stat(ptr, &sb) < 0) {
            mesg.data.string = "Could not stat reply directory";
            return -1;
        }

        if (!S_ISDIR(sb.st_mode)) {
            mesg.data.string = "Reply path is not a directory";
            return -1;
        }

        reply_dir = stralloc(ptr);
        if (!reply_dir) {
            mesg.data.string = "Could not allocate memory for reply dir";
            return -1;
        }
    }
    else if (strncmp("ssh://", url, ptr - url) == 0) {
        // Codegen path exists on another machine accessible via ssh.
        // First, create a local directory to temporarily hold outgoing
        // messages, and receive incoming messages.
        //
        reply_dir = sprintf_alloc("%s/codegen-%s.%lx",
                                  tmp_dir, sess_name, time(NULL));
        if (!reply_dir) {
            mesg.data.string = "Could not allocate memory for temp directory";
            return -1;
        }

        if (mkdir(reply_dir, 0770) < 0) {
            mesg.data.string = "Could not create temp directory";
            return -1;
        }
        reply_dir_created = 1;

        // Now prepare the strings we'll use to launch scp for file
        // transfer.
        //
        url = ptr;

        ptr = strchr(url, '/');
        if (ptr) {
            ssh_host = url;
            host_len = ptr - url;
            ssh_path = ptr + 1;
        }
        else {
            ssh_host = url;
            host_len = strlen(ssh_host);
            ssh_path = "";
        }

        ptr = memchr(ssh_host, ':', host_len);
        if (ptr) {
            ssh_port = ptr + 1;
            port_len = host_len - (ssh_port - ssh_host);
            host_len -= port_len;
        }
        else {
            ssh_port = NULL;
        }

        scp_cmd = sprintf_alloc("scp %s%.*s",
                                ssh_port ? "-P "    : "",
                                ssh_port ? port_len : 0,
                                ssh_port ? ssh_port : "");
        scp_dst = sprintf_alloc("%.*s:%s", host_len, ssh_host, ssh_path);

        if (!scp_cmd || !scp_dst) {
            mesg.data.string = "Could not allocate memory for scp strings.";
            return -1;
        }
    }
    else {
        mesg.data.string = "Invalid use of codegen-helper";
        return -1;
    }
    return 0;
}

int reply_url_build(void)
{
    struct addrinfo  hints;
    struct addrinfo* info;
    char host[1024];
    int retval;

    if (gethostname(host, sizeof(host)) < 0) {
        mesg.data.string = "Could not determine local hostname";
        return -1;
    }
    host[sizeof(host) - 1] = '\0'; // Safety first!

    memset(&hints, 0, sizeof(hints));
    hints.ai_flags = AI_CANONNAME;

    retval = getaddrinfo(host, NULL, &hints, &info);
    if (retval < 0) {
        mesg.data.string = gai_strerror(retval);
        return -1;
    }

    retval = snprintf_grow(&buf, &buflen, "ssh://%s/%s",
                           info->ai_canonname, reply_dir);
    freeaddrinfo(info);
    if (retval < 0) {
        mesg.data.string = "Internal snprintf_grow error building reply URL";
        return -1;
    }
    return 0;
}

int mesg_read(int id)
{
    static const char in_filename[] = "code_complete";

    int msglen, fd, retries, retval;
    char* fullpath;
    char* newbuf;
    struct stat sb;
    struct timeval polltime;
    const char* errmsg;

    fullpath = sprintf_alloc("%s/%s.%d", reply_dir, in_filename, id);
    if (!fullpath) {
        mesg.data.string = "Internal sprintf_alloc error for flag filename";
        return -1;
    }

    retval = 0;
    retries = 3;

  top:
    fd = open(fullpath, O_RDONLY);
    if (fd == -1 && errno == ENOENT)
        goto cleanup;

    if (fd == -1) {
        errmsg = strerror(errno);
        goto retry;
    }

    // Obtain file size.
    if (fstat(fd, &sb) != 0) {
        errmsg = strerror(errno);
        goto retry;
    }

    msglen = sb.st_size + 1;
    if (mesg.recv_len < msglen) {
        newbuf = realloc(mesg.recv_buf, msglen);
        if (!newbuf) {
            mesg.data.string = "Could not allocate memory for message data.";
            retval = -1;
            goto cleanup;
        }
        mesg.recv_buf = newbuf;
        mesg.recv_len = msglen;
    }
    mesg.recv_buf[sb.st_size] = '\0';

    if (read_loop(fd, mesg.recv_buf, sb.st_size) != 0) {
        errmsg = "Error reading message file.";
        goto retry;
    }

    if (close(fd) < 0) {
        mesg.data.string = "Error closing code completion message file.";
        retval = -1;
        goto cleanup;
    }
    fd = -1;

    // Need to unpack to make sure have read a valid (complete) message.
    if (hmesg_unpack(&mesg) < 0) {
        errmsg = "Error decoding message file.";
        goto retry;
    }

    if (mesg_forward(STDIN_FILENO, &mesg) < 0) {
        mesg.data.string = "Error sending code completion message to session.";
        retval = -1;
        goto cleanup;
    }
    retval = 1;

    if (remove(fullpath) < 0) {
        mesg.data.string = strerror(errno);
        retval = -1;
    }

  cleanup:
    if (fd >= 0)
        close(fd);
    free(fullpath);
    return retval;

  retry:
    // Avoid the race condition where we attempt to read a message
    // before the remote code server scp process has completely
    // written the file.
    //
    // At some point, this retry method should probably be replaced
    // with file locks as a more complete solution.
    //
    close(fd);
    if (--retries) {
        polltime.tv_sec = 0;
        polltime.tv_usec = POLL_TIME;
        select(0, NULL, NULL, NULL, &polltime);
        goto top;
    }
    mesg.data.string = errmsg;
    free(fullpath);
    return -1;
}

int mesg_write(int id)
{
    static const char out_filename[] = "candidate";
    int fd, count;

    count = snprintf_grow(&buf, &buflen, "%s/%s.%d",
                          reply_dir, out_filename, id);
    if (count < 0) {
        mesg.data.string = "Internal snprintf_grow error for msg filename.";
        return -1;
    }

    fd = open(buf, O_WRONLY | O_CREAT, 0666);
    if (!fd) {
        mesg.data.string = strerror(errno);
        return -1;
    }

    unsigned short pkt_len;
    memcpy(&pkt_len, mesg.recv_buf + HMESG_LEN_OFFSET, HMESG_LEN_SIZE);
    pkt_len = ntohs(pkt_len);

    for (int i = 0; i < pkt_len; ++i) {
        if (mesg.recv_buf[i] == '\0')
            mesg.recv_buf[i] = '"';
    }

    if (write_loop(fd, mesg.recv_buf, pkt_len) != 0) {
        mesg.data.string = strerror(errno);
        return -1;
    }

    if (close(fd) < 0) {
        mesg.data.string = strerror(errno);
        return -1;
    }

    if (scp_cmd) {
        // Call scp to transfer the file.
        count = snprintf_grow(&buf, &buflen, "%s %s/%s.%d %s",
                              scp_cmd, reply_dir, out_filename, id, scp_dst);
        if (count < 0) {
            mesg.data.string = "Internal snprintf_grow error for scp command.";
            return -1;
        }

        if (system(buf) < 0) {
            mesg.data.string = strerror(errno);
            return -1;
        }

        // Remove the file from the local filesystem.
        count = snprintf_grow(&buf, &buflen, "%s/%s.%d",
                              reply_dir, out_filename, id);
        if (count < 0) {
            mesg.data.string = "Internal snprintf_grow error for scp command.";
            return -1;
        }
        if (remove(buf) < 0) {
            mesg.data.string = strerror(errno);
            return -1;
        }
    }
    return 0;
}

int read_loop(int fd, char* readbuf, int readlen)
{
    int count;

    while (readlen > 0) {
        count = read(fd, readbuf, readlen);
        if (count < 0 && errno == EINTR)
            continue;

        if (count <= 0) {
            mesg.data.string = strerror(errno);
            return -1;
        }

        readbuf += count;
        readlen -= count;
    }
    return 0;
}

int write_loop(int fd, char* writebuf, int writelen)
{
    int count;

    while (writelen > 0) {
        count = write(fd, writebuf, writelen);
        if (count < 0 && errno == EINTR)
            continue;

        if (count <= 0) {
            mesg.data.string = strerror(errno);
            return -1;
        }

        writebuf += count;
        writelen -= count;
    }
    return 0;
}
