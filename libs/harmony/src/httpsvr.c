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
#include "hutil.h"
#include "hsockutil.h"
#include "hcfg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <fcntl.h>
#include <errno.h>
#include <limits.h>

#include <sys/socket.h>
#include <netinet/ip.h>
#include <netinet/tcp.h>

#define HTTP_ENDL "\r\n"

/*
 * The following opt_* definitions allow the compiler to possibly
 * optimize away a call to strlen().
 */
#define opt_sock_write(x, y) socket_write((x), (y), strlen(y))
#define opt_http_write(x, y) http_chunk_send((x), (y), strlen(y))
#define opt_strncmp(x, y) strncmp((x), (y), strlen(y))

unsigned int http_connection_limit = 32;

typedef enum {
    CONTENT_HTML,
    CONTENT_JAVASCRIPT,
    CONTENT_CSS
} content_t;

typedef struct {
    const char* filename;
    content_t type;
    char* buf;
    size_t buflen;
} memfile_t;

static memfile_t html_file[] = {
    // The overview html file must be first in this array.
    { "overview.cgi",                 CONTENT_HTML,       NULL, 0 },
    { "overview.js",                  CONTENT_JAVASCRIPT, NULL, 0 },
    { "session-view.cgi",             CONTENT_HTML,       NULL, 0 },
    { "session-view.js",              CONTENT_JAVASCRIPT, NULL, 0 },
    { "common.js",                    CONTENT_JAVASCRIPT, NULL, 0 },
    { "activeharmony.css",            CONTENT_CSS,        NULL, 0 },
    { "jquery.min.js",                CONTENT_JAVASCRIPT, NULL, 0 },
    { "jquery.flot.min.js",           CONTENT_JAVASCRIPT, NULL, 0 },
    { "jquery.flot.time.min.js",      CONTENT_JAVASCRIPT, NULL, 0 },
    { "jquery.flot.resize.min.js",    CONTENT_JAVASCRIPT, NULL, 0 },
    { "jquery.flot.selection.min.js", CONTENT_JAVASCRIPT, NULL, 0 },
    { "excanvas.min.js",              CONTENT_JAVASCRIPT, NULL, 0 },
    // A null entry must end this array.
    { NULL, CONTENT_HTML, NULL, 0 }
};

static const char status_200[] = "HTTP/1.1 200 OK" HTTP_ENDL;
static const char status_400[] = "HTTP/1.1 400 Bad Request" HTTP_ENDL;
static const char status_404[] = "HTTP/1.1 404 Not Found" HTTP_ENDL;
static const char status_500[] = "HTTP/1.1 500 Internal Error" HTTP_ENDL;
static const char status_501[] = "HTTP/1.1 501 Not Implemented" HTTP_ENDL;
static const char status_503[] = "HTTP/1.1 503 Service Unavailable" HTTP_ENDL;

static const char http_type_html[] = "Content-Type: text/html" HTTP_ENDL;
static const char http_type_text[] = "Content-Type: text/plain" HTTP_ENDL;
static const char http_type_js[] = "Content-Type: text/javascript" HTTP_ENDL;
static const char http_type_css[] = "Content-Type: text/css" HTTP_ENDL;

static const char http_headers[] =
    "Connection: Keep-Alive"     HTTP_ENDL
    "Keep-Alive: max=1024"       HTTP_ENDL
    "Cache-Control: no-cache"    HTTP_ENDL
    "Transfer-Encoding: chunked" HTTP_ENDL;

char sendbuf[8192];  // Static buffer used for outgoing HTTP data.
char recvbuf[10240]; // Static buffer used for incoming HTTP data.

/*
 * Internal helper function prototypes.
 */
const char* status_string(sinfo_t* sinfo);
char* uri_decode(char* buf);
sinfo_t* find_search(const char* name);
int http_request_handle(int fd, char* req);
char* http_request_recv(int fd, char* buf, int buflen, char** data);
int http_chunk_send(int fd, const char* data, int datalen);
int http_send_overview(int fd);
int http_send_init(int fd, sinfo_t* sinfo);
int http_send_refresh(int fd, sinfo_t* sinfo, const char* arg);
int report_append(char** buf, int* buflen, sinfo_t* sinfo,
                  struct timeval* tv, const hpoint_t* pt, const double perf);

int http_init(const char* basedir)
{
    int i;
    char* filename;

    for (i = 0; html_file[i].filename != NULL; ++i) {
        if (html_file[i].buf != NULL)
            file_unmap(html_file[i].buf, html_file[i].buflen);

        filename = sprintf_alloc("%s/libexec/http/%s", basedir,
                                 html_file[i].filename);
        if (!filename) {
            perror("Could not allocate temp memory for filename");
            goto error;
        }

        html_file[i].buf = file_map(filename, &html_file[i].buflen);
        if (!html_file[i].buf) {
            perror("Error mapping HTML server support file");
            goto error;
        }
        free(filename);
    }
    return 0;

  error:
    free(filename);
    return -1;
}

void http_send_error(int fd, int status, const char* arg)
{
    const char* status_line;
    const char* message = NULL;

    if (status == 400) {
        status_line = status_400;
        message = "The following request is malformed: ";

    } else if (status == 404) {
        status_line = status_404;
        message = "The requested URL could not be found: ";

    } else if (status == 501) {
        status_line = status_501;
        message = "The following request method has not been implemented: ";

    } else if (status == 503) {
        status_line = status_503;
        message = "The maximum number of HTTP connections has been exceeded.";

    } else {
        status_line = status_500;
        message = "An unknown status was passed to http_send_error().";
    }

    opt_sock_write(fd, status_line);
    opt_sock_write(fd, http_type_html);
    opt_sock_write(fd, http_headers);
    opt_sock_write(fd, HTTP_ENDL);
    opt_http_write(fd, "<html><head><title>");
    opt_http_write(fd, status_line +  9);
    opt_http_write(fd, "</title></head><body><h1>");
    opt_http_write(fd, status_line + 13);
    opt_http_write(fd, "</h1>");
    if (message || arg) {
        opt_http_write(fd, "<p>");
        if (message)
            opt_http_write(fd, message);
        if (arg)
            opt_http_write(fd, arg);
        opt_http_write(fd, "</p>");
    }
    opt_http_write(fd, "<hr>Active Harmony HTTP Server</body></html>");
    opt_http_write(fd, "");
}

int handle_http_socket(int fd)
{
    static char buf[4096];
    char* ptr = NULL;
    char* req;
    char* endp;

    errno = 0;
    req = http_request_recv(fd, buf, sizeof(buf), &ptr);
    while (req != NULL) {
        if (opt_strncmp(req, "GET ") == 0) {
            endp = strchr(req + 4, ' ');
            if (endp == NULL) {
                http_send_error(fd, 400, req); // "Bad Request" error.
            }
            else {
                *endp = '\0';
                req = uri_decode(req + 4);
                if (http_request_handle(fd, req) < 0)
                    http_send_error(fd, 404, req); // "Not Found" error.
            }
        }
        else if (opt_strncmp(req, "OPTIONS ") == 0 ||
                 opt_strncmp(req, "HEAD ")    == 0 ||
                 opt_strncmp(req, "POST ")    == 0 ||
                 opt_strncmp(req, "PUT ")     == 0 ||
                 opt_strncmp(req, "DELETE ")  == 0 ||
                 opt_strncmp(req, "TRACE ")   == 0 ||
                 opt_strncmp(req, "CONNECT ") == 0) {

            // "Unimplemented Function" error.
            http_send_error(fd, 501, req);
        }

        req = http_request_recv(fd, buf, sizeof(buf), &ptr);
    }

    if (errno == 0) {
        // No data returned, and no error reported.  Socket closed by peer.
        printf("[AH]: Closing socket %d (HTTP Socket)\n", fd);
        if (close(fd) < 0)
            printf("[AH]: Error closing HTTP socket\n");

        return -1;
    }

    return 0;
}

const char* status_string(sinfo_t* sinfo)
{
    if (sinfo->flags & FLAG_PAUSED)
        return "Paused";

    if (sinfo->flags & FLAG_CONVERGED)
        return "Converged";

    return "Searching";
}

/*
 * Internal helper function implementation.
 */
char* uri_decode(char* buf)
{
    char* head = buf;
    char* tail = buf;

    while (*head != '\0') {
        if (*head == '%') {
            if (isxdigit(head[1]) && isxdigit(head[2])) {
                unsigned int val;
                sscanf(head + 1, "%2x", &val);
                *(tail++) = (char)val;
                head += 3;
                continue;
            }
        }
        *(tail++) = *(head++);
    }
    *tail = '\0';

    return buf;
}

sinfo_t* find_search(const char* name)
{
    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i].id == -1)
            continue;

        if (strcmp(slist[i].space.name, name) == 0)
            return &slist[i];
    }
    return NULL;
}

int http_request_handle(int fd, char* req)
{
    sinfo_t* sinfo = NULL;
    char* search_name;
    char* arg = NULL;
    int i;

    search_name = strchr(req, '?');
    if (search_name) {
        *(search_name++) = '\0';

        arg = strchr(search_name, '&');
        if (arg)
            *(arg++) = '\0';

        sinfo = find_search(search_name);
    }

    if (strcmp(req, "/") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_html);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);
        http_chunk_send(fd, html_file[0].buf, html_file[0].buflen);
        opt_http_write(fd, "");

        return 0;
    }
    else if (strcmp(req, "/session-list") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        http_send_overview(fd);
        opt_http_write(fd, "");
        return 0;
    }
    else if (strcmp(req, "/init") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        if (sinfo) {
            http_send_init(fd, sinfo);
            opt_http_write(fd, "");
            return 0;
        }
        opt_http_write(fd, "FAIL");
        opt_http_write(fd, "");
        return 0;
    }
    else if (strcmp(req, "/refresh") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        if (sinfo && request_refresh(sinfo) == 0) {
            http_send_refresh(fd, sinfo, arg);
            opt_http_write(fd, "");
            return 0;
        }
        opt_http_write(fd, "FAIL");
        opt_http_write(fd, "");
        return 0;
    }
    else if (strcmp(req, "/kill") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        if (sinfo && request_command(sinfo, "kill") == 0) {
            opt_http_write(fd, "OK");
            opt_http_write(fd, "");
            return 0;
        }
        opt_http_write(fd, "FAIL");
        opt_http_write(fd, "");
        return 0;
    }
    else if (strcmp(req, "/pause") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        if (sinfo) {
            request_setcfg(sinfo, CFGKEY_PAUSED, "1");
            opt_http_write(fd, "OK");
            opt_http_write(fd, "");
            return 0;
        }
        opt_http_write(fd, "FAIL");
        opt_http_write(fd, "");
        return 0;
    }
    else if (strcmp(req, "/resume") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        if (sinfo) {
            request_setcfg(sinfo, CFGKEY_PAUSED, "0");
            opt_http_write(fd, "OK");
            opt_http_write(fd, "");
            return 0;
        }
        opt_http_write(fd, "FAIL");
        opt_http_write(fd, "");
        return 0;
    }
    else if (strcmp(req, "/restart") == 0) {
        opt_sock_write(fd, status_200);
        opt_sock_write(fd, http_type_text);
        opt_sock_write(fd, http_headers);
        opt_sock_write(fd, HTTP_ENDL);

        if (sinfo) {
            if (arg && strstr(arg, "_=") != arg) {
                request_setcfg(sinfo, CFGKEY_INIT_POINT, arg);
            }
            request_command(sinfo, "restart");
            opt_http_write(fd, "OK");
            opt_http_write(fd, "");
            return 0;
        }
        opt_http_write(fd, "FAIL");
        opt_http_write(fd, "");
        return 0;
    }

    // If request is not handled by any special cases above,
    // look for a known html file corresponding to the request.
    //
    for (i = 0; html_file[i].filename != NULL; ++i) {
        if (strcmp(req + 1, html_file[i].filename) == 0) {
            opt_sock_write(fd, status_200);
            switch (html_file[i].type) {
            case CONTENT_HTML:
                opt_sock_write(fd, http_type_html);
                break;
            case CONTENT_JAVASCRIPT:
                opt_sock_write(fd, http_type_js);
                break;
            case CONTENT_CSS:
                opt_sock_write(fd, http_type_css);
                break;
            }
            opt_sock_write(fd, http_headers);
            opt_sock_write(fd, HTTP_ENDL);

            http_chunk_send(fd, html_file[i].buf, html_file[i].buflen);
            opt_http_write(fd, "");
            return 0;
        }
    }

    return -1;
}

int http_chunk_send(int fd, const char* data, int datalen)
{
    int n;
    char buf[11];

    n = snprintf(buf, sizeof(buf), "%x" HTTP_ENDL, datalen);

    // Relying on TCP buffering to keep this somewhat efficient.
    socket_write(fd, buf, n);
    if (datalen > 0)
        n += socket_write(fd, data, datalen);
    return n + opt_sock_write(fd, HTTP_ENDL);
}

char* http_request_recv(int fd, char* buf, int buflen, char** data)
{
    const char delim[] = HTTP_ENDL HTTP_ENDL;
    char* retval;
    char* split;
    int len = 0, recvlen;

    if (!*data) {
        *data = buf;
        *buf = '\0';
    }

    while ( (split = strstr(*data, delim)) == NULL) {
        len = strlen(*data);

        if (*data == buf && len == (buflen-1)) {
            // Buffer overflow.  Return truncated request.
            break;
        }

        if (*data != buf) {
            if (len > 0) {
                // Move existing data to the front of buffer.
                if (len < *data - buf) {
                    strcpy(buf, *data);

                } else {
                    // Memory region overlaps.  Can't use strcpy().
                    for (len = 0; (*data)[len] != '\0'; ++len)
                        buf[len] = (*data)[len];
                    buf[len] = '\0';
                }
            }
            *data = buf;
        }

        // Read more data from socket into the remaining buffer space.
        recvlen = recv(fd, buf + len, buflen - (len + 1), MSG_DONTWAIT);
        if (recvlen < 0) {
            return NULL;
        } else if (recvlen == 0) {
            if (len == 0) {
                // No data in buffer, and connection was closed.
                return NULL;
            }
            break;
        }

        len += recvlen;
        buf[len] = '\0';
    }

    retval = *data;
    if (split) {
        memset(split, '\0', strlen(delim));
        *data = split + strlen(delim);

    } else {
        // Socket closed while buffer not empty.  Return remaining data
        *data += len;
    }
    // /* DEBUG */ printf("HTTP REQ: %s", retval);
    return retval;
}

int http_send_overview(int fd)
{
    char* buf;
    char* ptr;
    int buflen, count, total;

    sendbuf[0] = '\0';
    buf = sendbuf;
    buflen = sizeof(sendbuf);
    total = 0;

    for (int i = 0; i < slist_cap; ++i) {
        if (slist[i].id < 0)
            continue;

        ptr = buf;
        count = snprintf_serial(&buf, &buflen, "%s:%ld%03ld:%d:%d:",
                                slist[i].space.name,
                                slist[i].start.tv_sec,
                                slist[i].start.tv_usec/1000,
                                slist[i].client.len,
                                slist[i].reported);
        if (count < 0) return -1;
        total += count;

        if (!slist[i].best.id) {
            count = snprintf_serial(&buf, &buflen, "&lt;unknown&gt;");
            if (count < 0) return -1;
            total += count;
        }
        else {
            hrange_t* dim = slist[i].space.dim;
            for (int j = 0; j < slist[i].best.len; ++j) {
                const hval_t* v = &slist[i].best.term[j];

                if (j > 0) {
                    count = snprintf_serial(&buf, &buflen, " ");
                    if (count < 0) return -1;
                    total += count;
                }

                switch (dim[j].type) {
                case HVAL_INT:
                    count = snprintf_serial(&buf, &buflen, "%ld", v->value.i);
                    if (count < 0) return -1;
                    total += count;
                    break;
                case HVAL_REAL:
                    count = snprintf_serial(&buf, &buflen, "%lf", v->value.r);
                    if (count < 0) return -1;
                    total += count;
                    break;
                case HVAL_STR:
                    count = snprintf_serial(&buf, &buflen, "%s", v->value.s);
                    if (count < 0) return -1;
                    total += count;
                    break;
                default:
                    return -1;
                }
            }
        }

        count = snprintf_serial(&buf, &buflen, "|");
        if (count < 0) return -1;
        total += count;

        if (total >= sizeof(sendbuf)) {
            // Corner Case: No room for any search data.  Error out.
            if (ptr == sendbuf) {
                opt_http_write(fd, "FAIL");
                return -1;
            }

            // Make room in buffer by sending the contents thus far.
            *ptr = '\0';
            opt_http_write(fd, sendbuf);

            // Reset the loop variables and try again.
            buf = sendbuf;
            buflen = sizeof(sendbuf);
            total = 0;
            --i;
        }
    }

    if (sendbuf[0])
        opt_http_write(fd, sendbuf);
    return 0;
}

int http_send_init(int fd, sinfo_t* sinfo)
{
    char* buf = sendbuf;
    int buflen = sizeof(sendbuf), total = 0;
    int count;

    count = snprintf_serial(&buf, &buflen, "space:");
    if (count < 0)
        goto error;
    total += count;

    for (int i = 0; i < sinfo->space.len; ++i) {
        char type;
        hrange_t* range = &sinfo->space.dim[i];

        if (i > 0) {
            count = snprintf_serial(&buf, &buflen, ",");
            if (count < 0)
                goto error;
            total += count;
        }

        switch (range->type) {
        case HVAL_INT:  type = 'i'; break;
        case HVAL_REAL: type = 'r'; break;
        case HVAL_STR:  type = 's'; break;
        default: goto error;
        }

        count = snprintf_serial(&buf, &buflen, "%s;%c", range->name, type);
        if (count < 0)
            goto error;
        total += count;

        if (range->type == HVAL_STR) {
            for (int j = 0; j < range->bounds.e.len; ++j) {
                count = snprintf_serial(&buf, &buflen, ";%s",
                                        range->bounds.e.set[j]);
                if (count < 0)
                    goto error;
                total += count;
            }
        }
    }

    if (sinfo->strategy)
        count = snprintf_serial(&buf, &buflen, "|strat:%s", sinfo->strategy);
    else
        count = snprintf_serial(&buf, &buflen, "|strat:<unknown>");
    if (count < 0)
        goto error;
    total += count;

    if (total >= sizeof(sendbuf))
        goto error;

    opt_http_write(fd, sendbuf);
    return 0;

  error:
    opt_http_write(fd, "FAIL");
    return -1;
}

int http_send_refresh(int fd, sinfo_t* sinfo, const char* arg)
{
    char* ptr;
    char* buf = sendbuf;
    int idx = 0, buflen = sizeof(sendbuf), total = 0;
    int i, count;
    struct timeval tv;

    if (arg)
        idx = atoi(arg);

    if (gettimeofday(&tv, NULL) != 0)
        goto error;

    count = snprintf_serial(&buf, &buflen,
                            "time:%ld%03ld|status:%s|clients:%d|index:%d",
                            tv.tv_sec, tv.tv_usec/1000,
                            status_string(sinfo),
                            sinfo->client.len,
                            idx);
    if (count < 0)
        goto error;
    total += count;

    count = snprintf_serial(&buf, &buflen, "|best:");
    if (count < 0)
        goto error;
    total += count;

    count = report_append(&buf, &buflen, sinfo, NULL,
                          &sinfo->best, sinfo->best_perf);
    if (count < 0)
        goto error;
    total += count;

    for (i = idx; i < sinfo->log_len; ++i) {
        ptr = buf;

        count = snprintf_serial(&buf, &buflen, "|trial:");
        if (count < 0)
            goto error;
        total += count;

        count = report_append(&buf, &buflen, sinfo, &sinfo->log[i].stamp,
                              &sinfo->log[i].pt, sinfo->log[i].perf);
        if (count < 0)
            goto error;
        total += count;

        if (total >= sizeof(sendbuf)) {
            // Corner Case: No room for any search data.  Error out.
            if (ptr == sendbuf)
                goto error;

            // Make room in buffer by sending the contents thus far.
            *ptr = '\0';
            opt_http_write(fd, sendbuf);

            // Reset the loop variables and try again.
            buf = sendbuf;
            buflen = sizeof(sendbuf);
            total = 0;
            --i;
        }
    }

    opt_http_write(fd, sendbuf);
    return 0;

  error:
    opt_http_write(fd, "FAIL");
    return -1;
}

int report_append(char** buf, int* buflen, sinfo_t* sinfo,
                  struct timeval* tv, const hpoint_t* pt, const double perf)
{
    int count, total = 0;

    if (tv) {
        count = snprintf_serial(buf, buflen, "%ld%03ld,",
                                tv->tv_sec, tv->tv_usec/1000);
        if (count < 0)
            return -1;
        total += count;
    }

    if (!pt->id) {
        for (int i = 0; i < sinfo->space.len; ++i) {
            count = snprintf_serial(buf, buflen, ",");
            if (count < 0) return -1;
            total += count;
        }
    }
    else {
        hrange_t* dim = sinfo->space.dim;
        for (int i = 0; i < sinfo->space.len; ++i) {
            const hval_t* v = &pt->term[i];

            switch (dim[i].type) {
            case HVAL_INT:
                count = snprintf_serial(buf, buflen, "%ld,", v->value.i);
                if (count < 0) return -1;
                total += count;
                break;
            case HVAL_REAL:
                count = snprintf_serial(buf, buflen, "%lf,", v->value.r);
                if (count < 0) return -1;
                total += count;
                break;
            case HVAL_STR: {
                unsigned long index = hrange_index(&sinfo->space.dim[i], v);
                count = snprintf_serial(buf, buflen, "%ld,", index);
                if (count < 0) return -1;
                total += count;
                break;
            }
            default:
                return -1;
            }

            if (count < 0) return -1;
            total += count;
        }

        count = snprintf_serial(buf, buflen, "%lf", perf);
        if (count < 0) return -1;
        total += count;
    }

    return total;
}
