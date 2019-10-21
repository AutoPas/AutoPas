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

#include "hutil.h"

#include <stdio.h>
#include <fcntl.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <limits.h>
#include <errno.h>
#include <sys/mman.h>
#include <sys/stat.h>

int file_exists(const char* filename)
{
    struct stat sb;
    return (stat(filename, &sb) == 0 && S_ISREG(sb.st_mode));
}

void* file_map(const char* filename, size_t* size)
{
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Error on open(%s): %s\n", filename, strerror(errno));
        return NULL;
    }

    // Obtain file size.
    struct stat sb;
    void* retval = NULL;
    if (fstat(fd, &sb) != 0) {
        fprintf(stderr, "Error on fstat(%s): %s\n", filename, strerror(errno));
    }
    else {
        retval = mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE,
                      MAP_PRIVATE, fd, 0);
        if (retval == MAP_FAILED) {
            fprintf(stderr, "Error on mmap(%s): %s\n",
                    filename, strerror(errno));
        }
    }

    if (close(fd) != 0) {
        fprintf(stderr, "Warning: Ignoring error on close(%s): %s\n",
                filename, strerror(errno));
    }

    *size = sb.st_size;
    return retval;
}

void file_unmap(void* buf, size_t size)
{
    if (munmap(buf, size) != 0)
        fprintf(stderr, "Ignoring error on munmap(): %s\n", strerror(errno));
}

int valid_id(const char* key, int len)
{
    int span = -1;
    sscanf(key, "%*[0-9A-Z_a-z]%n", &span);
    return len > 0 && len == span && !isdigit(key[0]);
}

int line_unquote(char* buf)
{
    // Skip leading whitespace.
    char* src = buf;
    while (isspace(*src) && *src != '\n') ++src;

    // Find first separating '=' character.
    char* sep = strchr(buf, '=');
    while (buf < sep && isspace(*(sep - 1))) --sep;

    char* dst = buf;
    char  quote = '\0';
    while (*src) {
        // Remove spaces around the separating '=' character.
        if (src == sep) {
            *(dst++) = '=';
            while ((isspace(*src) && *src != '\n') || *src == '=') ++src;
            continue;
        }

        // Remove unquoted backslash characters.
        if (!quote) {
            if      (*src == '\'') quote = '\'';
            else if (*src ==  '"') quote =  '"';
            else if (*src ==  '#') break;
            else if (*src == '\n') break;
            else if (*src == '\\' && *(src + 1)) ++src;
        }
        else {
            if      (*src == '\\')  *(dst++) = *(src++);
            else if (*src == quote) quote = '\0';
        }

        if (*src)
            *(dst++) = *(src++);
    }
    if (quote)
        return -1;

    // Trim trailing whitespace.
    while (buf < dst && isspace(*(dst - 1))) --dst;

    // End line with '\0'.
    *dst = '\0';
    return 0;
}

int line_length(char* buf, int* linecnt)
{
    int  i = 0;
    char quote = '\0';

    *linecnt = 1;
    while (buf[i]) {
        if (buf[i] == '\\' && *(buf + 1)) ++i;
        else if (!quote) {
            if      (buf[i] == '\'') quote = '\'';
            else if (buf[i] ==  '"') quote =  '"';
            else if (buf[i] == '\n') break;
            else if (buf[i] ==  '#') { i += strcspn(&buf[i], "\n"); break; }
        }
        else if (buf[i] == quote) quote = '\0';

        if (buf[i] == '\n') ++(*linecnt);
        if (buf[i] != '\0') ++i;
    }


    return i;
}

int file_read_line(FILE* fp, char** buf, int* cap,
                   char** line, char** end, const char** errptr)
{
    int linecnt;
    const char* errstr;

    // Allocate an initial buffer if necessary.
    if (*cap == 0) {
        *cap = 1024;
        *buf = malloc(*cap * sizeof(**buf));

        if (!*buf) {
            errstr = "Could not allocate configuration parse buffer";
            goto error;
        }
        **buf = '\0';
    }

    *line = *buf;
    if (*end)
        *line = *end;

    // Loop until end of line is found (or end of file).
    int len;
    while (1) {
        len = line_length(*line, &linecnt);

        if ((*line)[len] == '\0' && !feof(fp)) {
            // Line is incomplete, but more data is available.

            // Move partial line to beginning of buffer.
            memmove(*buf, *line, len);

            int remain = *cap - len - 2;
            if (remain < 1) {
                // Extend the line buffer.
                if (array_grow(buf, cap, sizeof(**buf)) != 0) {
                    errstr = "Could not grow config parsing buffer";
                    goto error;
                }
                remain = *cap - len - 2;
            }

            // Read more data into the line buffer.
            int tail = len + fread(*buf + len, sizeof(**buf), remain, fp);
            (*buf)[tail] = '\0';
            *line = *buf;
        }
        else break; // Loop exit.
    }

    // Mark where line buffer processing has ended.
    if ((*line)[len] == '\n') ++len;
    *end = *line + len;

    // File and line buffer exhausted.
    if (*line == *end)
        return 0;

    if (line_unquote(*line) != 0) {
        errstr = "Non-terminated quote detected";
        goto error;
    }
    return linecnt;

  error:
    if (errptr)
        *errptr = errstr;
    return -1;
}

int array_grow(void* oldbuf, int* oldcap, int size)
{
    char* newbuf;
    int   newcap = 8;

    if (*oldcap >= newcap)
        newcap = *oldcap << 1;

    if (*(void**)oldbuf == NULL) {
        newbuf = calloc(newcap, size);
        if (!newbuf)
            return -1;
    }
    else {
        newbuf = realloc(*(void**)oldbuf, newcap * size);
        if (!newbuf)
            return -1;

        memset(newbuf + (*oldcap * size), 0, (newcap - *oldcap) * size);
    }
    *(void**)oldbuf = newbuf;
    *oldcap = newcap;

    return 0;
}

char* search_path(const char* filename)
{
    // XXX - Not re-entrant.
    static char* fullpath = NULL;
    static int pathlen = 0;

    char* path;
    char* pend;
    char* newbuf;
    int count;

    path = getenv("PATH");
    if (path == NULL)
        return NULL;

    pend = path;
    while (*pend != '\0') {
        pend = strchr(path, ':');
        if (!pend)
            pend = path + strlen(path);

      retry:
        count = snprintf(fullpath, pathlen, "%.*s/%s",
                         (int)(pend - path), path, filename);
        if (pathlen <= count) {
            newbuf = realloc(fullpath, count + 1);
            if (!newbuf) return NULL;

            fullpath = newbuf;
            pathlen = count + 1;
            goto retry;
        }

        if (file_exists(fullpath))
            return fullpath;

        path = pend + 1;
    }
    return NULL;
}

char* stralloc(const char* in)
{
    char* out;
    if (!in)
        return NULL;

    out = malloc(sizeof(char) * (strlen(in) + 1));
    if (out != NULL)
        strcpy(out, in);

    return out;
}

char* sprintf_alloc(const char* fmt, ...)
{
    va_list ap;
    int count;
    char* retval;

    va_start(ap, fmt);
    count = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);

    if (count < 0) return NULL;
    retval = malloc(count + 1);
    if (!retval) return NULL;

    va_start(ap, fmt);
    vsnprintf(retval, count + 1, fmt, ap);
    va_end(ap);

    return retval;
}

int snprintf_grow(char** buf, int* buflen, const char* fmt, ...)
{
    va_list ap;
    int count;
    char* newbuf;

  retry:
    va_start(ap, fmt);
    count = vsnprintf(*buf, *buflen, fmt, ap);
    va_end(ap);

    if (count < 0)
        return -1;

    if (*buflen <= count) {
        newbuf = realloc(*buf, count + 1);
        if (!newbuf)
            return -1;
        *buf = newbuf;
        *buflen = count + 1;
        goto retry;
    }

    return count;
}

int snprintf_serial(char** buf, int* buflen, const char* fmt, ...)
{
    va_list ap;
    int count;

    va_start(ap, fmt);
    count = vsnprintf(*buf, *buflen, fmt, ap);
    va_end(ap);

    if (count < 0)
        return -1;

    *buflen -= count;
    if (*buflen < 0)
        *buflen = 0;
    else
        *buf += count; // Only advance *buf if the write was successful.

    return count;
}

int printstr_serial(char** buf, int* buflen, const char* str)
{
    if (!str) return snprintf_serial(buf, buflen, " 0\"0\"");
    return snprintf_serial(buf, buflen, " %u\"%s\"", strlen(str), str);
}

int scanstr_serial(const char** str, char* buf)
{
    int count;
    unsigned int len;

    if (sscanf(buf, " %u\"%n", &len, &count) < 1)
        goto invalid;
    buf += count;

    if (len == 0 && *buf == '0') {
        len = 1;
        *str = NULL;
    }
    else if (len < strlen(buf) && buf[len] == '\"') {
        buf[len] = '\0';
        *str = buf;
    }
    else goto invalid;

    return count + len + 1;

  invalid:
    errno = EINVAL;
    return -1;
}

int unquote_string(const char* buf, char** token, const char** errptr)
{
    const char* src;
    int cap = 0;

    while (1) {
        src = buf;
        // Skip leading whitespace.
        while (isspace(*src)) ++src;

        char* dst = *token;
        int   len = cap;
        char  quote = '\0';
        while (*src) {
            if (!quote) {
                if      (*src == '\'')  { quote = '\''; ++src; continue; }
                else if (*src ==  '"')  { quote =  '"'; ++src; continue; }
                else if (*src ==  ',')  { break; }
                else if (isspace(*src)) { break; }
            }
            else {
                if      (*src == '\\')  { ++src; }
                else if (*src == quote) { quote = '\0'; ++src; continue; }
            }

            if (len-- > 0)
                *(dst++) = *(src++);
            else
                ++src;
        }
        if (len-- > 0)
            *dst = '\0';

        if (quote != '\0') {
            if (errptr)
                *errptr = "Non-terminated quote detected";
            return -1;
        }

        if (len < -1) {
            // Token buffer size is -len;
            cap += -len;
            *token = malloc(cap * sizeof(**token));
            if (!*token) {
                if (errptr)
                    *errptr = "Could not allocate memory in unquote_string()";
                return -1;
            }
        }
        else if (len == -1) {
            // Empty token.
            *token = NULL;
            break;
        }
        else break; // Loop exit.
    }
    return src - buf;
}
