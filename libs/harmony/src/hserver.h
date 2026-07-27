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
#ifndef __HSERVER_H__
#define __HSERVER_H__

#include "hmesg.h"
#include "hpoint.h"
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ilist {
    int* slot;
    int  len;
    int  cap;
} ilist_t;

typedef enum search_flags {
    FLAG_CONVERGED = 0x1,
    FLAG_PAUSED    = 0x2
} search_flag_t;

typedef struct http_log {
    struct timeval stamp;
    hpoint_t pt;
    double perf;
} http_log_t;

typedef struct sinfo {
    int id;

    // Best known search point and performance.
    hpoint_t best;
    double best_perf;

    // Client-related lists/
    ilist_t client;
    ilist_t request;

    // Fields used by the HTTP server.
    struct timeval start;
    hspace_t space;
    char* strategy;
    unsigned int flags;

    http_log_t* log;
    int log_len, log_cap;

    hpoint_t* fetched;
    int fetched_len, fetched_cap;
    int reported;
} sinfo_t;

extern sinfo_t* slist;
extern int slist_cap;

int  request_command(sinfo_t* sinfo, const char* command);
int  request_refresh(sinfo_t* sinfo);
int  request_setcfg(sinfo_t* sinfo, const char* key, const char* val);

#ifdef __cplusplus
}
#endif

#endif // __HSERVER_H__
