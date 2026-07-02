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

#ifndef __HMESG_H__
#define __HMESG_H__

#include "hspace.h"
#include "hcfg.h"
#include "hpoint.h"
#include "hperf.h"

/*
 * Message packet header layout.
 *
 * These fields use a binary encoding, whereas the message data
 * currently uses a text encoding.
 *
 *  0             15 16            31
 * |--------|--------|--------|--------|
 * |    HARMONY_MAGIC_BASE    |  Ver   |
 * |--------|--------|--------|--------|
 * | Message Length  |   Destination   |
 * |--------|--------|--------|--------|
 * |     Source      |  Message Data   |
 * |--------|--------|                 |
 * |       Message Data (cont.)        |
 * |                ...                |
 *
 */

// Offset and size of static header members in a byte-stream.
#define HMESG_MAGIC_OFFSET 0
#define HMESG_MAGIC_SIZE   4
#define HMESG_LEN_OFFSET   4
#define HMESG_LEN_SIZE     2
#define HMESG_PEEK_SIZE    6
#define HMESG_DEST_OFFSET  6
#define HMESG_DEST_SIZE    2
#define HMESG_SRC_OFFSET   8
#define HMESG_SRC_SIZE     2
#define HMESG_HEADER_SIZE 10

// Magic number for messages between the harmony server and its clients.
#define HMESG_OLDER_MAGIC  0x5261793a // Magic number for packets (pre v4.5).
#define HMESG_OLD_MAGIC    0x5261797c // Magic number for packets (pre v4.6.0).
#define HMESG_MAGIC_BASE   0x52617900 // Base for current magic number.
#define HMESG_MAGIC_VER          0x05 // Protocol version.
#define HMESG_MAGIC (HMESG_MAGIC_BASE | HMESG_MAGIC_VER)

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Data structures, definitions, and functions for Harmony messages.
 */
typedef enum {
    HMESG_UNKNOWN = 0x00,
    HMESG_SESSION, // Tuning session description.
    HMESG_JOIN,    // Client registration info.
    HMESG_GETCFG,  // Get search cfg key/value pair.
    HMESG_SETCFG,  // Set new search cfg key/value pair.
    HMESG_BEST,    // Retrieve best known point.
    HMESG_FETCH,   // Retrieve search space point to test.
    HMESG_REPORT,  // Report search space point performance.
    HMESG_COMMAND, // Request an action from the search task.

    HMESG_TYPE_MAX
} hmesg_type;

typedef enum {
    HMESG_STATUS_UNKNOWN = 0x00,
    HMESG_STATUS_REQ,  // Initial request message.
    HMESG_STATUS_OK,   // Request acknowledged.
    HMESG_STATUS_FAIL, // Request failed.
    HMESG_STATUS_BUSY, // Server could not respond in time.

    HMESG_STATUS_MAX
} hmesg_status;

/** \brief The hmesg_t structure.
 */
typedef struct hmesg {
    int dest;
    int src;
    hmesg_type type;
    hmesg_status status;

    // External state access pointers.
    struct hmesg_state {
        const hspace_t* space;
        const hpoint_t* best;
        const char*     client;
    } state;

    // External data access pointers.
    struct hmesg_data {
        const hcfg_t*   cfg;
        const hpoint_t* point;
        const hperf_t*  perf;
        const char*     string;
    } data;

    // Storage space for *_unpack() routines.
    hspace_t unpacked_space;
    hpoint_t unpacked_best;
    hcfg_t   unpacked_cfg;
    hpoint_t unpacked_point;
    hperf_t  unpacked_perf;

    char* recv_buf;
    int   recv_len;
    char* send_buf;
    int   send_len;
} hmesg_t;
#define HMESG_INITIALIZER {0}
extern const hmesg_t hmesg_zero;

void hmesg_fini(hmesg_t* mesg);
int  hmesg_forward(hmesg_t* mesg);
int  hmesg_pack(hmesg_t* mesg);
int  hmesg_unpack(hmesg_t* mesg);

#ifdef __cplusplus
}
#endif

#endif  // ifndef __HMESG_H__
