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

#include "hmesg.h"
#include "hval.h"
#include "hutil.h"
#include "hsockutil.h"

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <arpa/inet.h>

const hmesg_t hmesg_zero = HMESG_INITIALIZER;

/*
 * Internal helper function prototypes.
 */
static int pack_header(char* buf, int length, const hmesg_t* mesg);
static int unpack_header(hmesg_t* mesg);
static int pack_state(char** buf, int* buflen, const hmesg_t* mesg);
static int unpack_state(hmesg_t* mesg, char* buf);
static int pack_data(char** buf, int* buflen, const hmesg_t* mesg);
static int unpack_data(hmesg_t* mesg, char* buf);

/*
 * To avoid excessive memory allocation, the *_unpack() routines build
 * non-standard versions of their structures by taking advantage of the
 * hmesg_t internal buffers.
 *
 * As a result, we use the *_scrub() routines here instead of their
 * corresponding *_fini() routine.
 */
void hmesg_fini(hmesg_t* mesg)
{
    hspace_scrub(&mesg->unpacked_space);
    hpoint_scrub(&mesg->unpacked_best);
    hcfg_scrub(&mesg->unpacked_cfg);
    hpoint_scrub(&mesg->unpacked_point);
    hperf_fini(&mesg->unpacked_perf); // No scrub routine for hperf_t.

    free(mesg->recv_buf);
    free(mesg->send_buf);
}

int hmesg_forward(hmesg_t* mesg)
{
    return pack_header(mesg->recv_buf, 0, mesg);
}

int hmesg_pack(hmesg_t* mesg)
{
    int count, total;

    while (1) {
        char* buf = mesg->send_buf;
        int   buflen = mesg->send_len;

        // Leave room for the packet header.
        buf += HMESG_HEADER_SIZE;
        buflen -= HMESG_HEADER_SIZE;
        if (buflen < 0)
            buflen = 0;
        total = HMESG_HEADER_SIZE;

        const char* type_str;
        switch (mesg->type) {
        case HMESG_UNKNOWN: type_str = "UNK"; break;
        case HMESG_SESSION: type_str = "SES"; break;
        case HMESG_JOIN:    type_str = "JOI"; break;
        case HMESG_GETCFG:  type_str = "QRY"; break;
        case HMESG_SETCFG:  type_str = "INF"; break;
        case HMESG_BEST:    type_str = "BST"; break;
        case HMESG_FETCH:   type_str = "FET"; break;
        case HMESG_REPORT:  type_str = "REP"; break;
        case HMESG_COMMAND: type_str = "CMD"; break;
        default:
            fprintf(stderr, "Error during hmesg_pack():"
                    "Message type (%d) is invalid\n", (int) mesg->type);
            goto invalid;
        }

        const char* status_str;
        switch (mesg->status) {
        case HMESG_STATUS_REQ:  status_str = "REQ"; break;
        case HMESG_STATUS_OK:   status_str = "ACK"; break;
        case HMESG_STATUS_FAIL: status_str = "ERR"; break;
        case HMESG_STATUS_BUSY: status_str = "BSY"; break;
        default:
            fprintf(stderr, "Error during hmesg_pack():"
                    "Message status (%d) is invalid\n", (int) mesg->status);
            goto invalid;
        }
        count = snprintf_serial(&buf, &buflen, " %s %s", type_str, status_str);
        if (count < 0) goto error;
        total += count;

        if (mesg->status == HMESG_STATUS_FAIL) {
            count = printstr_serial(&buf, &buflen, mesg->data.string);
            if (count < 0) goto error;
            total += count;
        }
        else {
            count = pack_state(&buf, &buflen, mesg);
            if (count < 0) goto error;
            total += count;

            count = pack_data(&buf, &buflen, mesg);
            if (count < 0) goto error;
            total += count;
        }

        if (total >= mesg->send_len) {
            buf = realloc(mesg->send_buf, total + 1);
            if (!buf) goto error;

            mesg->send_buf = buf;
            mesg->send_len = total + 1;
        }
        else break;
    }

    // Now that message length is known, write packet header.
    if (pack_header(mesg->send_buf, total, mesg) != 0)
        goto error;

    return total;

  invalid:
    errno = EINVAL;
  error:
    return -1;
}

int hmesg_unpack(hmesg_t* mesg)
{
    int count, total;
    char* buf = mesg->recv_buf;

    total = unpack_header(mesg);
    if (total < 0)
        goto invalid;

    char type_str[4];
    char status_str[4];
    if (sscanf(buf + total, " %3s %3s%n", type_str, status_str, &count) < 2)
        goto invalid;
    total += count;

    if      (strcmp(type_str, "UNK") == 0) mesg->type = HMESG_UNKNOWN;
    else if (strcmp(type_str, "SES") == 0) mesg->type = HMESG_SESSION;
    else if (strcmp(type_str, "JOI") == 0) mesg->type = HMESG_JOIN;
    else if (strcmp(type_str, "QRY") == 0) mesg->type = HMESG_GETCFG;
    else if (strcmp(type_str, "INF") == 0) mesg->type = HMESG_SETCFG;
    else if (strcmp(type_str, "BST") == 0) mesg->type = HMESG_BEST;
    else if (strcmp(type_str, "FET") == 0) mesg->type = HMESG_FETCH;
    else if (strcmp(type_str, "REP") == 0) mesg->type = HMESG_REPORT;
    else if (strcmp(type_str, "CMD") == 0) mesg->type = HMESG_COMMAND;
    else goto invalid;

    if      (strcmp(status_str, "REQ") == 0) mesg->status = HMESG_STATUS_REQ;
    else if (strcmp(status_str, "ACK") == 0) mesg->status = HMESG_STATUS_OK;
    else if (strcmp(status_str, "ERR") == 0) mesg->status = HMESG_STATUS_FAIL;
    else if (strcmp(status_str, "BSY") == 0) mesg->status = HMESG_STATUS_BUSY;
    else goto invalid;

    if (mesg->status == HMESG_STATUS_FAIL) {
        count = scanstr_serial(&mesg->data.string, buf + total);
        if (count < 0) goto error;
        total += count;
    }
    else {
        count = unpack_state(mesg, buf + total);
        if (count < 0) goto error;
        total += count;

        count = unpack_data(mesg, buf + total);
        if (count < 0) goto error;
        total += count;
    }
    return total;

  invalid:
    errno = EINVAL;
  error:
    return -1;
}

/*
 * Internal helper function implementation.
 */
int pack_header(char* buf, int length, const hmesg_t* mesg)
{
    if (length && (length < HMESG_HEADER_SIZE || length > 0xFFFF)) {
        fprintf(stderr, "Error during pack_header():"
                "Message length (%d) is out of range [%d, 65535)\n",
                HMESG_HEADER_SIZE, length);
        return -1;
    }

    if (mesg->dest < -1 || mesg->dest >= 0xFFFF) {
        fprintf(stderr, "Error during pack_header():"
                "Destination (%d) is out of range [-1, 65534]\n", mesg->dest);
        return -1;
    }

    if (mesg->src < -1 || mesg->src >= 0xFFFF) {
        fprintf(stderr, "Error during pack_header():"
                "Source (%d) is out of range [-1, 65534]\n", mesg->src);
        return -1;
    }

    unsigned int   pkt_magic = htonl(HMESG_MAGIC);
    unsigned short pkt_len   = htons((unsigned short) length);
    unsigned short pkt_dest  = htons((unsigned short) mesg->dest);
    unsigned short pkt_src   = htons((unsigned short) mesg->src);

    if (length) {
        memcpy(buf + HMESG_MAGIC_OFFSET, &pkt_magic, HMESG_MAGIC_SIZE);
        memcpy(buf + HMESG_LEN_OFFSET,   &pkt_len,   HMESG_LEN_SIZE);
    }
    memcpy(buf + HMESG_DEST_OFFSET,  &pkt_dest,  HMESG_DEST_SIZE);
    memcpy(buf + HMESG_SRC_OFFSET,   &pkt_src,   HMESG_SRC_SIZE);

    return 0;
}

int unpack_header(hmesg_t* mesg)
{
    unsigned int   pkt_magic;
    unsigned short pkt_dest;
    unsigned short pkt_src;

    memcpy(&pkt_magic, mesg->recv_buf + HMESG_MAGIC_OFFSET, HMESG_MAGIC_SIZE);
    memcpy(&pkt_dest,  mesg->recv_buf + HMESG_DEST_OFFSET,  HMESG_DEST_SIZE);
    memcpy(&pkt_src,   mesg->recv_buf + HMESG_SRC_OFFSET,   HMESG_SRC_SIZE);

    if (ntohl(pkt_magic) != HMESG_MAGIC)
        return -1;

    mesg->dest = ntohs(pkt_dest);
    if (mesg->dest >= 0xFFFF)
        mesg->dest = -1;


    mesg->src = ntohs(pkt_src);
    if (mesg->src >= 0xFFFF)
        mesg->src = -1;

    return HMESG_HEADER_SIZE;
}

int pack_state(char** buf, int* buflen, const hmesg_t* mesg)
{
    int count, total = 0;

    if (mesg->type == HMESG_SESSION ||
        mesg->type == HMESG_JOIN ||
        mesg->type == HMESG_UNKNOWN)
    {
        // Session state doesn't exist for these messages.
        return 0;
    }

    switch (mesg->status) {
    case HMESG_STATUS_REQ:
        count = snprintf_serial(buf, buflen, " state:%u %u",
                                mesg->state.space->id, mesg->state.best->id);
        if (count < 0) return -1;
        total += count;

        count = printstr_serial(buf, buflen, mesg->state.client);
        if (count < 0) return -1;
        total += count;
        break;

    case HMESG_STATUS_OK:
    case HMESG_STATUS_BUSY:
        count = snprintf_serial(buf, buflen, " state:");
        if (count < 0) return -1;
        total += count;

        count = hspace_pack(buf, buflen, mesg->state.space);
        if (count < 0) return -1;
        total += count;

        count = hpoint_pack(buf, buflen, mesg->state.best);
        if (count < 0) return -1;
        total += count;
        break;

    case HMESG_STATUS_FAIL:
        break; // No need to send state.

    default:
        return -1;
    }
    return total;
}

int unpack_state(hmesg_t* mesg, char* buf)
{
    int count = -1, total = 0;

    if (mesg->type == HMESG_SESSION ||
        mesg->type == HMESG_JOIN)
    {
        // Session state doesn't exist just yet.
        return 0;
    }

    switch (mesg->status) {
    case HMESG_STATUS_REQ:
        if (sscanf(buf, " state:%u %u%n", &mesg->unpacked_space.id,
                   &mesg->unpacked_best.id, &count) < 2)
            return -1;
        total += count;
        mesg->state.space = &mesg->unpacked_space;
        mesg->state.best  = &mesg->unpacked_best;

        count = scanstr_serial(&mesg->state.client, buf + total);
        if (count < 0) return -1;
        total += count;
        break;

    case HMESG_STATUS_OK:
    case HMESG_STATUS_BUSY:
        sscanf(buf, " state:%n", &count);
        if (count < 0) return -1;
        total += count;

        count = hspace_unpack(&mesg->unpacked_space, buf + total);
        if (count < 0) return -1;
        total += count;
        mesg->state.space = &mesg->unpacked_space;

        count = hpoint_unpack(&mesg->unpacked_best, buf + total);
        if (count < 0) return -1;
        total += count;
        mesg->state.best = &mesg->unpacked_best;
        break;

    case HMESG_STATUS_FAIL:
        break; // No need to send state.

    default:
        return -1;
    }
    return total;
}

int pack_data(char** buf, int* buflen, const hmesg_t* mesg)
{
    int count, total = 0;

    // Pack message data based on message type and status.
    switch (mesg->type) {
    case HMESG_SESSION:
        if (mesg->status == HMESG_STATUS_REQ) {
            count = hspace_pack(buf, buflen, mesg->state.space);
            if (count < 0) return -1;
            total += count;

            count = hcfg_pack(buf, buflen, mesg->data.cfg);
            if (count < 0) return -1;
            total += count;
        }
        break;

    case HMESG_JOIN:
        if (mesg->status == HMESG_STATUS_REQ) {
            count = printstr_serial(buf, buflen, mesg->data.string);
            if (count < 0) return -1;
            total += count;
            break;
        }
        else if (mesg->status == HMESG_STATUS_OK) {
            count = hspace_pack(buf, buflen, mesg->state.space);
            if (count < 0) return -1;
            total += count;
        }

    case HMESG_GETCFG:
    case HMESG_SETCFG:
    case HMESG_COMMAND:
        count = printstr_serial(buf, buflen, mesg->data.string);
        if (count < 0) return -1;
        total += count;
        break;

    case HMESG_FETCH:
        if (mesg->status == HMESG_STATUS_OK) {
            count = hpoint_pack(buf, buflen, mesg->data.point);
            if (count < 0) return -1;
            total += count;
        }
        break;

    case HMESG_REPORT:
        if (mesg->status == HMESG_STATUS_REQ) {
            count = snprintf_serial(buf, buflen, " %d", mesg->data.point->id);
            if (count < 0) return -1;
            total += count;

            count = hperf_pack(buf, buflen, mesg->data.perf);
            if (count < 0) return -1;
            total += count;
        }
        break;

    case HMESG_UNKNOWN:
    case HMESG_BEST:
        break;

    default:
        return -1;
    }
    return total;
}

int unpack_data(hmesg_t* mesg, char* buf)
{
    int count, total = 0;

    switch (mesg->type) {
    case HMESG_SESSION:
        if (mesg->status == HMESG_STATUS_REQ) {
            count = hspace_unpack(&mesg->unpacked_space, buf + total);
            if (count < 0) return -1;
            total += count;
            mesg->state.space = &mesg->unpacked_space;

            count = hcfg_unpack(&mesg->unpacked_cfg, buf + total);
            if (count < 0) return -1;
            total += count;
            mesg->data.cfg = &mesg->unpacked_cfg;
        }
        break;

    case HMESG_JOIN:
        if (mesg->status == HMESG_STATUS_REQ) {
            count = scanstr_serial(&mesg->data.string, buf + total);
            if (count < 0) return -1;
            total += count;
        }
        else if (mesg->status == HMESG_STATUS_OK) {
            count = hspace_unpack(&mesg->unpacked_space, buf + total);
            if (count < 0) return -1;
            total += count;
            mesg->state.space = &mesg->unpacked_space;
        }
        break;

    case HMESG_GETCFG:
    case HMESG_SETCFG:
    case HMESG_COMMAND:
        count = scanstr_serial(&mesg->data.string, buf + total);
        if (count < 0) return -1;
        total += count;
        break;

    case HMESG_FETCH:
        if (mesg->status == HMESG_STATUS_OK) {
            count = hpoint_unpack(&mesg->unpacked_point, buf + total);
            if (count < 0) return -1;
            total += count;
            mesg->data.point = &mesg->unpacked_point;
        }
        break;

    case HMESG_REPORT:
        if (mesg->status == HMESG_STATUS_REQ) {
            if (sscanf(buf + total, " %d%n",
                       &mesg->unpacked_point.id, &count) < 1)
                return -1;
            total += count;
            mesg->data.point = &mesg->unpacked_point;

            count = hperf_unpack(&mesg->unpacked_perf, buf + total);
            if (count < 0) return -1;
            total += count;
            mesg->data.perf = &mesg->unpacked_perf;
        }
        break;

    case HMESG_BEST:
        break;

    default:
        return -1;
    }
    return total;
}
