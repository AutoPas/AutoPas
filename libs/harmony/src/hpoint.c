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

#include "hpoint.h"
#include "hspace.h"
#include "hrange.h"
#include "hutil.h"

#include <stdlib.h> // For realloc().
#include <string.h> // For memcpy() and memcmp().
#include <assert.h> // For assert().
#include <ctype.h>  // For isspace().

const hpoint_t hpoint_zero = HPOINT_INITIALIZER;

/*
 * Internal helper function prototypes.
 */
static int align_int(hval_t* val, const hrange_t* range);
static int align_real(hval_t* val, const hrange_t* range);
static int align_str(hval_t* val, const hrange_t* range);

/*
 * Base structure management implementation.
 */
int hpoint_init(hpoint_t* point, int newcap)
{
    if (point->cap < newcap) {
        hval_t* newbuf = realloc(point->term, newcap * sizeof(*point->term));
        if (!newbuf)
            return -1;

        // Initialize any newly created hval_t structures.
        memset(newbuf + point->cap, 0,
               (newcap - point->cap) * sizeof(*point->term));

        point->term = newbuf;
        point->cap  = newcap;
    }
    return 0;
}

int hpoint_copy(hpoint_t* dst, const hpoint_t* src)
{
    if (dst->cap < src->len) {
        if (hpoint_init(dst, src->len) != 0)
            return -1;
    }

    for (int i = 0; i < src->len; ++i) {
        if (hval_copy(&dst->term[i], &src->term[i]) != 0)
            return -1;
    }
    dst->len = src->len;
    dst->id  = src->id;

    return 0;
}

void hpoint_fini(hpoint_t* point)
{
    for (int i = 0; i < point->cap; ++i)
        hval_fini(&point->term[i]);
    free(point->term);
}

void hpoint_scrub(hpoint_t* point)
{
    free(point->term);
}

/*
 * Point manipulation implementation.
 */
int hpoint_align(hpoint_t* point, const hspace_t* space)
{
    if (!point->id)
        return 0;

    for (int i = 0; i < space->len; ++i) {
        int retval;
        hval_t* val = &point->term[i];

        switch (space->dim[i].type) {
        case HVAL_INT:  retval = align_int(val, &space->dim[i]); break;
        case HVAL_REAL: retval = align_real(val, &space->dim[i]); break;
        case HVAL_STR:  retval = align_str(val, &space->dim[i]); break;
        default:        retval = -1;
        }
        if (retval != 0)
            return -1;
    }
    return 0;
}

int hpoint_eq(const hpoint_t* a, const hpoint_t* b)
{
    if (a->len != b->len)
        return 0;
    else
        return memcmp(a->term, b->term, a->len * sizeof(*a->term)) == 0;
}

int hpoint_cmp(const hpoint_t* a, const hpoint_t* b)
{
    if (a->len != b->len)
        return a->len - b->len;
    else
        return memcmp(a->term, b->term, a->len * sizeof(*a->term));
}

/*
 * Data transmission implementation.
 */
int hpoint_pack(char** buf, int* buflen, const hpoint_t* point)
{
    int total = snprintf_serial(buf, buflen, " pt:%u", point->id);
    if (total < 0) return -1;

    if (point->id) {
        int count = snprintf_serial(buf, buflen, " %d", point->len);
        if (count < 0) return -1;
        total += count;

        for (int i = 0; i < point->len; ++i) {
            count = hval_pack(buf, buflen, &point->term[i]);
            if (count < 0) return -1;
            total += count;
        }
    }
    return total;
}

int hpoint_unpack(hpoint_t* point, char* buf)
{
    int total;
    if (sscanf(buf, " pt:%u%n", &point->id, &total) < 1)
        return -1;

    if (point->id) {
        int newlen, count;

        if (sscanf(buf + total, " %d%n", &newlen, &count) < 1)
            return -1;
        total += count;

        if (point->cap < newlen) {
            if (hpoint_init(point, newlen) != 0)
                return -1;
        }

        for (int i = 0; i < newlen; ++i) {
            count = hval_unpack(&point->term[i], buf + total);
            if (count < 0) return -1;
            total += count;
        }
        point->len = newlen;
    }
    return total;
}

int hpoint_parse(hpoint_t* point, const char* buf, const hspace_t* space)
{
    if (point->cap < space->len)
        hpoint_init(point, space->len);

    for (int i = 0; i < space->len; ++i) {
        hval_t* val = &point->term[i];
        int count = hval_parse(val, space->dim[i].type, buf);
        if (count < 0) return -1;
        buf += count;

        // Skip whitespace and a single separator character.
        while (isspace(*buf)) ++buf;
        if (*buf == ',' || *buf == ';')
            ++buf;
    }
    if (*buf != '\0')
        return -1;

    point->len = space->len;
    return 0;
}

/*
 * Internal helper function implementation.
 */
int align_int(hval_t* val, const hrange_t* range)
{
    if (val->value.i < range->bounds.i.min)
        val->value.i = range->bounds.i.min;
    if (val->value.i > range->bounds.i.max)
        val->value.i = range->bounds.i.max;

    long rank;
    rank  = val->value.i;
    rank += range->bounds.i.step / 2; // To support proper integer rounding.
    rank -= range->bounds.i.min;
    rank /= range->bounds.i.step;

    val->value.i  = range->bounds.i.step;
    val->value.i *= rank;
    val->value.i += range->bounds.i.min;

    return 0;
}

int align_real(hval_t* val, const hrange_t* range)
{
    if (val->value.r < range->bounds.r.min)
        val->value.r = range->bounds.r.min;
    if (val->value.r > range->bounds.r.max)
        val->value.r = range->bounds.r.max;

    if (range->bounds.r.step > 0.0) {
        double rank;
        rank  = val->value.r;
        rank -= range->bounds.r.min;
        rank /= range->bounds.r.step;
        rank += 0.5; // To support proper integer rounding.

        val->value.r  = range->bounds.r.step;
        val->value.r *= (long) rank;
        val->value.r += range->bounds.r.min;
    }
    return 0;
}

int align_str(hval_t* val, const hrange_t* range)
{
    // First pass compare by address.
    for (int i = 0; i < range->bounds.e.len; ++i)
        if (val->value.s == range->bounds.e.set[i])
            return 0;

    // Second pass compare by value.
    for (int i = 0; i < range->bounds.e.len; ++i) {
        if (strcmp(val->value.s, range->bounds.e.set[i]) == 0) {
            val->value.s = range->bounds.e.set[i];
            free(val->buf);
            val->buf = NULL;
            return 0;
        }
    }
    return -1;
}
