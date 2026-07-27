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

#include "hrange.h"
#include "hval.h"
#include "hmesg.h"
#include "hutil.h"

#include <string.h> // For strcmp().
#include <ctype.h>  // For isspace().
#include <math.h>   // For isnan() and NAN.
#include <limits.h> // For LONG_MIN.

const hrange_t hrange_zero = HRANGE_INITIALIZER;

/*
 * Internal helper function prototypes.
 */
static int copy_enum(range_enum_t* dst, const range_enum_t* src);
static int parse_int(range_int_t* bounds, const char* buf,
                     const char** errptr);
static int parse_real(range_real_t* bounds, const char* buf,
                      const char** errptr);
static int parse_enum(range_enum_t* bounds, const char* buf,
                      const char** errptr);
static unsigned long index_of_int(const range_int_t* bounds, long val);
static unsigned long index_of_real(const range_real_t* bounds, double val);
static unsigned long index_of_enum(const range_enum_t* bounds,
                                   const char* val);
static unsigned long limit_of_int(const range_int_t* bounds);
static unsigned long limit_of_real(const range_real_t* bounds);
static unsigned long limit_of_enum(const range_enum_t* bounds);
static hval_t random_int(const range_int_t* bounds, double entropy);
static hval_t random_real(const range_real_t* bounds, double entropy);
static hval_t random_enum(const range_enum_t* bounds, double entropy);
static hval_t value_of_int(const range_int_t* bounds, unsigned long idx);
static hval_t value_of_real(const range_real_t* bounds, unsigned long idx);
static hval_t value_of_enum(const range_enum_t* bounds, unsigned long idx);

/*
 * Enumerated domain range structure utility interface implementation.
 */
int range_enum_add_value(range_enum_t* bounds, char* str, const char** errptr)
{
    for (int i = 0; i < bounds->len; ++i) {
        if (strcmp(bounds->set[i], str) == 0) {
            *errptr = "Cannot add duplicate value to enumerated domain";
            return -1;
        }
    }

    if (bounds->len == bounds->cap) {
        if (array_grow(&bounds->set, &bounds->cap,
                       sizeof(*bounds->set)) != 0)
        {
            *errptr = "Could not extend enumerated domain's value list";
            return -1;
        }
    }

    bounds->set[ bounds->len++ ] = str;
    return 0;
}

/*
 * Base structure management implementation.
 */
int hrange_copy(hrange_t* dst, const hrange_t* src)
{
    // Free heap data allocated by the destination structure.
    hrange_fini(dst);

    dst->name = stralloc(src->name);
    if (!dst->name)
        return -1;

    dst->type = src->type;
    switch (src->type) {
    case HVAL_INT:  dst->bounds.i = src->bounds.i; break;
    case HVAL_REAL: dst->bounds.r = src->bounds.r; break;
    case HVAL_STR:  if (copy_enum(&dst->bounds.e, &src->bounds.e) == 0) break;
    default:        return -1;
    }
    return 0;
}

void hrange_fini(hrange_t* range)
{
    if (range->type == HVAL_STR) {
        for (int i = 0; i < range->bounds.e.len; ++i)
            free(range->bounds.e.set[i]);
        free(range->bounds.e.set);
    }
    free(range->name);
}

void hrange_scrub(hrange_t* range)
{
    if (range->type == HVAL_STR)
        free(range->bounds.e.set);
}

/*
 * Range value query implementation.
 */
int hrange_finite(const hrange_t* range)
{
    return (range->type != HVAL_REAL || range->bounds.r.step > 0.0);
}

unsigned long hrange_index(const hrange_t* range, const hval_t* val)
{
    switch (range->type) {
    case HVAL_INT:  return index_of_int(&range->bounds.i, val->value.i);
    case HVAL_REAL: return index_of_real(&range->bounds.r, val->value.r);
    case HVAL_STR:  return index_of_enum(&range->bounds.e, val->value.s);
    default:        return 0;
    }
}

unsigned long hrange_limit(const hrange_t* range)
{
    switch (range->type) {
    case HVAL_INT:  return limit_of_int(&range->bounds.i);
    case HVAL_REAL: return limit_of_real(&range->bounds.r);
    case HVAL_STR:  return limit_of_enum(&range->bounds.e);
    default:        return 0;
    }
}

hval_t hrange_random(const hrange_t* range, double entropy)
{
    switch (range->type) {
    case HVAL_INT:  return random_int(&range->bounds.i, entropy);
    case HVAL_REAL: return random_real(&range->bounds.r, entropy);
    case HVAL_STR:  return random_enum(&range->bounds.e, entropy);
    default:        return hval_zero;
    }
}

hval_t hrange_value(const hrange_t* range, unsigned long idx)
{
    switch (range->type) {
    case HVAL_INT:  return value_of_int(&range->bounds.i, idx);
    case HVAL_REAL: return value_of_real(&range->bounds.r, idx);
    case HVAL_STR:  return value_of_enum(&range->bounds.e, idx);
    default:        return hval_zero;
    }
}

/*
 * Data transmission implementation.
 */
int hrange_pack(char** buf, int* buflen, const hrange_t* range)
{
    int i, count, total;
    const char* type_str;

    count = snprintf_serial(buf, buflen, " range:");
    if (count < 0) return -1;
    total = count;

    count = printstr_serial(buf, buflen, range->name);
    if (count < 0) return -1;
    total += count;

    switch (range->type) {
    case HVAL_INT:  type_str = "INT"; break;
    case HVAL_REAL: type_str = "REAL"; break;
    case HVAL_STR:  type_str = "ENUM"; break;
    default: return -1;
    }

    count = snprintf_serial(buf, buflen, " %s", type_str);
    if (count < 0) return -1;
    total += count;

    switch (range->type) {
    case HVAL_INT:
        count = snprintf_serial(buf, buflen, " %ld %ld %ld",
                                range->bounds.i.min,
                                range->bounds.i.max,
                                range->bounds.i.step);
        if (count < 0) return -1;
        total += count;
        break;

    case HVAL_REAL:
        count = snprintf_serial(buf, buflen, " %la %la %la",
                                range->bounds.r.min,
                                range->bounds.r.max,
                                range->bounds.r.step);
        if (count < 0) return -1;
        total += count;
        break;

    case HVAL_STR:
        count = snprintf_serial(buf, buflen, " %d", range->bounds.e.len);
        if (count < 0) return -1;
        total += count;

        for (i = 0; i < range->bounds.e.len; ++i) {
            count = printstr_serial(buf, buflen, range->bounds.e.set[i]);
            if (count < 0) return -1;
            total += count;
        }
        break;

    default:
        return -1;
    }
    return total;
}

int hrange_unpack(hrange_t* range, char* buf)
{
    int count = -1, total;

    sscanf(buf, " range:%n", &count);
    if (count < 0)
        return -1;
    total = count;

    // Free heap data currently allocated for the structure.
    hrange_scrub(range);

    count = scanstr_serial((const char**)&range->name, buf + total);
    if (count < 0) return -1;
    total += count;

    char type_str[5];
    if (sscanf(buf + total, " %4s%n", type_str, &count) < 1)
        return -1;
    total += count;

    if      (strcmp(type_str, "INT") == 0) range->type = HVAL_INT;
    else if (strcmp(type_str, "REAL") == 0) range->type = HVAL_REAL;
    else if (strcmp(type_str, "ENUM") == 0) range->type = HVAL_STR;
    else return -1;

    switch (range->type) {
    case HVAL_INT:
        if (sscanf(buf + total, " %ld %ld %ld%n",
                   &range->bounds.i.min,
                   &range->bounds.i.max,
                   &range->bounds.i.step,
                   &count) < 3)
            return -1;
        total += count;
        break;

    case HVAL_REAL:
        if (sscanf(buf + total, " %la %la %la%n",
                   &range->bounds.r.min,
                   &range->bounds.r.max,
                   &range->bounds.r.step,
                   &count) < 3)
            return -1;
        total += count;
        break;

    case HVAL_STR:
        if (sscanf(buf + total, " %d%n", &range->bounds.e.len, &count) > 0) {
            total += count;

            char** newbuf = malloc(range->bounds.e.len * sizeof(*newbuf));
            if (!newbuf)
                return -1;

            for (int i = 0; i < range->bounds.e.len; ++i) {
                count = scanstr_serial((const char**) &newbuf[i], buf + total);
                if (count < 0) return -1;
                total += count;
            }

            range->bounds.e.set = newbuf;
            range->bounds.e.cap = range->bounds.e.len;
            break;
        }

    default:
        return -1;
    }
    return total;
}

int hrange_parse(hrange_t* range, const char* buf, const char** errptr)
{
    int id, idlen, bounds = 0, tail = 0;

    while (isspace(*buf)) ++buf;
    if (*buf == '\0')
        return 0;

    sscanf(buf, " int %n%*[^=]%n=%n", &id, &idlen, &bounds);
    if (bounds) {
        int len = parse_int(&range->bounds.i, buf + bounds, errptr);
        if (len == -1)
            goto error;

        range->type = HVAL_INT;
        tail += len;
        goto found;
    }

    sscanf(buf, " real %n%*[^=]%n=%n", &id, &idlen, &bounds);
    if (bounds) {
        int len = parse_real(&range->bounds.r, buf + bounds, errptr);
        if (len == -1)
            goto error;

        range->type = HVAL_REAL;
        tail += len;
        goto found;
    }

    sscanf(buf, " enum %n%*[^=]%n=%n", &id, &idlen, &bounds);
    if (bounds) {
        int len = parse_enum(&range->bounds.e, buf + bounds, errptr);
        if (len == -1)
            goto error;

        range->type = HVAL_STR;
        tail += len;
        goto found;
    }
    *errptr = "Unknown tuning variable type";
    goto error;

  found:
    while (isspace(buf[idlen - 1])) --idlen;
    if (!valid_id(&buf[id], idlen - id)) {
        *errptr = "Invalid variable name";
        goto error;
    }

    if (buf[bounds + tail] != '\0') {
        *errptr = "Invalid trailing characters";
        goto error;
    }

    range->name = sprintf_alloc("%.*s", idlen - id, &buf[id]);
    if (!range->name) {
        *errptr = "Cannot copy range name";
        goto error;
    }
    return 1;

  error:
    hrange_fini(range);
    return -1;
}

/*
 * Internal helper function implementation.
 */
int copy_enum(range_enum_t* dst, const range_enum_t* src)
{
    dst->set = malloc(src->cap * sizeof(*dst->set));
    if (!dst->set)
        return -1;

    dst->cap = src->cap;
    for (dst->len = 0; dst->len < src->len; ++dst->len) {
        dst->set[ dst->len ] = stralloc(src->set[ dst->len ]);
        if (!dst->set[ dst->len ])
            return -1;
    }
    return 0;
}

/*
 * Range definition parsing implementation.
 */
int parse_int(range_int_t* bounds, const char* buf, const char** errptr)
{
    *bounds = (range_int_t){LONG_MIN, LONG_MIN, 1};

    int tail = 0;
    sscanf(buf, " min: %ld max: %ld %nstep: %ld %n",
           &bounds->min,
           &bounds->max, &tail,
           &bounds->step, &tail);

    if (!tail) {
        if (bounds->min == LONG_MIN)
            *errptr = "Invalid integer-domain minimum range value";
        else
            *errptr = "Invalid integer-domain maximum range value";
        return -1;
    }
    return tail;
}

int parse_real(range_real_t* bounds, const char* buf, const char** errptr)
{
    *bounds = (range_real_t){NAN, NAN, NAN};

    int tail = 0;
    sscanf(buf, " min: %lf max: %lf step: %lf %n",
           &bounds->min,
           &bounds->max,
           &bounds->step, &tail);

    if (!tail) {
        if (isnan(bounds->min))
            *errptr = "Invalid real-domain minimum range value";
        else if (isnan(bounds->max))
            *errptr = "Invalid real-domain minimum range value";
        else
            *errptr = "Invalid real-domain step value";
        return -1;
    }
    return tail;
}

int parse_enum(range_enum_t* bounds, const char* buf, const char** errptr)
{
    *bounds = (range_enum_t){NULL, 0, 0};

    int tail = 0;
    while (buf[tail]) {
        char* token;
        int len = unquote_string(buf + tail, &token, errptr);
        if (len == -1)
            return -1;

        tail += len;
        if (token && range_enum_add_value(bounds, token, errptr) != 0) {
            free(token);
            return -1;
        }
        if (!token && buf[tail] != '\0') {
            *errptr = "Empty enumerated-domain value";
            return -1;
        }

        while (isspace(buf[tail])) ++tail;
        if (buf[tail] == ',') ++tail;
    }
    if (bounds->len == 0) {
        *errptr = "No enumerated string tokens found";
        return -1;
    }
    return tail;
}

/*
 * Value to range index implementation.
 */
unsigned long index_of_int(const range_int_t* bounds, long val)
{
    if (val < bounds->min) return 0;
    if (val > bounds->max) return limit_of_int(bounds) - 1;

    val += bounds->step / 2; // To support proper integer rounding.
    val -= bounds->min;
    val /= bounds->step;
    return (unsigned long) val;
}

unsigned long index_of_real(const range_real_t* bounds, double val)
{
    if (bounds->step > 0.0) {
        if (val < bounds->min) return 0;
        if (val > bounds->max) return limit_of_real(bounds) - 1;

        val -= bounds->min;
        val /= bounds->step;
        val += 0.5; // To support proper integer rounding.
        return (unsigned long) val;
    }
    return 0;
}

unsigned long index_of_enum(const range_enum_t* bounds, const char* val)
{
    int index;

    for (index = 0; index < bounds->len; ++index) {
        if (bounds->set[index] == val)
            break;
    }
    return (unsigned long) index;
}

/*
 * Count of range indexes implementation.
 */
unsigned long limit_of_int(const range_int_t* bounds)
{
    unsigned long bins;

    bins  = bounds->max;
    bins -= bounds->min;
    bins /= bounds->step;

    return 1UL + bins;
}

unsigned long limit_of_real(const range_real_t* bounds)
{
    if (bounds->step > 0.0) {
        double bins;

        bins  = bounds->max;
        bins -= bounds->min;
        bins /= bounds->step;

        return 1UL + (unsigned long) bins;
    }
    return 0;
}

unsigned long limit_of_enum(const range_enum_t* bounds)
{
    return (unsigned long) bounds->len;
}

/*
 * Random range value implementation.
 */
hval_t random_int(const range_int_t* bounds, double entropy)
{
    unsigned long idx = (unsigned long) (entropy * limit_of_int(bounds));
    return value_of_int(bounds, idx);
}

hval_t random_real(const range_real_t* bounds, double entropy)
{
    if (bounds->step > 0.0) {
        unsigned long idx = (unsigned long) (entropy * limit_of_real(bounds));
        return value_of_real(bounds, idx);
    }
    else {
        hval_t val;

        val.type     = HVAL_REAL;
        val.value.r  = bounds->max;
        val.value.r -= bounds->min;
        val.value.r *= entropy;
        val.value.r += bounds->min;
        val.buf      = NULL;

        return val;
    }
}

hval_t random_enum(const range_enum_t* bounds, double entropy)
{
    unsigned long idx = (unsigned long) (entropy * limit_of_enum(bounds));
    return value_of_enum(bounds, idx);
}

/*
 * Index to range value implementation.
 */
hval_t value_of_int(const range_int_t* bounds, unsigned long idx)
{
    hval_t val = hval_zero;

    val.type     = HVAL_INT;
    val.value.i  = bounds->step;
    val.value.i *= idx;
    val.value.i += bounds->min;
    val.buf      = NULL;

    return val;
}

hval_t value_of_real(const range_real_t* bounds, unsigned long idx)
{
    if (bounds->step > 0.0) {
        hval_t val = hval_zero;

        val.type     = HVAL_REAL;
        val.value.r  = bounds->step;
        val.value.r *= idx;
        val.value.r += bounds->min;
        val.buf      = NULL;

        return val;
    }
    return hval_zero;
}

hval_t value_of_enum(const range_enum_t* bounds, unsigned long idx)
{
    if (idx < bounds->len) {
        hval_t val = hval_zero;

        val.type    = HVAL_STR;
        val.value.s = bounds->set[ idx ];
        val.buf     = NULL;

        return val;
    }
    return hval_zero;
}
