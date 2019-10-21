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

#include "hval.h"
#include "hutil.h"

#include <ctype.h> // For isspace().

const hval_t hval_zero = HVAL_INITIALIZER;

/*
 * Internal helper function prototypes.
 */
static int parse_int(hval_t* val, const char* buf);
static int parse_real(hval_t* val, const char* buf);
static int parse_str(hval_t* val, const char* buf);

/*
 * Base structure management implementation.
 */
int hval_copy(hval_t* dst, const hval_t* src)
{
    // Free heap data allocated by the destination structure.
    free(dst->buf);

    *dst = *src;
    if (src->buf) {
        dst->buf = stralloc(src->buf);
        dst->value.s = dst->buf;
        if (!dst->buf)
            return -1;
    }
    return 0;
}

void hval_fini(hval_t* val)
{
    free(val->buf);
}

/*
 * Value utility interface implementation.
 */
int hval_eq(const hval_t* a, const hval_t* b)
{
    if (a->type != b->type)
        return 0;

    switch (a->type) {
    case HVAL_INT:  return (a->value.i == b->value.i);
    case HVAL_REAL: return (a->value.r == b->value.r);
    case HVAL_STR:  return (a->value.s == b->value.s);
    default:        return 0;
    }
}

/*
 * Data transmission implementation.
 */
int hval_pack(char** buf, int* buflen, const hval_t* val)
{
    switch (val->type) {
    case HVAL_INT:
        return snprintf_serial(buf, buflen, " i%ld", val->value.i);
    case HVAL_REAL:
        return snprintf_serial(buf, buflen, " r%a(%g)",
                               val->value.r, val->value.r);
    case HVAL_STR:
        return printstr_serial(buf, buflen, val->value.s);
    default:
        return -1;
    }
}

int hval_unpack(hval_t* val, char* buf)
{
    int count = 0;

    while (isspace( buf[count] ))
        ++count;

    if (buf[count] == 'i') {
        val->type = HVAL_INT;
        if (sscanf(buf, " i%ld%n", &val->value.i, &count) < 1)
            return -1;
    }
    else if (buf[count] == 'r') {
        val->type = HVAL_REAL;
        if (sscanf(buf, " r%la(%*f)%n", &val->value.r, &count) < 1)
            return -1;
    }
    else {
        val->type = HVAL_STR;
        count = scanstr_serial(&val->value.s, buf);
    }
    return count;
}

int hval_parse(hval_t* val, hval_type_t type, const char* buf)
{
    int span;

    switch (type) {
    case HVAL_INT:  span = parse_int(val, buf);  break;
    case HVAL_REAL: span = parse_real(val, buf); break;
    case HVAL_STR:  span = parse_str(val, buf);  break;
    default:        return -1;
    }

    val->type = type;
    return span;
}

/*
 * Internal helper function implementation.
 */
int parse_int(hval_t* val, const char* buf)
{
    int span = -1;
    sscanf(buf, " %ld%n", &val->value.i, &span);
    return span;
}

int parse_real(hval_t* val, const char* buf)
{
    int span = -1;
    sscanf(buf, " %lf%n", &val->value.r, &span);
    return span;
}

int parse_str(hval_t* val, const char* buf)
{
    int span = unquote_string(buf, &val->buf, NULL);
    val->value.s = val->buf;
    return span;
}
