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

#ifndef __HRANGE_H__
#define __HRANGE_H__

#include "hval.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Integer domain range type.
 */
typedef struct range_int {
    long min;
    long max;
    long step;
} range_int_t;

/*
 * Real domain range type.
 */
typedef struct range_real {
    double min;
    double max;
    double step;
} range_real_t;

/*
 * Enumerated domain range type.
 */
typedef struct range_enum {
    char** set;
    int    len;
    int    cap;
} range_enum_t;

int range_enum_add_value(range_enum_t* bounds, char* str, const char** errptr);

/*
 * Range type to represent a single dimension of the tuning search space.
 */
typedef struct hrange {
    char* name;
    hval_type_t type;
    union {
        range_int_t  i;
        range_real_t r;
        range_enum_t e;
    } bounds;
} hrange_t;
#define HRANGE_INITIALIZER {0}
extern const hrange_t hrange_zero;

/*
 * Base structure management interface.
 */
int  hrange_copy(hrange_t* dst, const hrange_t* src);
void hrange_fini(hrange_t* range);
void hrange_scrub(hrange_t* range);

/*
 * Range value query interface.
 */
int           hrange_finite(const hrange_t* range);
unsigned long hrange_index(const hrange_t* range, const hval_t* val);
unsigned long hrange_limit(const hrange_t* range);
hval_t        hrange_random(const hrange_t* range, double entropy);
hval_t        hrange_value(const hrange_t* range, unsigned long idx);

/*
 * Data transmission interface.
 */
int hrange_pack(char** buf, int* buflen, const hrange_t* range);
int hrange_unpack(hrange_t* range, char* buf);
int hrange_parse(hrange_t* range, const char* buf, const char** errptr);

#ifdef __cplusplus
}
#endif

#endif // __HRANGE_H__
