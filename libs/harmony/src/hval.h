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

#ifndef __HVAL_H__
#define __HVAL_H__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Harmony structures that encapsulate values within a search space.
 */
typedef enum hval_type {
    HVAL_UNKNOWN = 0,
    HVAL_INT,  // Integer domain value.
    HVAL_REAL, // Real domain value.
    HVAL_STR,  // String domain value (for enumerated types).

    HVAL_MAX
} hval_type_t;

typedef union hval_value {
    long        i;
    double      r;
    const char* s;
} hval_value_t;

typedef struct hval {
    hval_type_t  type;
    hval_value_t value;
    char*        buf;
} hval_t;

#define HVAL_INITIALIZER {HVAL_UNKNOWN}
extern const hval_t hval_zero;

/*
 * Base structure management interface.
 */
int  hval_copy(hval_t* dst, const hval_t* src);
void hval_fini(hval_t* v);

/*
 * Value comparison interface.
 */
int hval_eq(const hval_t* a, const hval_t* b);

/*
 * Data transmission interface.
 */
int hval_pack(char** buf, int* buflen, const hval_t* v);
int hval_unpack(hval_t* v, char* buf);
int hval_parse(hval_t* v, hval_type_t type, const char* buf);

#ifdef __cplusplus
}
#endif

#endif
