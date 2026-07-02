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

#ifndef __HPOINT_H__
#define __HPOINT_H__

#include "hspace.h"
#include "hval.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Harmony structure that represents a point within a search space.
 */
typedef struct hpoint {
    unsigned id;
    hval_t*  term;
    int      len;
    int      cap;
} hpoint_t;

#define HPOINT_INITIALIZER {0}
extern const hpoint_t hpoint_zero;

/*
 * Base structure management interface.
 */
int  hpoint_init(hpoint_t* point, int newcap);
int  hpoint_copy(hpoint_t* dst, const hpoint_t* src);
void hpoint_fini(hpoint_t* point);
void hpoint_scrub(hpoint_t* point);

/*
 * Point comparison interface.
 */
int hpoint_align(hpoint_t* point, const hspace_t* space);
int hpoint_eq(const hpoint_t* a, const hpoint_t* b);
int hpoint_cmp(const hpoint_t* a, const hpoint_t* b);

/*
 * Data transmission interface.
 */
int hpoint_pack(char** buf, int* buflen, const hpoint_t* point);
int hpoint_unpack(hpoint_t* point, char* buf);
int hpoint_parse(hpoint_t* point, const char* buf, const hspace_t* space);

#ifdef __cplusplus
}
#endif

#endif
