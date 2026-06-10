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

#ifndef __HPERF_H__
#define __HPERF_H__

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Harmony structure that represents a (possibly multi-objective)
 * performance value.
 */
typedef struct hperf {
    double* obj;
    int     len;
    int     cap;
} hperf_t;
#define HPERF_INITIALIZER {0}
extern const hperf_t hperf_zero;

/*
 * Base structure management interface.
 */
int    hperf_init(hperf_t* perf, int newcap);
void   hperf_reset(hperf_t* perf);
int    hperf_copy(hperf_t* src, const hperf_t* dst);
void   hperf_fini(hperf_t* perf);

/*
 * Performance utility interface.
 */
int    hperf_cmp(const hperf_t* a, const hperf_t* b);
double hperf_unify(const hperf_t* perf);

/*
 * Data transmission interface.
 */
int    hperf_pack(char** buf, int* buflen, const hperf_t* perf);
int    hperf_unpack(hperf_t* perf, char* buf);

#ifdef __cplusplus
}
#endif

#endif
