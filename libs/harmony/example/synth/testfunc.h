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

#ifndef __TESTFUNC_H__
#define __TESTFUNC_H__

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef double benchfunc_t(int n, double x[], double option[]);

typedef enum ftype {
    FTYPE_UNKNOWN,
    FTYPE_INT,
    FTYPE_REAL,

    FTYPE_MAX
} ftype_t;

typedef struct finfo {
    const char* name;
    char* title;
    int n_max;
    ftype_t type;
    double b_min;
    double b_max;
    int has_optimal;
    double optimal;
    benchfunc_t* f;
    char* description;
} finfo_t;

extern finfo_t flist[];

void flist_print(FILE* fd, int verbose);
finfo_t* flist_find(const char* name);

#ifdef __cplusplus
}
#endif

#endif
