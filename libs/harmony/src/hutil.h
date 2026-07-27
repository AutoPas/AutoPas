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

#ifndef __HUTIL_H__
#define __HUTIL_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

int   file_exists(const char* filename);
void* file_map(const char* filename, size_t* size);
void  file_unmap(void* buf, size_t size);
int   valid_id(const char* key, int len);
int   file_read_line(FILE* fp, char** buf, int* cap,
                     char** line, char** end, const char** errstr);
char* search_path(const char* filename);
int   array_grow(void* buf, int* cap, int size);
char* stralloc(const char* in);
char* sprintf_alloc(const char* fmt, ...);
int   snprintf_grow(char** buf, int* buflen, const char* fmt, ...);
int   snprintf_serial(char** buf, int* buflen, const char* fmt, ...);
int   printstr_serial(char** buf, int* buflen, const char* str);
int   scanstr_serial(const char** str, char* buf);
int   unquote_string(const char* buf, char** token, const char** errptr);

#ifdef __cplusplus
}
#endif

#endif /* __HUTIL_H__ */
