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

/**
 * \file hclient.h
 * \brief Harmony client application function header.
 *
 * All clients must include this file to participate in a Harmony
 * tuning session.
 */

#ifndef __HCLIENT_H__
#define __HCLIENT_H__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DOXYGEN_SKIP

typedef struct hdesc hdesc_t;
typedef struct hdef  hdef_t;
typedef struct htask htask_t;

#endif

/*
 * Harmony descriptor interface.
 */
hdesc_t* ah_alloc(void);
int      ah_args(hdesc_t* hdesc, int* argc, char** argv);
int      ah_id(hdesc_t* hdesc, const char* id);
int      ah_connect(hdesc_t* hdesc, const char* host, int port);
int      ah_close(hdesc_t* hdesc);
void     ah_free(hdesc_t* hdesc);

/*
 * Search definition interface.
 */
hdef_t* ah_def_alloc(void);
hdef_t* ah_def_load(const char* filename);
int     ah_def_name(hdef_t* hdef, const char* name);
int     ah_def_int(hdef_t* hdef, const char* name,
                   long min, long max, long step, long* ptr);
int     ah_def_real(hdef_t* hdef, const char* name,
                    double min, double max, double step, double* ptr);
int     ah_def_enum(hdef_t* hdef, const char* name, const char** ptr);
int     ah_def_enum_value(hdef_t* hdef, const char* name, const char* value);
int     ah_def_strategy(hdef_t* hdef, const char* strategy);
int     ah_def_layers(hdef_t* hdef, const char* list);
int     ah_def_cfg(hdef_t* hdef, const char* key, const char* val);
void    ah_def_free(hdef_t* hdef);

/*
 * Search task control interface.
 */
htask_t* ah_start(hdesc_t* hdesc, hdef_t* hdef);
htask_t* ah_join(hdesc_t* hdesc, const char* name);
int      ah_pause(htask_t* htask);
int      ah_resume(htask_t* htask);
int      ah_restart(htask_t* htask);
int      ah_leave(htask_t* htask);
int      ah_kill(htask_t* htask);

/*
 * Search task interaction interface.
 */
int         ah_bind_int(htask_t* htask, const char* name, long* ptr);
int         ah_bind_real(htask_t* htask, const char* name, double* ptr);
int         ah_bind_enum(htask_t* htask, const char* name, const char** ptr);
long        ah_get_int(htask_t* htask, const char* name);
double      ah_get_real(htask_t* htask, const char* name);
const char* ah_get_enum(htask_t* htask, const char* name);
const char* ah_get_cfg(htask_t* htask, const char* key);
const char* ah_set_cfg(htask_t* htask, const char* key, const char* val);
int         ah_fetch(htask_t* htask);
int         ah_report(htask_t* htask, double* perf);
int         ah_report_one(htask_t* htask, int index, double value);
int         ah_best(htask_t* htask);
int         ah_converged(htask_t* htask);

/*
 * Harmony error reporting interface.
 */
const char* ah_error(void);
void        ah_error_clear(void);

/*
 * Deprecated API.
 *
 * These functions are slated for removal in a future release.
 */
#ifndef DOXYGEN_SKIP

#if defined(__INTEL_COMPILER)
#define DEPRECATED(message) __attribute__((__deprecated__))
#elif (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5)
#define DEPRECATED(message) __attribute__((__deprecated__(message)))
#elif defined(__GNUC__)
#define DEPRECATED(message) __attribute__((__deprecated__))
#else
#define DEPRECATED(message)
#endif

DEPRECATED("Use ah_alloc() instead")
hdesc_t* harmony_init(int* argc, char*** argv);

DEPRECATED("Use ah_args() instead")
int harmony_parse_args(hdesc_t* hdesc, int argc, char** argv);

DEPRECATED("Use ah_free() instead")
void harmony_fini(hdesc_t* hdesc);

DEPRECATED("Use ah_def_int() instead")
int harmony_int(hdesc_t* hdesc, const char* name,
                long min, long max, long step);

DEPRECATED("Use ah_def_real() instead")
int harmony_real(hdesc_t* hdesc, const char* name,
                 double min, double max, double step);

DEPRECATED("Use ah_def_enum() instead")
int harmony_enum(hdesc_t* hdesc, const char* name, const char* value);

DEPRECATED("Use ah_def_name() instead")
int harmony_session_name(hdesc_t* hdesc, const char* name);

DEPRECATED("Use ah_def_strategy() instead")
int harmony_strategy(hdesc_t* hdesc, const char* strategy);

DEPRECATED("Use ah_def_layers() instead")
int harmony_layers(hdesc_t* hdesc, const char* list);

DEPRECATED("Use ah_join() or ah_start() instead")
int harmony_launch(hdesc_t* hdesc, const char* host, int port);

DEPRECATED("Use ah_id() instead")
int harmony_id(hdesc_t* hdesc, const char* id);

DEPRECATED("Use ah_bind_int() instead")
int harmony_bind_int(hdesc_t* hdesc, const char* name, long* ptr);

DEPRECATED("Use ah_bind_real() instead")
int harmony_bind_real(hdesc_t* hdesc, const char* name, double* ptr);

DEPRECATED("Use ah_bind_enum() instead")
int harmony_bind_enum(hdesc_t* hdesc, const char* name, const char** ptr);

DEPRECATED("Use ah_join() instead")
int harmony_join(hdesc_t* hdesc, const char* host, int port, const char* name);

DEPRECATED("Use ah_detach() instead")
int harmony_leave(hdesc_t* hdesc);

DEPRECATED("Use ah_get_int() instead")
long harmony_get_int(hdesc_t* hdesc, const char* name);

DEPRECATED("Use ah_get_real() instead")
double harmony_get_real(hdesc_t* hdesc, const char* name);

DEPRECATED("Use ah_get_enum() instead")
const char* harmony_get_enum(hdesc_t* hdesc, const char* name);

DEPRECATED("Use ah_get_cfg() instead")
char* harmony_getcfg(hdesc_t* hdesc, const char* key);

DEPRECATED("Use ah_set_cfg() instead")
char* harmony_setcfg(hdesc_t* hdesc, const char* key, const char* val);

DEPRECATED("Use ah_fetch() instead")
int harmony_fetch(hdesc_t* hdesc);

DEPRECATED("Use ah_report() instead")
int harmony_report(hdesc_t* hdesc, double perf);

DEPRECATED("Use ah_report_one() instead")
int harmony_report_one(hdesc_t* hdesc, int index, double value);

DEPRECATED("Use ah_best() instead")
int harmony_best(hdesc_t* hdesc);

DEPRECATED("Use ah_converged() instead")
int harmony_converged(hdesc_t* hdesc);

DEPRECATED("Use ah_error() instead")
const char* harmony_error_string(hdesc_t* hdesc);

DEPRECATED("Use ah_error_clear() instead")
void harmony_error_clear(hdesc_t* hdesc);

#endif

#ifdef __cplusplus
}
#endif

#endif /* __HCLIENT_H__ */
