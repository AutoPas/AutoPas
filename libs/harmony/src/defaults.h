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
#ifndef __DEFAULTS_H__
#define __DEFAULTS_H__

#define SESSION_CORE_EXECFILE "session-core"
#define DEFAULT_HOST "localhost"
#define DEFAULT_PORT 1979

/*
 * The following definitions allow the compiler to find typos for us.
 * Plug-in developers need not augment this file until their strategy
 * or layer is incorporated into the Active Harmony codebase.
 */

// Session-wide configuration variables.
#define CFGKEY_HARMONY_HOME       "HARMONY_HOME"
#define CFGKEY_HARMONY_HOST       "HARMONY_HOST"
#define CFGKEY_HARMONY_PORT       "HARMONY_PORT"
#define CFGKEY_RANDOM_SEED        "RANDOM_SEED"
#define CFGKEY_PERF_COUNT         "PERF_COUNT"
#define CFGKEY_GEN_COUNT          "GEN_COUNT"
#define CFGKEY_CLIENT_COUNT       "CLIENT_COUNT"
#define CFGKEY_STRATEGY           "STRATEGY"
#define CFGKEY_LAYERS             "LAYERS"

// Search state configuration variables.
#define CFGKEY_PAUSED             "PAUSED"
#define CFGKEY_CONVERGED          "CONVERGED"
#define CFGKEY_CURRENT_CLIENT     "CURRENT_CLIENT"

// Plug-in strategy configuration variables.
#define CFGKEY_PASSES             "PASSES"
#define CFGKEY_INIT_METHOD        "INIT_METHOD"
#define CFGKEY_INIT_RADIUS        "INIT_RADIUS"
#define CFGKEY_INIT_POINT         "INIT_POINT"
#define CFGKEY_SIMPLEX_SIZE       "SIMPLEX_SIZE"
#define CFGKEY_REFLECT            "REFLECT"
#define CFGKEY_EXPAND             "EXPAND"
#define CFGKEY_CONTRACT           "CONTRACT"
#define CFGKEY_SHRINK             "SHRINK"
#define CFGKEY_FVAL_TOL           "FVAL_TOL"
#define CFGKEY_SIZE_TOL           "SIZE_TOL"
#define CFGKEY_DIST_TOL           "DIST_TOL"
#define CFGKEY_TOL_CNT            "TOL_CNT"
#define CFGKEY_REJECT_METHOD      "REJECT_METHOD"
#define CFGKEY_ANGEL_LOOSE        "ANGEL_LOOSE"
#define CFGKEY_ANGEL_MULT         "ANGEL_MULT"
#define CFGKEY_ANGEL_ANCHOR       "ANGEL_ANCHOR"
#define CFGKEY_ANGEL_SAMESIMPLEX  "ANGEL_SAMESIMPLEX"
#define CFGKEY_ANGEL_LEEWAY       "ANGEL_LEEWAY"
#define CFGKEY_ANGEL_PHASE        "ANGEL_PHASE"

// Plug-in layer configuration variables.
#define CFGKEY_AGG_FUNC           "AGG_FUNC"
#define CFGKEY_AGG_TIMES          "AGG_TIMES"
#define CFGKEY_CACHE_FILE         "CACHE_FILE"
#define CFGKEY_GROUP_LIST         "GROUP_LIST"
#define CFGKEY_LOG_FILE           "LOG_FILE"
#define CFGKEY_LOG_MODE           "LOG_MODE"
#define CFGKEY_OC_BIN             "OC_BIN"
#define CFGKEY_OC_CONSTRAINTS     "OC_CONSTRAINTS"
#define CFGKEY_OC_FILE            "OC_FILE"
#define CFGKEY_OC_QUIET           "OC_QUIET"
#define CFGKEY_TAUDB_NAME         "TAUDB_NAME"
#define CFGKEY_TAUDB_STORE_METHOD "TAUDB_STORE_METHOD"
#define CFGKEY_TAUDB_STORE_NUM    "TAUDB_STORE_NUM"
#define CFGKEY_XML_FILE           "XML_FILE"

// Codegen server configuration variables.
#define CFGKEY_SERVER_URL         "SERVER_URL"
#define CFGKEY_TARGET_URL         "TARGET_URL"
#define CFGKEY_REPLY_URL          "REPLY_URL"
#define CFGKEY_TMPDIR             "TMPDIR"
#define CFGKEY_SLAVE_LIST         "SLAVE_LIST"
#define CFGKEY_SLAVE_PATH         "SLAVE_PATH"

#endif
