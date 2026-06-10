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

#ifdef __HSTRATEGY_H__
#error "Multiple inclusion of hstrategy.h detected."
#endif
#define __HSTRATEGY_H__

#include "session-core.h"
#include "hplugin.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Plug-in configuration variables read by Harmony plug-in loader.
 */
const        char        hplugin_type[] = "strategy";
extern const char        hplugin_name[];
extern const hcfg_info_t hplugin_keyinfo[];

/*
 * The following functions are used for plug-in setup.
 *
 * They are called once along with each search task that uses this
 * plug-in.
 */
hplugin_data_t* strategy_alloc(void);

/*
 * The following functions are required.
 *
 * Active Harmony will not recognize shared objects as search
 * strategies unless these functions exist.
 */
int strategy_analyze(hplugin_data_t* data, htrial_t* trial);
int strategy_best(hplugin_data_t* data, hpoint_t* point);
int strategy_generate(hplugin_data_t* data, hflow_t* flow, hpoint_t* point);
int strategy_rejected(hplugin_data_t* data, hflow_t* flow, hpoint_t* point);

/*
 * The following functions are optional.
 *
 * They will be invoked at the appropriate time if and only if they
 * exist.
 */
int strategy_init(hplugin_data_t* data, hspace_t* space);
int strategy_join(hplugin_data_t* data, const char* client);
int strategy_setcfg(hplugin_data_t* data, const char* key, const char* val);
int strategy_fini(hplugin_data_t* data);

#ifdef __cplusplus
}
#endif
