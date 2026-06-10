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

#ifdef __HLAYER_H__
#error "Multiple inclusion of hlayer.h detected."
#endif
#define __HLAYER_H__

#include "session-core.h"
#include "hplugin.h"
#include "hspace.h"
#include "hpoint.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Plug-in configuration variables read by Harmony plug-in loader.
 */
const        char        hplugin_type[] = "layer";
extern const char        hplugin_name[];
extern const hcfg_info_t hplugin_keyinfo[];

/*
 * Commented prototypes for layer plug-in functions.
 *
 * Actual function prototypes from this header file are not possible
 * because the hplugin_name should be different for each layer.  The
 * following section merely gives the prospective plug-in developer an
 * idea of what could be defined.
 */

/*
hplugin_data_t* <name>_alloc(void);

int <name>_init(hplugin_data_t* data, hspace_t* space);
int <name>_join(hplugin_data_t* data, const char* client);
int <name>_setcfg(hplugin_data_t* data, const char* key, const char* val);
int <name>_generate(hplugin_data_t* data, hflow_t* flow, htrial_t* trial);
int <name>_analyze(hplugin_data_t* data, hflow_t* flow, htrial_t* trial);
int <name>_fini(hplugin_data_t* data);
*/

#ifdef __cplusplus
}
#endif
