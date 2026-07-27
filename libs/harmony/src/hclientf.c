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
#include "hclient.h"
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#if 0
int harmony_connect_(char* host, int* port, int* use_sigs,
                     int* relocated, int* hdesc)
{
    *hdesc = harmony_connect(host, *port, *use_sigs, *relocated);
    return (*hdesc >= 0 ? 0 : -1);
}

int harmony_disconnect_(int* hdesc)
{
    return harmony_disconnect(*hdesc);
}

int harmony_getcfg_(int* hdesc, const char* key, char* value, long len)
{
    int retval;
    char* newstr;
    retval = harmony_getcfg(*hdesc, key, &newstr);
    strncpy(value, newstr, len);
    free(newstr);

    return retval;
}

int harmony_application_setup_(int* hdesc, char* description)
{
    return harmony_application_setup(*hdesc, description);
}

int harmony_application_setup_file_(int* hdesc, char* fname)
{
    return harmony_application_setup_file(*hdesc, fname);
}

long harmony_add_variable_(int* hdesc, char* appName, char* bundleName,
                           int* type, int* local)
{
    return (long)harmony_add_variable(*hdesc, appName, bundleName,
                                      *type, *local);
}

int harmony_set_variable_int_(int* hdesc, long* var)
{
    return harmony_set_variable(*hdesc, (void*)(*var));
}

int harmony_set_all_(int* hdesc)
{
    return harmony_set_all(*hdesc);
}

void harmony_request_variable_int_(int* hdesc, char* name, int* var)
{
    *var = *(int*)harmony_request_variable(*hdesc, name);
}

void harmony_request_variable_double_(int* hdesc, char* name, double* var)
{
    *var = *(double*)harmony_request_variable(*hdesc, name);
}

int harmony_request_all_(int* hdesc, int* pull)
{
    return harmony_request_all(*hdesc, *pull);
}

int harmony_performance_update_int_(int* hdesc, int* value)
{
    return harmony_performance_update_int(*hdesc, *value);
}

int harmony_performance_update_double_(int* hdesc, double* value)
{
    return harmony_performance_update_double(*hdesc, *value);
}

int harmony_get_best_configuration_(int* hdesc, char* ret, long len)
{
    char* newstr = harmony_get_best_configuration(*hdesc);
    if (newstr) {
        strncpy(ret, newstr, len);
        free(newstr);
        return 0;
    }
    return -1;
}

void harmony_check_convergence_(int* hdesc, int* retval)
{
    *retval = harmony_check_convergence(*hdesc);
}

int harmony_code_generation_complete_(int* hdesc)
{
    return harmony_code_generation_complete(*hdesc);
}

long harmony_database_lookup_(int* hdesc)
{
    return (long)harmony_database_lookup(*hdesc);
}

int harmony_pseudo_barrier_(int* hdesc)
{
    return harmony_pseudo_barrier(*hdesc);
}

#endif

#ifdef __cplusplus
}
#endif
