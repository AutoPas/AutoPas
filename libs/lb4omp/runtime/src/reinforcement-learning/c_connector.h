#ifndef AAA_C_CONNECTOR_H
#define AAA_C_CONNECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "kmp_loopdata.h"

int search(const char *loop_id, int agent_type, LoopData *stats, int dimension);

void setEnvVar(const char* var_name, int value);

void setLoopData(const char *loop_id, LoopData *loop_data);

LoopData *getLoopData(const char *loop_id);

#ifdef __cplusplus
}
#endif


#endif