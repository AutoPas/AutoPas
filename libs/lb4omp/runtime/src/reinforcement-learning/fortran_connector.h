#ifndef FORTRAN_CONNECTOR_H
#define FORTRAN_CONNECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "kmp_loopdata.h"

int search(const char *loop_id, int agent_type, LoopData *stats, int dimension);

int my_cpp_function(int arg);

#ifdef __cplusplus
}
#endif

#endif