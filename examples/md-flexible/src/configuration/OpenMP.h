/**
 * @file OpenMP.h
 * @date 25.09.2026
 * @author R. Horn
 */

#pragma once

#include "autopas/utils/WrapOpenMP.h"

#ifdef MD_FLEXIBLE_USE_TUNED_THREADS
#define MD_FLEXIBLE_NUM_THREADS num_threads(autopas::autopas_get_tuned_num_threads())
#else
// Do not specify number of threads in md-flexible's OpenMP macros and use the maximum (default) instead.
// This does not impact internal AutoPas thread count tuning.
#define MD_FLEXIBLE_NUM_THREADS
#endif
