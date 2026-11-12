/**
 * @file OpenMP.h
 * @date 25.09.2026
 * @author R. Horn
 */

#pragma once

#include "autopas/utils/WrapOpenMP.h"

#ifdef MD_FLEXIBLE_USE_TUNED_THREADS
#define MD_FLEXIBLE_NUM_THREADS AUTOPAS_NUM_THREADS
#else
/**
 * Empty macro to disable using tuned number of threads.
 * This does not impact internal AutoPas thread count tuning.
 */
#define MD_FLEXIBLE_NUM_THREADS
#endif
