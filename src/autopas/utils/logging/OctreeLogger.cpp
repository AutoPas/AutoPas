/**
 * @file OctreeLogger.cpp
 * @author Johannes Spies
 * @date 21.04.2021
 */

#include <spdlog/spdlog.h>
#include "autopas/utils/logging/OctreeLogger.h"
#include "autopas/utils/Timer.h"

autopas::OctreeLogger::OctreeLogger() {
//#ifdef AUTOPAS_LOG_OCTREE
//#endif
}

autopas::OctreeLogger::~OctreeLogger() {
//#ifdef AUTOPAS_LOG_OCTREE
    //spdlog::drop(_loggerName);
//#endif
}
