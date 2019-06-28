/**
 * @file MemoryProfiler.h
 * @author F. Gratl
 * @date 9/18/18
 */

#pragma once

#include <cstdlib>

namespace autopas::memoryProfiler {

/**
 * path to status file of current process
 */
constexpr const char *statusFileName = "/proc/self/status";

/**
 * Reads the current RAM (VmRSS) usage from the operating system.
 * @return Memory usage in kilobytes.
 */
size_t currentMemoryUsage();

}  // namespace autopas::memoryProfiler