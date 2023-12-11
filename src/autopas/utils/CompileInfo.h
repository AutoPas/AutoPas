/**
 * @file CompileInfo.h
 * @author F. Gratl
 * @date 20.04.2023
 */

#pragma once

#include <string>

namespace autopas::utils::CompileInfo {

/**
 * Get name and version number of a list of known compilers.
 * Currently known compilers:
 *   - clang
 *   - gcc
 *   - icpc
 *   - icpx
 *   - fcc (trad and clang mode)
 * @return String : "Name Version.Number"
 */
std::string getCompilerInfo();
}  // namespace autopas::utils::CompileInfo
