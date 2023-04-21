/**
 * @file CompileInfo.cpp
 * @author F. Gratl
 * @date 20.04.2023
 */

#include "CompileInfo.h"

std::string utils::CompileInfo::getCompilerInfo() {
  std::stringstream ss;
#if defined(__GNUC__) and not(defined(__clang__) or defined(__INTEL_COMPILER))
  ss << "GCC ";
#endif
  // defined for:
  //   - icpc (without version)
  //   - icpx
  //   - gcc (without name)
  //   - clang
#if defined(__VERSION__)
  ss << __VERSION__;
#endif
  // icpc does not specify the version in __VERSION__
#if defined(__INTEL_COMPILER)
  ss << " (icpc version: " << __INTEL_COMPILER << "." << __INTEL_COMPILER_UPDATE << ")";
#endif
  // Base case if nothing could be found
  if (ss.str().empty()) {
    ss << "Unknown";
  }
  return ss.str();
}
