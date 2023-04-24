/**
 * @file CompileInfo.cpp
 * @author F. Gratl
 * @date 20.04.2023
 */

#include "CompileInfo.h"

std::string utils::CompileInfo::getCompilerInfo() {
  std::stringstream ss;
  // prefix if not included in VERSION
#if defined(__FUJITSU) or defined(__CLANG_FUJITSU)
  ss << "Fujitsu ";
#elif defined(__GNUC__) and not(defined(__clang__) or defined(__INTEL_COMPILER))
  ss << "GCC ";
#endif

  // Main version info. Defined for:
  //   - icpc (without version)
  //   - icpx
  //   - gcc (without name)
  //   - clang
  //   - fcc (weird version number)
#if defined(__VERSION__)
  ss << __VERSION__;
#endif

  // suffix for compilers where version is not clear yet
  // icpc does not specify the version in __VERSION__
#if defined(__INTEL_COMPILER)
  ss << " (icpc version: " << __INTEL_COMPILER << "." << __INTEL_COMPILER_UPDATE << ")";
#elif defined(__FCC_version__)
  ss << " (FCC version: " << __FCC_version__ << ")";
#endif
  // Base case if nothing could be found
  if (ss.str().empty()) {
    ss << "Unknown";
  }

  return ss.str();
}
