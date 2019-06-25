/**
 * @file CudaSoAType.h
 * Provides a struct to easily generate the desired SoAType.
 * @author jspahl
 * @date 10.02.19
 */

#pragma once

#include <tuple>
#include <vector>
#include "autopas/utils/CudaDeviceVector.h"

namespace autopas {
namespace utils {

/**
 * Helper struct to get a the CudaSoAType.
 * The type is defined as CudaSoAType<size_t, double, double, double>::Type;
 * @tparam soatypes the template parameter list of types.
 */
template <typename... soatypes>
struct CudaSoAType {
  /**
   * This is the Type of the SoAType.
   * It is a tuple of aligned vectors.
   */
  typedef std::tuple<autopas::utils::CudaDeviceVector<soatypes>...> Type;
};

}  // namespace utils
}  // namespace autopas
