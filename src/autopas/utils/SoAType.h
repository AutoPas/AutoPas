/**
 * @file SoAType.h
 * Provides a struct to easily generate the desired SoAPartitionType.
 * @author seckler
 * @date 12.06.18
 */

#pragma once

#include <tuple>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"

namespace autopas::utils {

/**
 * Helper struct to get a the SoAType.
 * The type is defined as SoAType<size_t, double, double, double>::Type;
 * @tparam soatypes the template parameter list of types.
 */
template <typename... soatypes>
struct SoAType {
  /**
   * This is the Type of the SoAType.
   * It is a tuple of aligned vectors.
   */
  using Type = std::tuple<std::vector<soatypes, autopas::AlignedAllocator<soatypes>>...>;
};

}  // namespace autopas::utils
