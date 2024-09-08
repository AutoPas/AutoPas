/**
 * @file SoAType.h
 * Provides a struct to easily generate the desired SoAType.
 * @author seckler
 * @date 12.06.18
 */

#pragma once

#include <tuple>
#include <vector>

#include "autopas/utils/SoAPartitionType.h"

namespace autopas::utils {

/**
 * Helper struct to get a the SoAType.
 * The type is defined as SoAType<SoAPartitionType<size_t, double, double, double>, SoAPartitionType<double>>::Type;
 * @tparam ArrayTypes the template parameter list of types.
 */
template <typename mainPartitionType, typename... additionalPartitionTypes>
struct SoAType {
  /**
   * Enum for more verbose accesses to main and additional kinds of partitions.
   */
  enum SoATypePartitionKind : int {main, additional};

  /**
   * This is the Type of the SoAType.
   * It is a pair of the main SoAPartitionType and a tuple of vectors of any additional SoAPartitionTypes
   */
  using Type = std::pair<mainPartitionType, std::tuple<std::vector<additionalPartitionTypes>...>>;
};

}  // namespace autopas::utils
