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
 * Helper struct for getting the types of the main and additional parts of the SoA and passing the type of the SoA.
 *
 * @tparam ArrayTypes the template parameter list of types.
 */
template <typename mainPartitionType, typename... additionalPartitionTypes>
struct SoAType {
  /**
   * Type of the main SoA Partition.
   */
  using MainType = mainPartitionType;

  /**
   * Type of the structure of additional SoA Partitions
   */
  using AdditionalType = typename std::tuple<std::vector<additionalPartitionTypes>...>;
};

}  // namespace autopas::utils
