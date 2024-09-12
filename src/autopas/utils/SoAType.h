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
  using MainType = typename mainPartitionType::Type;

  /**
   * Type of the structure of additional SoA Partitions that the data is stored in.
   */
  using AdditionalType = typename std::tuple<std::vector<typename additionalPartitionTypes::Type>...>;

  /**
   * Types of the additional SoA Partitions. Data is not intended to be stored using this type, this type is to provide
   * a convenient lookup reference for the types of the additional partitions .
   */
  using AdditionalPartitionTypesReference = typename std::tuple<typename additionalPartitionTypes::Type...>;
};

}  // namespace autopas::utils
