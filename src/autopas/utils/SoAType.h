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
#include "autopas/utils/SoAPartition.h"

namespace autopas::utils {

/**
 * Helper struct for handling the data structure of SoAPartitions for the SoA class, allowing for simpler definitions of
 * SoAs in particle classes.
 *
 * @tparam ArrayTypes the template parameter list of types.
 */
template <typename mainPartitionType, typename... additionalPartitionTypes>
struct SoAType {
  /**
   * Type of the main SoA Partition.
   */
  using MainPartitionType = SoAPartition<mainPartitionType>;

  /**
   * Type of the structure of additional SoA Partitions that the data is stored in.
   */
  using AdditionalPartitionsType = std::tuple<std::vector<SoAPartition<additionalPartitionTypes>>...>;

  using MainPartitionViewType = SoAPartitionView<mainPartitionType>;

  using AdditionalPartitionsViewType = std::tuple<std::vector<SoAPartitionView<additionalPartitionTypes>>...>;

  /**
   * Types of the additional SoA Partitions. Data is not intended to be stored using this type, this type is to provide
   * a convenient lookup reference for the types of the additional partitions .
   */
  using AdditionalPartitionTypesReference = typename std::tuple<additionalPartitionTypes...>;
};

}  // namespace autopas::utils
