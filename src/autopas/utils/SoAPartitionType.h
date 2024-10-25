/**
 * @file SoAPartitionType.h
 * Provides a struct to easily generate the desired SoAPartitionType.
 * @author S. J. Newcome (Created from original SoAType.h by seckler)
 * @date 06/09/2024 (Original SoAType created on 12/06/2018)
 */

#pragma once

#include <tuple>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"

namespace autopas::utils {

/**
 * Helper struct to get a the SoAPartitionType.
 * The type is defined as SoAPartitionType<size_t, double, double, double>::Type;
 * @tparam ArrayTypes the template parameter list of types.
 */
template <typename AttributeTypes, size_t... Attributes>
struct SoAPartitionType {
  /**
   * This is the Type of the SoAPartitionType.
   * It is a tuple of aligned vectors.
   */
  using Type = std::tuple<std::vector<std::tuple_element_t<Attributes, AttributeTypes>, autopas::AlignedAllocator<std::tuple_element_t<Attributes, AttributeTypes>>>...>;
};

}  // namespace autopas::utils
