/**
 * @file PsVLTraversalInterface.h
 * @date 08.12.2025
 * @author Lars Doll
 */

#pragma once

#include "autopas/cells/SortedCellView.h"

namespace autopas {

/**
 * This class provides the Traversal Interface for the pseudoVerletLists container.
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 */
template <class ParticleCell_T>
class PsVLTraversalInterface {
 public:
  /**
   * Destructor.
   */
  virtual ~PsVLTraversalInterface() = default;

  /**
   * Sets the orientationList.
   * @param orientationList
   */
  virtual void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell_T>>> &orientationList) = 0;
};

}  // namespace autopas