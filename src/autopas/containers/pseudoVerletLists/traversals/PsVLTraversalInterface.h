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
 * The container only accepts traversals in its computeInteractions() method that implements this interface.
 */
template <class ParticleCell>
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
  virtual void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell>>> &orientationList) {
    _orientationList = &orientationList;
  }

 protected:
  /**
   * Orientation List: For each cell, 13 sortedCellViews are stored, each of which sorts in the direction of the
   * neighboring cell.
   */
  std::vector<std::vector<SortedCellView<ParticleCell>>> *_orientationList = nullptr;
};

}  // namespace autopas