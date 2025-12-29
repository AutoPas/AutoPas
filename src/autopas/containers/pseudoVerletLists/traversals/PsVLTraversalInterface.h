/**
* @file PsVLTraversalInterface.h
 * @date 08.12.2025
 * @author Lars Doll
 */

#pragma once

#include "autopas/cells/SortedCellView.h"

namespace autopas {

/**
 * This class provides the Traversal Interface for the verlet lists container.
 *
 * The container only accepts traversals in its computeInteractions() method that implements this interface.
 */
template <class ParticleCell>
class PsVLTraversalInterface {
public:
  /**
   * Destructor
   */
  virtual ~PsVLTraversalInterface() = default;


  virtual void setOrientationLists(std::vector<std::vector<SortedCellView<ParticleCell>>> &orientationLists){
    _orientationLists = &orientationLists;
  }

protected:
/**
 * Orientation Lists
 */
std::vector<std::vector<SortedCellView<ParticleCell>>> *_orientationLists = nullptr;
};

}  // namespace autopas