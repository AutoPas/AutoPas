/**
* @file TriwiseTraversalInterface.h
 * @author S. Newcome
 * @date 17/07/2023
*/

#pragma once

#include "TraversalInterface.h"

namespace autopas {

/**
 * This interface serves as a common parent class for all triwise traversals
 */
class TriwiseTraversalInterface : public TraversalInterface  {
 public:
  /**
   * Traverses all particle triplets.
   */
  virtual void traverseParticleTriplets() = 0;
};
}