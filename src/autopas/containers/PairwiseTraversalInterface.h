/**
* @file PairwiseTraversalInterface.h
 * @author S. Newcome
 * @date 17/07/2023
*/

#pragma once

#include "TraversalInterface.h"

namespace autopas {

/**
 * This interface serves as a common parent class for all pairwise traversals
 */
class PairwiseTraversalInterface : public TraversalInterface {
 public:
  /**
   * Traverses all particle pairs.
   */
  virtual void traverseParticlePairs() = 0;
};
}