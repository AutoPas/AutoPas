/**
 * @file PairwiseTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.2024
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/TraversalOption.h"
#include "TraversalInterface.h"

namespace autopas {

/**
 * This interface serves as a common parent class for all traversals specific to pairwise interactions.
 */
class PairwiseTraversalInterface : virtual public TraversalInterface {
 public:
  /**
   * Destructor of TraversalInterface
   */
  virtual ~PairwiseTraversalInterface() = default;

  /**
   * Traverses all particle pairs.
   */
  virtual void traverseParticlePairs() {
    utils::ExceptionHandler::exception(
        "Error: TraversalInterface::traverseParticlePairs() is "
        "not implemented for this traversal: {}!",
        typeid(*this).name());
  };
};
}  // namespace autopas
