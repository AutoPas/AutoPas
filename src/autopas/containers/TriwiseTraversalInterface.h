/**
 * @file TriwiseTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.2024
 */

#pragma once

#include "TraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/TraversalOption.h"

namespace autopas {

/**
 * This interface serves as a common parent class for all traversals specific to triwise interactions.
 */
class TriwiseTraversalInterface : virtual public TraversalInterface {
 public:
  /**
   * Destructor of TraversalInterface
   */
  virtual ~TriwiseTraversalInterface() = default;

  /**
   * Traverses all particle triplets.
   */
  virtual void traverseParticleTriplets() {
    utils::ExceptionHandler::exception(
        "Error: TraversalInterface::traverseParticleTriplets() is "
        "not implemented for this traversal: {}!",
        typeid(*this).name());
  }
};
}  // namespace autopas
