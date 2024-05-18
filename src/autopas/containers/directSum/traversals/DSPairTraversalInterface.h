/**
 * @file DSTraversalInterface.h
 * @author seckler
 * @date 09.01.19
 */

#pragma once

#include <vector>
#include "DSTraversalInterface.h"

namespace autopas {

/**
 * Interface for traversals used by the DirectSum container.
 *
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 */
class DSPairTraversalInterface : public DSTraversalInterface,
                                 public PairwiseTraversalInterface {};

}  // namespace autopas