/**
 * @file DSPairTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.24
 */

#pragma once

#include <vector>

#include "DSTraversalInterface.h"
#include "autopas/containers/PairwiseTraversalInterface.h"

namespace autopas {

/**
 * Interface for pairwise traversals used by the DirectSum container.
 *
 * The container only accepts pairwise traversals in its iteratePairwise() method that implement this interface.
 */
class DSPairTraversalInterface : public DSTraversalInterface, public PairwiseTraversalInterface {};

}  // namespace autopas