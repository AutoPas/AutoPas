/**
 * @file LCPairTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.2024
 */

#pragma once

#include <vector>

#include "autopas/containers/PairwiseTraversalInterface.h"
#include "LCTraversalInterface.h"

namespace autopas {

/**
 * Interface for traversals used by the LinkedCell class.
 *
 * The container only accepts traversals in its iteratePairwise() method that implement this interface.
 */
 template <class ParticleCell>
class LCPairTraversalInterface : public PairwiseTraversalInterface,
                                 public LCTraversalInterface {};
}  // namespace autopas