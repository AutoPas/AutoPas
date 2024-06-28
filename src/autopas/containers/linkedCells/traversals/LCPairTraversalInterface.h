/**
 * @file LCPairTraversalInterface.h
 * @author muehlhaeusser
 * @date 18.05.2024
 */

#pragma once

#include <vector>

#include "LCTraversalInterface.h"
#include "autopas/containers/PairwiseTraversalInterface.h"

namespace autopas {

/**
 * Interface for pairwise traversals used by the LinkedCell class.
 *
 * The container only accepts pairwise traversals in its iterateInteractions() method that implement this interface.
 */
template <class ParticleCell>
class LCPairTraversalInterface : public PairwiseTraversalInterface, public LCTraversalInterface {};
}  // namespace autopas