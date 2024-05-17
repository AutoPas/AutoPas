/**
 * @file LCTraversalInterface.h
 * @author seckler
 * @date 09.01.19
 */

#pragma once

#include <vector>

#include "autopas/containers/cellTraversals/CellTraversal.h"

namespace autopas {

/**
 * Interface for traversals used by the LinkedCell class.
 *
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class LCPairTraversalInterface : public PairwiseTraversalInterface {
 public:
//  LCPairTraversalInterface(DataLayoutOption dataLayout, bool useNewton3) : TraversalInterface(dataLayout, useNewton3), PairwiseTraversalInterface(dataLayout, useNewton3) {}

};

template <class ParticleCell>
class LCTriTraversalInterface : public TriwiseTraversalInterface {
 public:
//  LCTriTraversalInterface(DataLayoutOption dataLayout, bool useNewton3) : TriwiseTraversalInterface(dataLayout, useNewton3) {}

};

}  // namespace autopas