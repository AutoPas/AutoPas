/**
 * @file LCTraversalInterface.h
 * @author seckler
 * @date 09.01.19
 */

#pragma once

#include <vector>

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"

namespace autopas {

/**
 * Interface for traversals used by the LinkedCell class.
 *
 * The container only accepts traversals in its iteratePairwise() method that implement this interface.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class LCTraversalInterface {};

}  // namespace autopas