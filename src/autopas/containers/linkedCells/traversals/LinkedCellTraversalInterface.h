/**
 * @file LinkedCellTraversalInterface.h
 * @author seckler
 * @date 09.01.19
 */

#pragma once

#include <vector>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"

namespace autopas {

/**
 * Interface for traversals used by the LinkedCell class.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class LinkedCellTraversalInterface {};

}  // namespace autopas