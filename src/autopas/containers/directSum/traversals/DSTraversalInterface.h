/**
 * @file DSTraversalInterface.h
 * @author seckler
 * @date 09.01.19
 */

#pragma once

#include <vector>

namespace autopas {

/**
 * Interface for traversals used by the DirectSum container.
 *
 * The container only accepts traversals in its iteratePairwise() method that implement this interface.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class DSTraversalInterface {};

}  // namespace autopas