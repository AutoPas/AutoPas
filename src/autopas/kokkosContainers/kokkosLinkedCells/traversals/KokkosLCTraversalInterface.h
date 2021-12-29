/**
 * @file KokkosLCTraversalInterface.h
 * @author lgaertner
 * @date 29.12.21
 */

#pragma once

namespace autopas {

/**
 * Interface for traversals used by the KokkosLinkedCells container.
 *
 * The container only accepts traversals in its iteratePairwise() method that implement this interface.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class KokkosLCTraversalInterface {};

}  // namespace autopas
