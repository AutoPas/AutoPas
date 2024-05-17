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
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class DSPairTraversalInterface : public PairwiseTraversalInterface {};

template <class ParticleCell>
class DSTriTraversalInterface : public TriwiseTraversalInterface{};

}  // namespace autopas