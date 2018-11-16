/**
 * @file DummyTraversal.h
 * @author F. Gratl
 * @date 11/9/18
 */

#pragma once

#include "CellPairTraversal.h"

namespace autopas {

/**
 * Dummy traversal doing nothing.
 *
 * This traversal can be used as a workaround for containers which do not yet use traversals.
 *
 * @tparam ParticleCell
 */
template <class ParticleCell>
class DummyTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor
   * @param dims
   */
  DummyTraversal(const std::array<unsigned long, 3> &dims) : CellPairTraversal<ParticleCell>(dims){};

  ~DummyTraversal() = default;

  TraversalOptions getTraversalType() override;

  bool isApplicable() override;

  void traverseCellPairs(std::vector<ParticleCell> &cells) override;
};

template <class ParticleCell>
TraversalOptions DummyTraversal<ParticleCell>::getTraversalType() {
  return TraversalOptions::dummyTraversal;
}

template <class ParticleCell>
bool DummyTraversal<ParticleCell>::isApplicable() {
  return true;
}
template <class ParticleCell>
void DummyTraversal<ParticleCell>::traverseCellPairs(std::vector<ParticleCell> &cells) {}

}  // namespace autopas
