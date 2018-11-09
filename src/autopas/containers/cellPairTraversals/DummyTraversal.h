/**
 * @file DummyTraversal.h
 * @author F. Gratl
 * @date 11/9/18
 */

#pragma once

#include "CellPairTraversal.h"

namespace autopas {
template <class ParticleCell>
class DummyTraversal : public CellPairTraversal<ParticleCell> {

 public:
  DummyTraversal(const std::array<unsigned long, 3> &dims) : CellPairTraversal<ParticleCell>(dims) {};

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
template<class ParticleCell>
void DummyTraversal<ParticleCell>::traverseCellPairs(std::vector<ParticleCell> &cells) {
}

}  // namespace autopas