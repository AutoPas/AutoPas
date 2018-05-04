/**
 * @file C08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include <utils/WrapOpenMP.h>
#include "CellPairTraversal.h"

namespace autopas {

template <class ParticleCell, class CellFunctor>
class C08Traversal : public CellPairTraversals<ParticleCell, CellFunctor> {
 public:
  explicit C08Traversal(std::vector<ParticleCell> &cells,
  const std::array<unsigned long, 3> &dims,
      CellFunctor *cellfunctor)
  : CellPairTraversals<ParticleCell, CellFunctor>(cells, dims,
      cellfunctor) {
//    rebuild(cells, dims);
//    computeOffsets();
  }
 private:
};
}  // namespace autopas