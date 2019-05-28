/**
 * @file VarVerletTraversalAsBuild.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include "autopas/utils/WrapOpenMP.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/VerletNeighborListAsBuild.h"
#include "VarVerletTraversalInterface.h"

namespace autopas {

template<class Particle, class PairwiseFunctor, bool useNewton3>
class VarVerletTraversalAsBuild : public VarVerletTraversalInterface<VerletNeighborListAsBuild<Particle>> {
 public:
  explicit VarVerletTraversalAsBuild(PairwiseFunctor *pairwiseFunctor)
      : _functor(pairwiseFunctor) {}

  bool usesNewton3() override {
    return useNewton3;
  }

  void iterateVerletLists(VerletNeighborListAsBuild <Particle> &neighborList) override {
    const auto &list = neighborList.getInternalNeighborList();

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel num_threads(list[0].size())
#endif
    {
      constexpr int numColors = 8;
      for (int c = 0; c < numColors; c++) {
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(static)
#endif
        for (unsigned int thread = 0; thread < list[c].size(); thread++) {
          for (auto &pair : list[c][thread]) {
            _functor->AoSFunctor(*(pair.first), *(pair.second), useNewton3);
          }
        }
      }

    }
  }

 private:
  PairwiseFunctor *_functor;
};

} // namespace autopas
