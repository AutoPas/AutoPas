/**
 * @file VarVerletTraversalAsBuild.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include "VarVerletTraversalInterface.h"
#include "autopas/containers/verletListsCellBased/verletLists/neighborLists/VerletNeighborListAsBuild.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 *
 * @tparam ParticleCell Needed because every traversal has to be a CellPairTraversal at the moment.
 * @tparam Particle
 * @tparam PairwiseFunctor
 * @tparam useNewton3
 */
template <class ParticleCell, class Particle, class PairwiseFunctor, bool useNewton3>
class VarVerletTraversalAsBuild
    : public VarVerletTraversalInterface<ParticleCell, VerletNeighborListAsBuild<Particle>> {
 public:
  explicit VarVerletTraversalAsBuild(PairwiseFunctor *pairwiseFunctor) : _functor(pairwiseFunctor) {}

  bool usesNewton3() override { return useNewton3; }

  void iterateVerletLists(VerletNeighborListAsBuild<Particle> &neighborList) override {
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

  bool isApplicable() override { return true; }

  TraversalOption getTraversalType() override { return TraversalOption::varVerletTraversalAsBuild; }

 private:
  PairwiseFunctor *_functor;
};

}  // namespace autopas
