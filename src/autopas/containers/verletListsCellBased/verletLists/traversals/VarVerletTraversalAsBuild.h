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
template <class ParticleCell, class Particle, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class VarVerletTraversalAsBuild
    : public VarVerletTraversalInterface<ParticleCell, VerletNeighborListAsBuild<Particle>> {
 private:
  void iterateAoS(VerletNeighborListAsBuild<Particle> &neighborList);

  void iterateSoA(VerletNeighborListAsBuild<Particle> &neighborList);

 public:
  explicit VarVerletTraversalAsBuild(PairwiseFunctor *pairwiseFunctor) : _functor(pairwiseFunctor), _soa{nullptr} {}

  bool usesNewton3() override { return useNewton3; }

  DataLayoutOption getDataLayout() override { return DataLayout; }

  void initVerletTraversal(VerletNeighborListAsBuild<Particle> &neighborList) override {
    if (DataLayout == DataLayoutOption::soa) {
      _soa = neighborList.loadSoA(_functor);
    }
  }

  void endVerletTraversal(VerletNeighborListAsBuild<Particle> &neighborList) override {
    if (DataLayout == DataLayoutOption::soa) {
      neighborList.extractSoA(_functor);
      _soa = nullptr;
    }
  }

  void iterateVerletLists(VerletNeighborListAsBuild<Particle> &neighborList) override {
    switch (DataLayout) {
      case DataLayoutOption::aos:
        iterateAoS(neighborList);
        break;
      case DataLayoutOption::soa:
        iterateSoA(neighborList);
        break;
      default:
        autopas::utils::ExceptionHandler::exception("VarVerletTraversalAsBuild does not know this data layout!");
    }
  }

  bool isApplicable() override { return true; }

  TraversalOption getTraversalType() override { return TraversalOption::varVerletTraversalAsBuild; }

 private:
  PairwiseFunctor *_functor;
  SoA<typename Particle::SoAArraysType> *_soa;
};

template <class ParticleCell, class Particle, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VarVerletTraversalAsBuild<ParticleCell, Particle, PairwiseFunctor, DataLayout, useNewton3>::iterateAoS(
    VerletNeighborListAsBuild<Particle> &neighborList) {
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

template <class ParticleCell, class Particle, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VarVerletTraversalAsBuild<ParticleCell, Particle, PairwiseFunctor, DataLayout, useNewton3>::iterateSoA(
    VerletNeighborListAsBuild<Particle> &neighborList) {
  const auto &soaNeighborList = neighborList.getInternalSoANeighborList();

  {
    constexpr int numColors = 8;
    for (int color = 0; color < numColors; color++) {
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel num_threads(soaNeighborList[color].size())
#endif
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(static)
#endif
      for (unsigned int thread = 0; thread < soaNeighborList[color].size(); thread++) {
        const auto &threadNeighborList = soaNeighborList[color][thread];
        _functor->SoAFunctor(*_soa, threadNeighborList, 0, threadNeighborList.size(), useNewton3);
      }
    }
  }
}

}  // namespace autopas
