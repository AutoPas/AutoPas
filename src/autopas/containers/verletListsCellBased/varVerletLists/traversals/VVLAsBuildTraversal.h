/**
 * @file VVLAsBuildTraversal.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include "VVLTraversalInterface.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/asBuild/VerletNeighborListAsBuild.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * Traversal for VarVerletLists with VerletNeighborListAsBuild as neighbor list. Every particle pair will be processed
 * by the same thread in the same color as during the build of the neighbor list.
 *
 * @tparam ParticleCell
 * @tparam Particle The particle type used by the neighbor list.
 * @tparam PairwiseFunctor The type of the functor to use for the iteration.
 * @tparam dataLayout The data layout to use.
 * @tparam useNewton3 Whether or not this traversal uses newton 3.
 */
template <class ParticleCell, class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3>
class VVLAsBuildTraversal : public VVLTraversalInterface<VerletNeighborListAsBuild<Particle>>,
                            public TraversalInterface {
 private:
  /**
   * Internal iterate method for AoS.
   * @param neighborList The neighbor list to iterate over.
   */
  void iterateAoS(VerletNeighborListAsBuild<Particle> &neighborList);

  /**
   * Internal iterate method for SoA.
   * @param neighborList The neighbor list to iterate over.
   */
  void iterateSoA(VerletNeighborListAsBuild<Particle> &neighborList);

 public:
  /**
   * The Constructor of VVLAsBuildTraversal.
   * @param pairwiseFunctor The functor to use for the iteration.
   */
  explicit VVLAsBuildTraversal(PairwiseFunctor *pairwiseFunctor) : _functor(pairwiseFunctor), _soa{nullptr} {}

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  void initTraversal() override {
    auto &neighborList = *(this->_neighborList);
    if (dataLayout == DataLayoutOption::soa) {
      _soa = neighborList.loadSoA(_functor);
    }
  }

  void endTraversal() override {
    auto &neighborList = *(this->_neighborList);
    if (dataLayout == DataLayoutOption::soa) {
      neighborList.extractSoA(_functor);
      _soa = nullptr;
    }
  }

  void traverseParticlePairs() override {
    auto &neighborList = *(this->_neighborList);
    switch (dataLayout) {
      case DataLayoutOption::aos:
        iterateAoS(neighborList);
        break;
      case DataLayoutOption::soa:
        iterateSoA(neighborList);
        break;
      default:
        autopas::utils::ExceptionHandler::exception("VVLAsBuildTraversal does not know this data layout!");
    }
  }

  [[nodiscard]] bool isApplicable() const override {
    return dataLayout == DataLayoutOption::soa || dataLayout == DataLayoutOption::aos;
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vvl_as_built; }

 private:
  /**
   * The functor to use for the iteration.
   */
  PairwiseFunctor *_functor;
  /**
   * A pointer to the SoA to iterate over if DataLayout is soa.
   */
  SoA<typename Particle::SoAArraysType> *_soa;
};

template <class ParticleCell, class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3>
void VVLAsBuildTraversal<ParticleCell, Particle, PairwiseFunctor, dataLayout, useNewton3>::iterateAoS(
    VerletNeighborListAsBuild<Particle> &neighborList) {
  const auto &list = neighborList.getAoSNeighborList();

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel num_threads(list[0].size())
#endif
  {
    constexpr int numColors = 8;
    for (int color = 0; color < numColors; color++) {
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(static)
#endif
      for (unsigned int thread = 0; thread < list[color].size(); thread++) {
        const auto &particleToNeighborMap = list[color][thread];
        for (const auto &[particlePtr, neighborPtrList] : particleToNeighborMap) {
          for (auto neighborPtr : neighborPtrList) {
            _functor->AoSFunctor(*particlePtr, *neighborPtr, useNewton3);
          }
        }
      }
    }
  }
}

template <class ParticleCell, class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3>
void VVLAsBuildTraversal<ParticleCell, Particle, PairwiseFunctor, dataLayout, useNewton3>::iterateSoA(
    VerletNeighborListAsBuild<Particle> &neighborList) {
  const auto &soaNeighborList = neighborList.getSoANeighborList();

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel num_threads(soaNeighborList[0].size())
#endif
  {
    constexpr int numColors = 8;
    for (int color = 0; color < numColors; color++) {
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(static)
#endif
      for (unsigned int thread = 0; thread < soaNeighborList[color].size(); thread++) {
        const auto &threadNeighborList = soaNeighborList[color][thread];
        for (const auto &[indexFirst, neighbors] : threadNeighborList) {
          _functor->SoAFunctorVerlet(*_soa, indexFirst, neighbors, useNewton3);
        }
      }
    }
  }
}

}  // namespace autopas
