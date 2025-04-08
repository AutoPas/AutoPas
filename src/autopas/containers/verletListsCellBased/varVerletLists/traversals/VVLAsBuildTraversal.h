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
 * @tparam Particle_T The particle type used by the neighbor list.
 * @tparam PairwiseFunctor The type of the functor to use for the iteration.
 */
template <class ParticleCell, class Particle_T, class PairwiseFunctor>
class VVLAsBuildTraversal : public VVLTraversalInterface<VerletNeighborListAsBuild<Particle_T>>,
                            public TraversalInterface {
 private:
  /**
   * Internal iterate method for AoS.
   * @param neighborList The neighbor list to iterate over.
   */
  void iterateAoS(VerletNeighborListAsBuild<Particle_T> &neighborList);

  /**
   * Internal iterate method for SoA.
   * @param neighborList The neighbor list to iterate over.
   */
  void iterateSoA(VerletNeighborListAsBuild<Particle_T> &neighborList);

 public:
  /**
   * The Constructor of VVLAsBuildTraversal.
   * @param pairwiseFunctor The functor to use for the iteration.
   * @param dataLayout The data layout to use.
   * @param useNewton3 Whether or not this traversal uses newton 3.
   */
  explicit VVLAsBuildTraversal(PairwiseFunctor *pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(pairwiseFunctor), _soa{nullptr} {}

  void initTraversal() override {
    auto &neighborList = *(this->_neighborList);
    if (_dataLayout == DataLayoutOption::soa) {
      _soa = neighborList.loadSoA(_functor);
    }
  }

  void endTraversal() override {
    auto &neighborList = *(this->_neighborList);
    if (_dataLayout == DataLayoutOption::soa) {
      neighborList.extractSoA(_functor);
      _soa = nullptr;
    }
  }

  void traverseParticles() override {
    auto &neighborList = *(this->_neighborList);
    switch (this->_dataLayout) {
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
    return _dataLayout == DataLayoutOption::soa || _dataLayout == DataLayoutOption::aos;
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
  SoA<typename Particle_T::SoAArraysType> *_soa;
};

template <class ParticleCell, class Particle_T, class PairwiseFunctor>
void VVLAsBuildTraversal<ParticleCell, Particle_T, PairwiseFunctor>::iterateAoS(
    VerletNeighborListAsBuild<Particle_T> &neighborList) {
  const auto &list = neighborList.getAoSNeighborList();

  AUTOPAS_OPENMP(parallel num_threads(list[0].size())) {
    constexpr int numColors = 8;
    for (int color = 0; color < numColors; color++) {
      AUTOPAS_OPENMP(for schedule(static))
      for (unsigned int thread = 0; thread < list[color].size(); thread++) {
        const auto &particleToNeighborMap = list[color][thread];
        for (const auto &[particlePtr, neighborPtrList] : particleToNeighborMap) {
          for (auto neighborPtr : neighborPtrList) {
            _functor->AoSFunctor(*particlePtr, *neighborPtr, _useNewton3);
          }
        }
      }
    }
  }
}

template <class ParticleCell, class Particle_T, class PairwiseFunctor>
void VVLAsBuildTraversal<ParticleCell, Particle_T, PairwiseFunctor>::iterateSoA(
    VerletNeighborListAsBuild<Particle_T> &neighborList) {
  const auto &soaNeighborList = neighborList.getSoANeighborList();

  AUTOPAS_OPENMP(parallel num_threads(soaNeighborList[0].size())) {
    constexpr int numColors = 8;
    for (int color = 0; color < numColors; color++) {
      AUTOPAS_OPENMP(for schedule(static))
      for (unsigned int thread = 0; thread < soaNeighborList[color].size(); thread++) {
        const auto &threadNeighborList = soaNeighborList[color][thread];
        for (const auto &[indexFirst, neighbors] : threadNeighborList) {
          _functor->SoAFunctorVerlet(*_soa, indexFirst, neighbors, _useNewton3);
        }
      }
    }
  }
}

}  // namespace autopas
