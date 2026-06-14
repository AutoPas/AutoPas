/**
 * @file VLPairListIterationTraversal.h
 *
 * @date 18.12.2023
 * @author Alexander-Haberl-TUM
 */

#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a triwise traversal for the verlet lists container.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam TriwiseFunctor_T The functor that defines the interaction of three particles.
 */
template <class ParticleCell_T, class TriwiseFunctor_T>
class VLPairListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell_T> {
  using ParticleType = ParticleCell_T::ParticleType;

 public:
  /**
   * Constructor for the Verlet Pair List Iteration Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLPairListIterationTraversal(TriwiseFunctor_T &triwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_pair_list_iteration; }

  [[nodiscard]] bool isApplicable() const override {
    return (not _useNewton3) and _dataLayout == DataLayoutOption::aos;
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception(
          "VLPairListIterationTraversal::initTraversal(): SoA dataLayout not implemented yet for "
          "VLPairListIterationTraversal.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception(
          "VLPairListIterationTraversal::endTraversal(): SoA dataLayout not implemented yet for "
          "VLPairListIterationTraversal.");
    }
  }

  void traverseParticles() override {
    if constexpr (utils::isTriwiseFunctor<TriwiseFunctor_T>()) {
      traverseParticleTriplets();
    } else {
      utils::ExceptionHandler::exception(
          "VLPairListIterationTraversal::traverseParticles(): Functor {} is not of type TriwiseFunctor. The "
          "VLPairListIterationTraversal is implemented only for triwise interactions.",
          _functor.getName());
    }
  }

  void traverseParticleTriplets() {
    auto &aosNeighborPairsLists = *(this->_aosNeighborPairsLists);
    auto &particles = *(this->_indexToParticle);
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          const size_t numParticles = aosNeighborPairsLists.size();
          AUTOPAS_OPENMP(parallel for schedule(static))
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particle = *particles[i];
            if (not particle.isOwned()) {
              continue;
            }

            for (const auto &neighborPairPtr : aosNeighborPairsLists.neighborPairsOf(i)) {
              auto &[neighborPtr1, neighborPtr2] = neighborPairPtr;
              _functor.AoSFunctor(particle, *neighborPtr1, *neighborPtr2, false);
            }
          }
        } else {
          utils::ExceptionHandler::exception(
              "VLPairListIterationTraversal::traverseParticleTriplets(): VLPairListIterationTraversal does not support "
              "Newton3.");
        }
        return;
      }
      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception(
            "VLPairListIterationTraversal::traverseParticleTriplets(): SoA dataLayout not implemented yet for "
            "VLPairListIterationTraversal.");
        return;
      }
      default: {
        utils::ExceptionHandler::exception(
            "VLPairListIterationTraversal::traverseParticleTriplets(): VerletList dataLayout {} not available",
            _dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  TriwiseFunctor_T &_functor;
};

}  // namespace autopas