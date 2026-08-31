/**
 * @file VLListIntersectionTraversal.h
 *
 * @date 18.12.2023
 * @author Alexander-Haberl-TUM
 */

#pragma once

#include <algorithm>
#include <memory>

#include "VLTraversalInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 * It requires the NeighborLists to be sorted.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam TriwiseFunctor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class TriwiseFunctor_T>
class VLListIntersectionTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell_T> {
  using ParticleType = typename ParticleCell_T::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLListIntersectionTraversal(TriwiseFunctor_T &triwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_intersection; }

  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception(
          "VLListIntersectionTraversal::initTraversal(): SoA dataLayout not implemented yet for "
          "VLListIntersectionTraversal.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception(
          "VLListIntersectionTraversal::endTraversal(): SoA dataLayout not implemented yet for "
          "VLListIntersectionTraversal.");
    }
  }

  void traverseParticles() override {
    auto &neighborList = *(this->_neighborList);
    const size_t numParticles = neighborList.size();
    const auto &indexToParticle = *this->_indexToParticle;
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          AUTOPAS_OPENMP(parallel) {
            // create a buffer per Thread for all intersections
            std::vector<size_t> intersectingNeighbors;

            AUTOPAS_OPENMP(for schedule(dynamic))
            for (size_t i = 0; i < numParticles; ++i) {
              ParticleType &particle = *indexToParticle[i];
              if (not particle.isOwned()) {
                // skip Halo particles as N3 is disabled
                continue;
              }
              const size_t numNeighborsI = neighborList.count(i);
              const size_t *neighborsI = neighborList.begin(i);
              for (size_t j = 0; j < numNeighborsI; ++j) {
                const size_t neighbor1Idx = neighborsI[j];
                ParticleType &neighbor1 = *indexToParticle[neighbor1Idx];
                const size_t numNeighborsJ = neighborList.count(neighbor1Idx);
                const size_t *neighborsJ = neighborList.begin(neighbor1Idx);

                size_t maxIntersectionSize = std::min(numNeighborsI - j - 1, numNeighborsJ);
                intersectingNeighbors.clear();
                intersectingNeighbors.reserve(maxIntersectionSize);

                std::set_intersection(neighborsI + j + 1, neighborsI + numNeighborsI, neighborsJ,
                                      neighborsJ + numNeighborsJ, std::back_inserter(intersectingNeighbors));

                for (size_t kIdx : intersectingNeighbors) {
                  ParticleType &neighbor2 = *indexToParticle[kIdx];
                  _functor.AoSFunctor(particle, neighbor1, neighbor2, false);
                }
              }
            }
          }
        } else {
          // list intersection does not work with the current way neighborlists are built for N3 case
          utils::ExceptionHandler::exception(
              "VLListIntersectionTraversal::traverseParticles(): VLListIntersectionTraversal does not "
              "support Newton3.");
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception(
            "VLListIntersectionTraversal::traverseParticles(): SoA dataLayout not implemented yet for "
            "VLListIntersectionTraversal.");
        return;
      }
      default: {
        utils::ExceptionHandler::exception(
            "VLListIntersectionTraversal::traverseParticles(): VerletList dataLayout {} not available", _dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  TriwiseFunctor_T &_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;
};

}  // namespace autopas