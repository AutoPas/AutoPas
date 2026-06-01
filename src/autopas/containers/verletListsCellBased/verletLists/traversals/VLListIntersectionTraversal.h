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
  explicit VLListIntersectionTraversal(TriwiseFunctor_T *triwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_intersection; }

  [[nodiscard]] bool isApplicable() const override {
    return (not _useNewton3) and _dataLayout == DataLayoutOption::aos;
  }

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
    auto &crsNeighborList = *(this->_crsNeighborList);
    auto &particles = *(this->_indexToParticle);

    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (this->_useNewton3) {
          // Same limitation as the old implementation.
          utils::ExceptionHandler::exception(
              "VLListIntersectionTraversal::traverseParticles(): "
              "VLListIntersectionTraversal does not support Newton3.");
        }

        const auto &offsets = crsNeighborList.offsets();
        const auto &neighbors = crsNeighborList.neighbors();
        const size_t numParticles = crsNeighborList.size();

        /*
         * IMPORTANT:
         * This comparator matches the old AoS behavior if your AoS lists were sorted
         * as ParticleType* lists. If you instead sort CRS rows numerically by index,
         * remove this comparator and use the default std::set_intersection ordering.
         */
        const auto indexLessByParticlePtr = [&particles](size_t lhs, size_t rhs) {
          return particles[lhs] < particles[rhs];
        };

        AUTOPAS_OPENMP(parallel) {
          // Thread-local buffer for intersection results.
          std::vector<size_t> intersectingNeighbors;

        AUTOPAS_OPENMP(for schedule(dynamic))
        for (size_t particleIndex = 0; particleIndex < numParticles; ++particleIndex) {
          ParticleType &particle = *particles[particleIndex];

          if (not particle.isOwned()) {
            continue;
          }

          const size_t particleNeighborBegin = offsets[particleIndex];
          const size_t particleNeighborEnd = offsets[particleIndex + 1];

          for (size_t p1 = particleNeighborBegin; p1 < particleNeighborEnd; ++p1) {
            const size_t neighbor1Index = neighbors[p1];
            ParticleType &neighbor1 = *particles[neighbor1Index];

            const size_t neighbor1Begin = offsets[neighbor1Index];
            const size_t neighbor1End = offsets[neighbor1Index + 1];

            const size_t suffixBegin = p1 + 1;
            const size_t suffixEnd = particleNeighborEnd;

            const size_t maxIntersectionSize = std::min(suffixEnd - suffixBegin, neighbor1End - neighbor1Begin);

            intersectingNeighbors.clear();
            intersectingNeighbors.reserve(maxIntersectionSize);

            std::set_intersection(neighbors.begin() + suffixBegin, neighbors.begin() + suffixEnd,
                                  neighbors.begin() + neighbor1Begin, neighbors.begin() + neighbor1End,
                                  std::back_inserter(intersectingNeighbors), indexLessByParticlePtr);

            for (const size_t neighbor2Index : intersectingNeighbors) {
              ParticleType &neighbor2 = *particles[neighbor2Index];
              _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
            }
          }
        }
        }

        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception(
            "VLListIntersectionTraversal::traverseParticles(): "
            "SoA dataLayout not implemented yet for VLListIntersectionTraversal.");
        return;
      }

      default: {
        utils::ExceptionHandler::exception(
            "VLListIntersectionTraversal::traverseParticles(): "
            "VerletList dataLayout {} not available",
            this->_dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  TriwiseFunctor_T *_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;
};

}  // namespace autopas