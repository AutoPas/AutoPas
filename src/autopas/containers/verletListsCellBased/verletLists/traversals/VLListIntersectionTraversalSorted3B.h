/**
 * @file VLListIntersectionTraversalSorted3B.h
 *
 * @date 18.12.2023
 * @author Alexander-Haberl-TUM
 */

#pragma once

#include <algorithm>
#include <memory>

#include "VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 * It requires the NeighborLists to be sorted.
 *
 * @tparam ParticleCell the type of cells
 * @tparam TriwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class TriwiseFunctor>
class VLListIntersectionTraversalSorted3B : public TraversalInterface, public VLTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLListIntersectionTraversalSorted3B(TriwiseFunctor *triwiseFunctor, DataLayoutOption dataLayout,
                                               bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override {
    return TraversalOption::vl_list_intersection_sorted_3b;
  }

  [[nodiscard]] bool isApplicable() const override {
    return (not _useNewton3) and _dataLayout == DataLayoutOption::aos;
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception(
          "VLListIntersectionTraversalSorted3B::initTraversal(): SoA dataLayout not implemented yet for "
          "VLListIntersectionTraversalSorted3B.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception(
          "VLListIntersectionTraversalSorted3B::endTraversal(): SoA dataLayout not implemented yet for "
          "VLListIntersectionTraversalSorted3B.");
    }
  }

  void traverseParticles() override {
    auto &aosNeighborLists = *(this->_aosNeighborLists);
    auto &soaNeighborLists = *(this->_soaNeighborLists);
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        /// @todo add parelelization
        if (not _useNewton3) {
          size_t buckets = aosNeighborLists.bucket_count();

          AUTOPAS_OPENMP(parallel) {
            // create a buffer per Thread for all intersections
            auto intersectingNeighbors = std::vector<Particle *>();

            AUTOPAS_OPENMP(for schedule(dynamic))
            for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
              auto endIter = aosNeighborLists.end(bucketId);
              for (auto bucketIter = aosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
                Particle &particle = *(bucketIter->first);
                if (not particle.isOwned()) {
                  // skip Halo particles as N3 is disabled
                  continue;
                }
                auto &neighborList = bucketIter->second;
                auto neighborPtrIter1 = neighborList.begin();
                for (; neighborPtrIter1 != neighborList.end(); ++neighborPtrIter1) {
                  Particle &neighbor1 = *(*neighborPtrIter1);
                  auto &neighborList1 = (aosNeighborLists.find(&neighbor1))->second;

                  std::size_t maxIntersectionSize = std::min(neighborList1.size(), neighborList.size());
                  //  make sure the buffer has enough space
                  intersectingNeighbors.reserve(maxIntersectionSize);

                  auto intersectionIter = neighborPtrIter1;
                  auto intersectionEndIter =
                      std::set_intersection(++intersectionIter, neighborList.end(), neighborList1.begin(),
                                            neighborList1.end(), std::back_inserter(intersectingNeighbors));

                  for (auto neighborPtrIter2 = intersectingNeighbors.begin();
                       neighborPtrIter2 != intersectingNeighbors.end(); ++neighborPtrIter2) {
                    Particle &neighbor2 = *(*neighborPtrIter2);
                    _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
                  }

                  // clear buffer for next loop-iteration
                  intersectingNeighbors.clear();
                }
              }
            }
          }
        } else {
          // list intersection does not work with the current way neighborlists are built for N3 case
          utils::ExceptionHandler::exception(
              "VLListIntersectionTraversalSorted3B::traverseParticles(): VLListIntersectionTraversalSorted3B does not "
              "support Newton3.");
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception(
            "VLListIntersectionTraversalSorted3B::traverseParticles(): SoA dataLayout not implemented yet for "
            "VLListIntersectionTraversalSorted3B.");
        return;
      }
      default: {
        utils::ExceptionHandler::exception(
            "VLListIntersectionTraversalSorted3B::traverseParticles(): VerletList dataLayout {} not available",
            _dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  TriwiseFunctor *_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename Particle::SoAArraysType> _soa;
};

}  // namespace autopas