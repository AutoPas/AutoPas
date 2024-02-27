/**
 * @file VLListIterationTraversal3B.h
 *
 * @date 18.12.2023
 * @author Alexander-Haberl-TUM
 */

#pragma once

#include <algorithm>
#include <memory>

#include "VLTraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 * It requires the NeighborLists to be sorted.
 *
 * @tparam ParticleCell the type of cells
 * @tparam TriwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class TriwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VLListIntersectionTraversalSorted3B : public TraversalInterface<InteractionTypeOption::threeBody>,
                                            public VLTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   */
  explicit VLListIntersectionTraversalSorted3B(TriwiseFunctor *triwiseFunctor) : _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override {
    return TraversalOption::vl_list_intersection_sorted_3b;
  }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override { return (not useNewton3) and dataLayout == DataLayoutOption::aos; }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLListIntersectionTraversalSorted3B.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLListIntersectionTraversalSorted3B.");
    }
  }

  void traverseParticleTriplets() override {
    auto &aosNeighborLists = *(this->_aosNeighborLists);
    auto &soaNeighborLists = *(this->_soaNeighborLists);
    switch (dataLayout) {
      case DataLayoutOption::aos: {
        /// @todo add parelelization
        if (not useNewton3) {
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
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception(
            "SoA dataLayout not implemented yet for VLListIntersectionTraversalSorted3B.");
      }
      default: {
        utils::ExceptionHandler::exception("VerletList dataLayout {} not available", dataLayout);
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