/**
 * @file VLListIterationTraversal3B.h
 *
 * @date 18.12.2023
 * @author Alexander-Haberl-TUM
 */

#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam TriwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class TriwiseFunctor>
class VLPairListIterationTraversal3B : public TraversalInterface,
                                       public VLTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLPairListIterationTraversal3B(TriwiseFunctor *triwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_pair_list_iteration_3b; }

  [[nodiscard]] bool isApplicable() const override {
    return (not _useNewton3) and _dataLayout == DataLayoutOption::aos;
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLPairListIterationTraversal3B.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLPairListIterationTraversal3B.");
    }
  }

  void traverseParticles() override {
    auto &aosNeighborPairsLists = *(this->_aosNeighborPairsLists);
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          size_t buckets = aosNeighborPairsLists.bucket_count();
          /// @todo find a sensible chunk size
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = aosNeighborPairsLists.end(bucketId);
            for (auto bucketIter = aosNeighborPairsLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              Particle &particle = *(bucketIter->first);
              if (not particle.isOwned()) {
                // skip Halo paraticles, as N3 is disabled
                continue;
              }

              for (auto &neighborPtrPair : bucketIter->second) {
                Particle &neighbor1 = *(neighborPtrPair.first);
                Particle &neighbor2 = *(neighborPtrPair.second);
                _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
              }
            }
          }
        } else {
          utils::ExceptionHandler::exception("VLPairListIterationTraversal3B does not support Newton3.");
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLPairListIterationTraversal3B.");
        return;
      }
      default: {
        utils::ExceptionHandler::exception("VerletList dataLayout {} not available", _dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  TriwiseFunctor *_functor;
};

}  // namespace autopas