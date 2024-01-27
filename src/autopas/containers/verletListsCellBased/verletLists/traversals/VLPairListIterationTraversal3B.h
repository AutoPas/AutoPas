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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class TriwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VLPairListIterationTraversal3B : public TraversalInterface<InteractionTypeOption::threeBody>,
                                 public VLTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   */
  explicit VLPairListIterationTraversal3B(TriwiseFunctor *triwiseFunctor) : _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_pair_list_iteration_3b; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override {
    return (not useNewton3) and dataLayout == DataLayoutOption::aos;
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLPairListIterationTraversal3B.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLPairListIterationTraversal3B.");
    }
  }

  void traverseParticleTriplets() override {
    auto &pairwiseAosNeighborLists = *(this->_pairwiseAosNeighborLists);
    switch (dataLayout) {
      case DataLayoutOption::aos: {
        if (not useNewton3) {
          size_t buckets = pairwiseAosNeighborLists.bucket_count();
          /// @todo find a sensible chunk size
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = pairwiseAosNeighborLists.end(bucketId);
            for (auto bucketIter = pairwiseAosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              Particle &particle = *(bucketIter->first);
              if(not particle.isOwned()){
                // skip Halo paraticles, as N3 is disabled
                continue;
              }
              auto neighborPtrPairIter = bucketIter->second.begin();
              for (; neighborPtrPairIter != bucketIter->second.end(); ++neighborPtrPairIter) {
                Particle &neighbor1 = *(neighborPtrPairIter->first);
                Particle &neighbor2 = *(neighborPtrPairIter->second);
                _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
              }
            }
          }
        } else{
          utils::ExceptionHandler::exception("VLPairListIterationTraversal3B does not support Newton3.");
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLPairListIterationTraversal3B.");
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
};

}  // namespace autopas