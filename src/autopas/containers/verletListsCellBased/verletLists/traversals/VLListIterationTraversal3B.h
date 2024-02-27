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
class VLListIterationTraversal3B : public TraversalInterface<InteractionTypeOption::threeBody>,
                                   public VLTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   */
  explicit VLListIterationTraversal3B(TriwiseFunctor *triwiseFunctor) : _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration_3b; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override {
    return /*(not useNewton3) and*/ dataLayout == DataLayoutOption::aos;
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLListIterationTraversal3B.");
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLListIterationTraversal3B.");
    }
  }

  void traverseParticleTriplets() override {
    auto &aosNeighborLists = *(this->_aosNeighborLists);
    auto &soaNeighborLists = *(this->_soaNeighborLists);
    switch (dataLayout) {
      case DataLayoutOption::aos: {
        if (not useNewton3) {
          size_t buckets = aosNeighborLists.bucket_count();
          /// @todo find a sensible chunk size
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = aosNeighborLists.end(bucketId);
            for (auto bucketIter = aosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              Particle &particle = *(bucketIter->first);
              if (not particle.isOwned()) {
                // skip Halo paraticles, as N3 is disabled
                continue;
              }
              auto neighborPtrIter1 = bucketIter->second.begin();
              for (; neighborPtrIter1 != bucketIter->second.end(); ++neighborPtrIter1) {
                auto neighborPtrIter2 = neighborPtrIter1;
                for (++neighborPtrIter2; neighborPtrIter2 != bucketIter->second.end(); ++neighborPtrIter2) {
                  Particle &neighbor1 = *(*neighborPtrIter1);
                  Particle &neighbor2 = *(*neighborPtrIter2);
                  _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
                }
              }
            }
          }
        } else {
          for (auto &[particlePtr, neighborPtrList] : aosNeighborLists) {
            Particle &particle = *particlePtr;
            if ((not useNewton3) and (not particle.isOwned())) {
              // skip Halo Particles for N3 disabled
              continue;
            }
            auto neighborPtrIter1 = neighborPtrList.begin();
            for (; neighborPtrIter1 != neighborPtrList.end(); ++neighborPtrIter1) {
              auto neighborPtrIter2 = neighborPtrIter1;
              for (++neighborPtrIter2; neighborPtrIter2 != neighborPtrList.end(); ++neighborPtrIter2) {
                Particle &neighbor1 = *(*neighborPtrIter1);
                Particle &neighbor2 = *(*neighborPtrIter2);
                _functor->AoSFunctor(particle, neighbor1, neighbor2, useNewton3);
              }
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception("SoA dataLayout not implemented yet for VLListIterationTraversal3B.");
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