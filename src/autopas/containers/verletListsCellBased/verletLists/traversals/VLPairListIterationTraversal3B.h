/**
 * @file VLPairListIterationTraversal3B.h
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
 * @tparam ParticleCell the type of cells
 * @tparam TriwiseFunctor The functor that defines the interaction of three particles.
 */
template <class ParticleCell, class TriwiseFunctor>
class VLPairListIterationTraversal3B : public TraversalInterface, public VLTraversalInterface<ParticleCell> {
  using ParticleType = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for the Verlet Pair List Iteration Traversal
   * @param triwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLPairListIterationTraversal3B(TriwiseFunctor *triwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(triwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_pair_list_iteration_3b; }

  [[nodiscard]] bool isApplicable() const override { return (not _useNewton3); }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      // First resize the SoA to the required number of elements to store. This avoids resizing successively the SoA in
      // SoALoader.
      std::vector<size_t> offsets(cells.size() + 1);
      std::inclusive_scan(
          cells.begin(), cells.end(), offsets.begin() + 1,
          [](const size_t &partialSum, const auto &cell) { return partialSum + cell.size(); }, 0);

      _soa.resizeArrays(offsets.back());

      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cells.size(); ++i) {
        _functor->SoALoader(cells[i], _soa, offsets[i], /*skipSoAResize*/ true);
      }
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor->SoAExtractor(cell, _soa, offset);
        offset += cell.size();
      }
    }
  }

  void traverseParticles() override {
    auto &aosNeighborPairsLists = *(this->_aosNeighborPairsLists);
    auto &soaNeighborPairsLists = *(this->_soaNeighborPairsLists);
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          size_t buckets = aosNeighborPairsLists.bucket_count();
          /// @todo find a sensible chunk size
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = aosNeighborPairsLists.end(bucketId);
            for (auto bucketIter = aosNeighborPairsLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              ParticleType &particle = *(bucketIter->first);
              if (not particle.isOwned()) {
                // skip Halo particles, as N3 is disabled
                continue;
              }

              for (auto &neighborPairPtr : bucketIter->second) {
                ParticleType &neighbor1 = *(neighborPairPtr.first);
                ParticleType &neighbor2 = *(neighborPairPtr.second);
                _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
              }
            }
          }
        } else {
          utils::ExceptionHandler::exception(
              "VLPairListIterationTraversal3B::traverseParticles(): VLPairListIterationTraversal3B does not support "
              "Newton3.");
        }
        return;
      }
      case DataLayoutOption::soa: {
        if (not _useNewton3) {
          AUTOPAS_OPENMP(parallel for schedule(dynamic, std::max(soaNeighborPairsLists.size() / (autopas::autopas_get_max_threads() * 10), 1ul)))
          for (size_t particleIndex = 0; particleIndex < soaNeighborPairsLists.size(); particleIndex++) {
            _functor->SoAFunctorVerletPair(_soa, particleIndex, soaNeighborPairsLists[particleIndex], _useNewton3);
          }
        } else {
          utils::ExceptionHandler::exception(
              "VLPairListIterationTraversal3B::traverseParticles(): VLPairListIterationTraversal3B does not support "
              "Newton3.");
        }
        return;
      }
      default: {
        utils::ExceptionHandler::exception(
            "VLPairListIterationTraversal3B::traverseParticles(): VerletList dataLayout {} not available", _dataLayout);
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
  SoA<typename ParticleType::SoAArraysType> _soa;
};

}  // namespace autopas