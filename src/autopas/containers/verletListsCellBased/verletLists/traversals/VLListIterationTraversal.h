/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @author jspahl
 */

#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VLListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param pairwiseFunctor Functor to be used with this Traversal
   */
  explicit VLListIterationTraversal(PairwiseFunctor *pairwiseFunctor) : _functor(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override {
    return (not useNewton3) and (dataLayout == DataLayoutOption::aos or dataLayout == DataLayoutOption::soa);
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      // First resize the SoA to the required number of elements to store. This avoids resizing successively the SoA in
      // SoALoader.
      size_t numParticles = 0;
#ifdef AUTOPAS_OPENMP
      // this heuristik was taken from CellBasedParticleContainer::size()
      const int numThreads = std::clamp(static_cast<int>(cells.size() / 100000), 1, omp_get_max_threads());
#pragma omp parallel for num_threads(numThreads) reduction(+ : numParticles)
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        numParticles += cells[i].size();
      }
      _soa.resizeArrays(numParticles);

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for num_threads(numThreads)
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        const auto offset = std::accumulate(cells.begin(), cells.begin() + i, 0,
                                            [](auto &acc, auto &cell) { return acc + cell.size(); });
        _functor->SoALoader(cells[i], _soa, offset, /*skipSoAResize*/ true);
      }
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (dataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor->SoAExtractor(cell, _soa, offset);
        offset += cell.size();
      }
    }
  }

  void traverseParticlePairs() override {
    auto &aosNeighborLists = *(this->_aosNeighborLists);
    auto &soaNeighborLists = *(this->_soaNeighborLists);
    switch (dataLayout) {
      case DataLayoutOption::aos: {
#if defined(AUTOPAS_OPENMP)
        if (not useNewton3) {
          size_t buckets = aosNeighborLists.bucket_count();
          /// @todo find a sensible chunk size
#pragma omp parallel for schedule(dynamic)
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = aosNeighborLists.end(bucketId);
            for (auto bucketIter = aosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              Particle &particle = *(bucketIter->first);
              for (auto neighborPtr : bucketIter->second) {
                Particle &neighbor = *neighborPtr;
                _functor->AoSFunctor(particle, neighbor, false);
              }
            }
          }
        } else
#endif
        {
          for (auto &[particlePtr, neighborPtrList] : aosNeighborLists) {
            Particle &particle = *particlePtr;
            for (auto neighborPtr : neighborPtrList) {
              Particle &neighbor = *neighborPtr;
              _functor->AoSFunctor(particle, neighbor, useNewton3);
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
#if defined(AUTOPAS_OPENMP)
        if (not useNewton3) {
          /// @todo find a sensible chunk size
          const size_t chunkSize = std::max(soaNeighborLists.size() / (omp_get_max_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
          for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
            _functor->SoAFunctorVerlet(_soa, particleIndex, soaNeighborLists[particleIndex], useNewton3);
          }
        } else
#endif
        {
          // iterate over SoA
          for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
            _functor->SoAFunctorVerlet(_soa, particleIndex, soaNeighborLists[particleIndex], useNewton3);
          }
        }
        return;
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
  PairwiseFunctor *_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename Particle::SoAArraysType> _soa;
};

}  // namespace autopas