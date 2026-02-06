/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @author jspahl
 */
#pragma once

#include <likwid-marker.h>

#include "MortonIndexTraverslInterface.h"
#include "VLTraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/VerletListsLJCompactSoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class VLListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell>, public MortonIndexTraversalInterface {
  using ParticleType = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param pairwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLListIterationTraversal(PairwiseFunctor *pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  [[nodiscard]] bool isApplicable() const override {
    // No parallel version with N3 and no data races is available, hence no N3 is completely disabled.
    return (not _useNewton3) and (_dataLayout == DataLayoutOption::aos or _dataLayout == DataLayoutOption::soa);
  }

  void initTraversal() override {
    auto &cells = *(this->_cells);
      // LIKWID_MARKER_START("AoS to SoA - data copy");
    if (_dataLayout == DataLayoutOption::soa) {
      // First resize the SoA to the required number of elements to store. This avoids resizing successively the SoA in
      // SoALoader.
      _offsets.resize(cells.size() + 1);
      for (size_t i = 0; i < cells.size(); ++i) {
        const size_t cellId = this->_cellsByMortonIndex ? (*this->_cellsByMortonIndex)[i] : i;
        _offsets[i + 1] = _offsets[i] + cells[cellId].size();
      }

      if (this->_useCompactSoA) {
        _compactSoA.resize(_offsets.back());
        AUTOPAS_OPENMP(parallel for)
        for (size_t i = 0; i < cells.size(); ++i) {
          const size_t cellId = this->_cellsByMortonIndex ? (*this->_cellsByMortonIndex)[i] : i;
          _functor->CompactSoALoader(cells[cellId], _compactSoA, _offsets[i], /*skipSoAResize*/ true);
        }

      } else {
        _soa.resizeArrays(_offsets.back());
        AUTOPAS_OPENMP(parallel for)
        for (size_t i = 0; i < cells.size(); ++i) {
          const size_t cellId = this->_cellsByMortonIndex ? (*this->_cellsByMortonIndex)[i] : i;
          _functor->SoALoader(cells[cellId], _soa, _offsets[i], /*skipSoAResize*/ true);
        }
      }
    }
      // LIKWID_MARKER_STOP("AoS to SoA - data copy");
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      if (this->_useCompactSoA) {
        AUTOPAS_OPENMP(parallel for)
        for (size_t i = 0; i < cells.size(); ++i) {
          const size_t cellId = this->_cellsByMortonIndex ? (*this->_cellsByMortonIndex)[i] : i;
          _functor->SoAExtractor(cells[cellId], _compactSoA._soa, _offsets[i]);
        }
      } else {
        AUTOPAS_OPENMP(parallel for)
        for (size_t i = 0; i < cells.size(); ++i) {
          const size_t cellId = this->_cellsByMortonIndex ? (*this->_cellsByMortonIndex)[i] : i;
          _functor->SoAExtractor(cells[cellId], _soa, _offsets[i]);
        }
      }
    }
  }

  void traverseParticles() override {
    auto &aosNeighborLists = *(this->_aosNeighborLists);
    auto &soaNeighborLists = *(this->_soaNeighborLists);
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        // If we use parallelization,
        if (not _useNewton3) {
          size_t buckets = aosNeighborLists.bucket_count();
          /// @todo find a sensible chunk size
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = aosNeighborLists.end(bucketId);
            for (auto bucketIter = aosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              ParticleType &particle = *(bucketIter->first);
              for (auto neighborPtr : bucketIter->second) {
                ParticleType &neighbor = *neighborPtr;
                _functor->AoSFunctor(particle, neighbor, false);
              }
            }
          }
        } else {
          for (auto &[particlePtr, neighborPtrList] : aosNeighborLists) {
            ParticleType &particle = *particlePtr;
            for (auto neighborPtr : neighborPtrList) {
              ParticleType &neighbor = *neighborPtr;
              _functor->AoSFunctor(particle, neighbor, _useNewton3);
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
        // LIKWID_MARKER_START("force calculation SoA");
        if (this->_useOptimizedLJFunctor) {
          if (not _useNewton3) {
            if (this->_useCompactSoA) {
              AUTOPAS_OPENMP(parallel for schedule(dynamic, std::max(soaNeighborLists.size() / (autopas::autopas_get_max_threads() * 10), 1ul)))
              for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
                _functor->SoAFunctorVerletOptimizedCompactSoA(_compactSoA, particleIndex,
                                                              soaNeighborLists[particleIndex], _useNewton3);
              }
            } else {
              /// @todo find a sensible chunk size
              // LIKWID_MARKER_START("force calculation");
              AUTOPAS_OPENMP(parallel for schedule(dynamic, std::max(soaNeighborLists.size() / (autopas::autopas_get_max_threads() * 10), 1ul)))
              for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
                _functor->SoAFunctorVerletOptimized(_soa, particleIndex, soaNeighborLists[particleIndex], _useNewton3);
              }
            }
            // LIKWID_MARKER_STOP("force calculation");
          } else {
            // iterate over SoA
            for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
              _functor->SoAFunctorVerletOptimized(_soa, particleIndex, soaNeighborLists[particleIndex],
                                                        _useNewton3);
            }
          }
        } else {
          if (not _useNewton3) {
            /// @todo find a sensible chunk size
            AUTOPAS_OPENMP(parallel for schedule(dynamic, std::max(soaNeighborLists.size() / (autopas::autopas_get_max_threads() * 10), 1ul)))
            for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
              _functor->SoAFunctorVerlet(_soa, particleIndex, soaNeighborLists[particleIndex], _useNewton3);
            }
          } else {
            // iterate over SoA
            for (size_t particleIndex = 0; particleIndex < soaNeighborLists.size(); particleIndex++) {
              _functor->SoAFunctorVerlet(_soa, particleIndex, soaNeighborLists[particleIndex], _useNewton3);
            }
          }
        }
        // LIKWID_MARKER_STOP("force calculation SoA");
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
  PairwiseFunctor *_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;

  VerletListsLJCompactSoA<ParticleType> _compactSoA;

  std::vector<size_t> _offsets;
};

}  // namespace autopas