/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @authors jspahl, Alexander-Haberl-TUM
 */

#pragma once

#include "VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of two or three particles.
 */
template <class ParticleCell, class Functor>
class VLListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell> {
  using ParticleType = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param functor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLListIterationTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(functor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  [[nodiscard]] bool isApplicable() const override {
    // No parallel version with N3 and no data races is available, hence no N3 is completely disabled.
    return (not _useNewton3) and (_dataLayout == DataLayoutOption::aos or _dataLayout == DataLayoutOption::soa);
  }

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
    if constexpr (utils::isPairwiseFunctor<Functor>()) {
      traverseParticlePairs();
    } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
      traverseParticleTriplets();
    } else {
      utils::ExceptionHandler::exception(
          "VLListIterationTraversal::traverseParticles(): Functor {} is not of type PairwiseFunctor or TriwiseFunctor.",
          _functor->getName());
    }
  }

  void traverseParticlePairs() {
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
        return;
      }
      default: {
        utils::ExceptionHandler::exception(
            "VLListIterationTraversal::traverseParticlePairs(): VerletList dataLayout {} not available", _dataLayout);
      }
    }
  }

  void traverseParticleTriplets() {
    auto &aosNeighborLists = *(this->_aosNeighborLists);
    auto &soaNeighborLists = *(this->_soaNeighborLists);
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          size_t buckets = aosNeighborLists.bucket_count();
          /// @todo find a sensible chunk size
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t bucketId = 0; bucketId < buckets; bucketId++) {
            auto endIter = aosNeighborLists.end(bucketId);
            for (auto bucketIter = aosNeighborLists.begin(bucketId); bucketIter != endIter; ++bucketIter) {
              ParticleType &particle = *(bucketIter->first);
              if (not particle.isOwned()) {
                // skip Halo paraticles, as N3 is disabled
                continue;
              }
              auto neighborPtrIter1 = bucketIter->second.begin();
              for (; neighborPtrIter1 != bucketIter->second.end(); ++neighborPtrIter1) {
                auto neighborPtrIter2 = neighborPtrIter1;
                for (++neighborPtrIter2; neighborPtrIter2 != bucketIter->second.end(); ++neighborPtrIter2) {
                  ParticleType &neighbor1 = *(*neighborPtrIter1);
                  ParticleType &neighbor2 = *(*neighborPtrIter2);
                  _functor->AoSFunctor(particle, neighbor1, neighbor2, false);
                }
              }
            }
          }
        } else {
          for (auto &[particlePtr, neighborPtrList] : aosNeighborLists) {
            ParticleType &particle = *particlePtr;
            if ((not _useNewton3) and (not particle.isOwned())) {
              // skip Halo Particles for N3 disabled
              continue;
            }
            auto neighborPtrIter1 = neighborPtrList.begin();
            for (; neighborPtrIter1 != neighborPtrList.end(); ++neighborPtrIter1) {
              auto neighborPtrIter2 = neighborPtrIter1;
              for (++neighborPtrIter2; neighborPtrIter2 != neighborPtrList.end(); ++neighborPtrIter2) {
                ParticleType &neighbor1 = *(*neighborPtrIter1);
                ParticleType &neighbor2 = *(*neighborPtrIter2);
                _functor->AoSFunctor(particle, neighbor1, neighbor2, _useNewton3);
              }
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
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
        return;
      }
      default: {
        utils::ExceptionHandler::exception(
            "VLListIterationTraversal::traverseParticleTriplets(): VerletList dataLayout {} not available",
            _dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  Functor *_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;
};

}  // namespace autopas