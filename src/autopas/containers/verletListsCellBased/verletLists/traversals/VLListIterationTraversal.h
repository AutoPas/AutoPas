/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @authors jspahl, Alexander-Haberl-TUM
 */

#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam Functor_T The functor that defines the interaction of two or three particles.
 */
template <class ParticleCell_T, class Functor_T>
class VLListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell_T> {
  using ParticleType = typename ParticleCell_T::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param functor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLListIterationTraversal(Functor_T &functor, const DataLayoutOption dataLayout, const bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(functor) {
    if (useNewton3) {
      AutoPasLog(WARN,
                 "VLListIterationTraversal: This traversal is not parallelized for newton3 enabled and will utilize "
                 "only one thread.");
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  /**
   * VL List iteration is always applicable to the domain.
   * @return true
   */
  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  void initTraversal() override {
    if (_dataLayout == DataLayoutOption::soa) {
      auto &cells = *(this->_cells);
      // First, resize the SoA to the required number of elements to store. This avoids resizing successively the SoA in
      // SoALoader.
      std::vector<size_t> offsets(cells.size() + 1);
      std::inclusive_scan(
          cells.begin(), cells.end(), offsets.begin() + 1,
          [](const size_t &partialSum, const auto &cell) { return partialSum + cell.size(); }, 0);

      _soa.resizeArrays(offsets.back());

      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cells.size(); ++i) {
        _functor.SoALoader(cells[i], _soa, offsets[i], /*skipSoAResize*/ true);
      }
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor.SoAExtractor(cell, _soa, offset);
        offset += cell.size();
      }
    }
  }

  void traverseParticles() override {
    if constexpr (utils::isPairwiseFunctor<Functor_T>()) {
      traverseParticlePairs();
    } else if constexpr (utils::isTriwiseFunctor<Functor_T>()) {
      traverseParticleTriplets();
    } else {
      utils::ExceptionHandler::exception(
          "VLListIterationTraversal::traverseParticles(): Functor {} is not of type PairwiseFunctor or TriwiseFunctor.",
          _functor.getName());
    }
  }

  /**
   *  Iterate over all pairs of particles.
   */
  void traverseParticlePairs() {
    auto &neighborList = *(this->_neighborList);
    const size_t numParticles = neighborList.size();
    const auto &indexToParticle = *this->_indexToParticle;

    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          // Each particle i owns its own list slice — no write-conflict between iterations.
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *indexToParticle[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighborsIPtr = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              _functor.AoSFunctor(particleI, *indexToParticle[neighborsIPtr[j]], false);
            }
          }
        } else {
          // Newton3 cannot be parallelized here
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *indexToParticle[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighborsIPtr = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              _functor.AoSFunctor(particleI, *indexToParticle[neighborsIPtr[j]], true);
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
        if (not _useNewton3) {
          AUTOPAS_OPENMP(parallel for schedule(dynamic, std::max(numParticles / (autopas::autopas_get_max_threads() * 10), 1ul)))
          for (size_t i = 0; i < numParticles; ++i) {
            _functor.SoAFunctorVerlet(_soa, i, neighborList.getNeighbors(i), false);
          }
        } else {
          // Newton3 cannot be parallelized here
          for (size_t i = 0; i < numParticles; ++i) {
            _functor.SoAFunctorVerlet(_soa, i, neighborList.getNeighbors(i), true);
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

  /**
   *  Iterate over all triplets of particles.
   */
  void traverseParticleTriplets() {
    auto &neighborList = *(this->_neighborList);
    const size_t numParticles = neighborList.size();
    const auto &indexToParticle = *this->_indexToParticle;
    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          AUTOPAS_OPENMP(parallel for schedule(dynamic))
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particle = *indexToParticle[i];
            if (not particle.isOwned()) {
              // skip Halo particles, as N3 is disabled
              continue;
            }
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighbors = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              ParticleType &neighbor1 = *indexToParticle[neighbors[j]];
              for (size_t k = j + 1; k < numNeighbors; ++k) {
                ParticleType &neighbor2 = *indexToParticle[neighbors[k]];
                _functor.AoSFunctor(particle, neighbor1, neighbor2, false);
              }
            }
          }
        } else {
          // Newton3 cannot be parallelized here
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *indexToParticle[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighbors = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              ParticleType &neighbor1 = *indexToParticle[neighbors[j]];
              for (size_t k = j + 1; k < numNeighbors; ++k) {
                ParticleType &neighbor2 = *indexToParticle[neighbors[k]];
                _functor.AoSFunctor(particleI, neighbor1, neighbor2, true);
              }
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
        utils::ExceptionHandler::exception(
            "VLListIterationTraversal::traverseParticleTriplets(): SoA dataLayout not implemented yet for "
            "triwise interactions.");
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
  Functor_T &_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;
};

}  // namespace autopas