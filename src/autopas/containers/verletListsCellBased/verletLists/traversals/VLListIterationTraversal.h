/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @authors jspahl, Alexander-Haberl-TUM
 */

#pragma once

#include "VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
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
   * @param cellsPerDim Dimensions of the cell grid (needed for coloring)
   */
  explicit VLListIterationTraversal(Functor_T &functor, DataLayoutOption dataLayout, bool useNewton3,
                                    const std::array<unsigned long, 3> &cellsPerDim = {0, 0, 0})
      : TraversalInterface(dataLayout, useNewton3), _functor(functor), _cellsPerDim(cellsPerDim) {
    if (useNewton3 and std::ranges::all_of(_cellsPerDim, [](auto x) { return x == 0; })) {
      utils::ExceptionHandler::exception(
          "VLListIterationTraversal: cellsPerDim must be set if newton3 is enabled to avoid race conditions.");
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  /**
   * VL List iteration is always applicable to the domain.
   * @return true
   */
  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  void initTraversal() override {
    auto &cells = *(this->_cells);

    // Pre-compute per-cell offsets so each SoALoader call can be parallelized
    // independently without knowing the preceding cells' sizes. Also needed for coloring.
    std::vector<size_t> offsets(cells.size() + 1);
    offsets[0] = 0;
    for (size_t c = 0; c < cells.size(); ++c) {
      offsets[c + 1] = offsets[c] + cells[c].size();
    }

    if (_useNewton3) {
      // Initialize/clear color cells lists
      for (auto &colorGroup : _colorCells) {
        colorGroup.clear();
      }

      for (size_t c = 0; c < cells.size(); ++c) {
        if (offsets[c] < offsets[c + 1]) {
          auto c3D = utils::ThreeDimensionalMapping::oneToThreeD(c, _cellsPerDim);
          const size_t color = (c3D[0] % 3) + 3 * (c3D[1] % 3) + 9 * (c3D[2] % 3);
          _colorCells[color].push_back({offsets[c], offsets[c + 1]});
        }
      }
    }

    if (_dataLayout == DataLayoutOption::soa) {
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
          AUTOPAS_OPENMP(parallel for schedule(static))
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *indexToParticle[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighborsIPtr = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              _functor.AoSFunctor(particleI, *indexToParticle[neighborsIPtr[j]], false);
            }
          }
        } else {
          // Parallelized AoS with Newton3 using C27 coloring
          for (int color = 0; color < 27; ++color) {
            const auto &cellsOfColor = _colorCells[color];
            AUTOPAS_OPENMP(parallel for schedule(dynamic))
            for (const auto &[pFirst, pEnd] : cellsOfColor) {
              for (size_t i = pFirst; i < pEnd; ++i) {
                ParticleType &particleI = *indexToParticle[i];
                const size_t numNeighbors = neighborList.count(i);
                const size_t *neighborsIPtr = neighborList.begin(i);
                for (size_t j = 0; j < numNeighbors; ++j) {
                  _functor.AoSFunctor(particleI, *indexToParticle[neighborsIPtr[j]], true);
                }
              }
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
          // Parallelized SoA with Newton3 using C27 coloring
          for (int color = 0; color < 27; ++color) {
            const auto &cellsOfColor = _colorCells[color];
            AUTOPAS_OPENMP(parallel for schedule(dynamic))
            for (const auto &[pFirst, pEnd] : cellsOfColor) {
              for (size_t i = pFirst; i < pEnd; ++i) {
                _functor.SoAFunctorVerlet(_soa, i, neighborList.getNeighbors(i), true);
              }
            }
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
          // Parallelized AoS with Newton3 using C27 coloring
          for (int color = 0; color < 27; ++color) {
            const auto &cellsOfColor = _colorCells[color];
            AUTOPAS_OPENMP(parallel for schedule(dynamic))
            for (const auto &[pFirst, pEnd] : cellsOfColor) {
              for (size_t i = pFirst; i < pEnd; ++i) {
                ParticleType &particle = *indexToParticle[i];
                const size_t numNeighbors = neighborList.count(i);
                const size_t *neighbors = neighborList.begin(i);
                for (size_t j = 0; j < numNeighbors; ++j) {
                  ParticleType &neighbor1 = *indexToParticle[neighbors[j]];
                  for (size_t k = j + 1; k < numNeighbors; ++k) {
                    ParticleType &neighbor2 = *indexToParticle[neighbors[k]];
                    _functor.AoSFunctor(particle, neighbor1, neighbor2, true);
                  }
                }
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

  /**
   * Cell block dimensions (needed for C27 coloring).
   */
  std::array<unsigned long, 3> _cellsPerDim;

  /**
   * Particle indices grouped by their cell and then by the cell's C27 color.
   */
  std::array<std::vector<std::pair<size_t, size_t>>, 27> _colorCells;
};

}  // namespace autopas