/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @author jspahl
 */
#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class VLListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell> {
  using ParticleType = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param pairwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   * @param cellsPerDim Dimensions of the cell grid (needed for coloring)
   */
  explicit VLListIterationTraversal(PairwiseFunctor &pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3,
                                    const std::array<unsigned long, 3> &cellsPerDim = {0, 0, 0})
      : TraversalInterface(dataLayout, useNewton3), _functor(pairwiseFunctor), _cellsPerDim(cellsPerDim) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  /**
   * VL List iteration is always applicable to the domain.
   * @return true
   */
  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  void initTraversal() override {
    auto &cells = *(this->_cells);

    // Build a flat index→pointer array so the AoS path can resolve CRS indices
    // back to particle references without a hash-map lookup per neighbor.
    // The iteration order must match buildParticleIndex() in VerletLists exactly:
    // iterate over all cells in order, all particles within each cell in order.
    _particlePtrs.clear();
    _particlePtrs.reserve(this->_neighborList->size());
    for (auto &cell : cells) {
      for (auto &p : cell) {
        _particlePtrs.push_back(&p);
      }
    }

    // Pre-compute per-cell offsets so each SoALoader call can be parallelized
    // independently without knowing the preceding cells' sizes. Also needed for coloring.
    std::vector<size_t> offsets(cells.size() + 1);
    offsets[0] = 0;
    for (size_t c = 0; c < cells.size(); ++c) {
      offsets[c + 1] = offsets[c] + cells[c].size();
    }

    // Initialize/clear color indices lists
    for (auto &colorGroup : _colorIndices) {
      colorGroup.clear();
    }

    const bool hasCellsPerDim = (_cellsPerDim[0] > 0 && _cellsPerDim[1] > 0 && _cellsPerDim[2] > 0);
    if (hasCellsPerDim) {
      for (size_t c = 0; c < cells.size(); ++c) {
        auto c3D = utils::ThreeDimensionalMapping::oneToThreeD(c, _cellsPerDim);
        size_t color = (c3D[0] % 2) + 2 * (c3D[1] % 2) + 4 * (c3D[2] % 2);
        for (size_t pIdx = offsets[c]; pIdx < offsets[c + 1]; ++pIdx) {
          _colorIndices[color].push_back(pIdx);
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
    auto &neighborList = *(this->_neighborList);
    const size_t numParticles = neighborList.size();

    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          // Each particle i owns its own list slice — no write conflict between iterations.
          AUTOPAS_OPENMP(parallel for schedule(static))
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *_particlePtrs[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighborsIPtr = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              _functor.AoSFunctor(particleI, *_particlePtrs[neighborsIPtr[j]], false);
            }
          }
        } else {
          // Parallelized AoS with Newton3 using C08 coloring
          const bool hasCellsPerDim = (_cellsPerDim[0] > 0 && _cellsPerDim[1] > 0 && _cellsPerDim[2] > 0);
          if (hasCellsPerDim) {
            for (int color = 0; color < 8; ++color) {
              const auto &indices = _colorIndices[color];
              AUTOPAS_OPENMP(parallel for schedule(static))
              for (size_t k = 0; k < indices.size(); ++k) {
                const size_t i = indices[k];
                ParticleType &particleI = *_particlePtrs[i];
                const size_t numNeighbors = neighborList.count(i);
                const size_t *neighborsIPtr = neighborList.begin(i);
                for (size_t j = 0; j < numNeighbors; ++j) {
                  _functor.AoSFunctor(particleI, *_particlePtrs[neighborsIPtr[j]], true);
                }
              }
            }
          } else {
            // Fallback to serial if cell grid dimensions are not available
            for (size_t i = 0; i < numParticles; ++i) {
              ParticleType &particleI = *_particlePtrs[i];
              const size_t numNeighbors = neighborList.count(i);
              const size_t *neighborsIPtr = neighborList.begin(i);
              for (size_t j = 0; j < numNeighbors; ++j) {
                _functor.AoSFunctor(particleI, *_particlePtrs[neighborsIPtr[j]], true);
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
            _functor.SoAFunctorVerlet(_soa, i, neighborList.begin(i), neighborList.count(i), false);
          }
        } else {
          // Parallelized SoA with Newton3 using C08 coloring
          const bool hasCellsPerDim = (_cellsPerDim[0] > 0 && _cellsPerDim[1] > 0 && _cellsPerDim[2] > 0);
          if (hasCellsPerDim) {
            for (int color = 0; color < 8; ++color) {
              const auto &indices = _colorIndices[color];
              AUTOPAS_OPENMP(parallel for schedule(static))
              for (size_t k = 0; k < indices.size(); ++k) {
                const size_t i = indices[k];
                _functor.SoAFunctorVerlet(_soa, i, neighborList.begin(i), neighborList.count(i), true);
              }
            }
          } else {
            // Fallback to serial if cell grid dimensions are not available
            for (size_t i = 0; i < numParticles; ++i) {
              _functor.SoAFunctorVerlet(_soa, i, neighborList.begin(i), neighborList.count(i), true);
            }
          }
        }
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
  PairwiseFunctor &_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;

  /**
   * Flat array mapping SoA index i -> pointer to particle i.
   * Built in initTraversal() in the same cell/particle iteration order as VerletLists::buildParticleIndex().
   * Used during AoS force computations to resolve CRS neighbor indices to particles.
   */
  std::vector<ParticleType *> _particlePtrs;

  /**
   * Cell block dimensions (needed for C08 coloring).
   */
  std::array<unsigned long, 3> _cellsPerDim;

  /**
   * Particle indices grouped by their cell's C08 color.
   */
  std::array<std::vector<size_t>, 8> _colorIndices;
};

}  // namespace autopas