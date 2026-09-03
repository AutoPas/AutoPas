/**
 * @file VLListIterationC27Traversal.h
 *
 * @date 03.09.2026
 * @author muehlhaeusser
 */
#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides a colored Traversal for the verlet lists container.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor_T>
class VLListIterationC27Traversal : public TraversalInterface, public VLTraversalInterface<ParticleCell_T> {
  using ParticleType = ParticleCell_T::ParticleType;

 public:
  /**
   * Constructor for colored Verlet Traversal
   * @param pairwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   * @param cellsPerDim Dimensions of the cell grid (needed for coloring)
   */
  explicit VLListIterationC27Traversal(PairwiseFunctor_T &pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3,
                                       const std::array<unsigned long, 3> &cellsPerDim)
      : TraversalInterface(dataLayout, useNewton3), _functor(pairwiseFunctor), _cellsPerDim(cellsPerDim) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration_c27; }

  /**
   * VL List iteration C27 is always applicable to the domain.
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
    const auto &indexToParticle = *this->_indexToParticle;

    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        // Parallelized AoS with Newton3 using C27 coloring
        for (int color = 0; color < 27; ++color) {
          const auto &cellsOfColor = _colorCells[color];
            AUTOPAS_OPENMP(parallel for schedule(dynamic))
            for (size_t c = 0; c < cellsOfColor.size(); ++c) {
              const auto &range = cellsOfColor[c];
              for (size_t i = range.first; i < range.second; ++i) {
                ParticleType &particleI = *indexToParticle[i];
                const size_t numNeighbors = neighborList.count(i);
                const size_t *neighborsIPtr = neighborList.begin(i);
                for (size_t j = 0; j < numNeighbors; ++j) {
                  _functor.AoSFunctor(particleI, *indexToParticle[neighborsIPtr[j]], _useNewton3);
                }
              }
            }
        }
        return;
      }

      case DataLayoutOption::soa: {
        // Parallelized SoA with Newton3 using C27 coloring
        for (int color = 0; color < 27; ++color) {
          const auto &cellsOfColor = _colorCells[color];
            AUTOPAS_OPENMP(parallel for schedule(dynamic))
            for (size_t c = 0; c < cellsOfColor.size(); ++c) {
              const auto &range = cellsOfColor[c];
              for (size_t i = range.first; i < range.second; ++i) {
                _functor.SoAFunctorVerlet(_soa, i, neighborList.getNeighbors(i), _useNewton3);
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
  PairwiseFunctor_T &_functor;

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