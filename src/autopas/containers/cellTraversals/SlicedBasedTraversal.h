/**
 * @file SlicedBasedTraversal.h
 *
 * @date 09 Jan 2019
 * @author seckler
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides base for locked- and colored sliced traversals.
 *
 * These traversals find the longest dimension of the simulation domain and cut
 * the domain into multiple slices along this dimension. Slices are
 * assigned to the threads in a round robin fashion.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam Functor The functor that defines the interaction between particles.
 */
template <class ParticleCell, class Functor>
class SlicedBasedTraversal : public CellTraversal<ParticleCell>, public TraversalInterface {
 public:
  /**
   * Constructor of the sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param functor The functor that defines the interaction between particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @param spaciallyForward Whether the base step only covers neighboring cells that are spacially forward (for example
   * c08).
   */
  explicit SlicedBasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor,
                                const double interactionLength, const std::array<double, 3> &cellLength,
                                DataLayoutOption dataLayout, bool useNewton3, bool spaciallyForward)
      : CellTraversal<ParticleCell>(dims),
        TraversalInterface(dataLayout, useNewton3),
        _overlap{},
        _dimsSortedByLength{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlapLongestAxis(0),
        _sliceThickness{},
        _spaciallyForward(spaciallyForward),
        _dataLayoutConverter(functor, dataLayout) {
    this->init();
  }

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  [[nodiscard]] bool isApplicable() const override {
    auto minSliceThickness = _overlapLongestAxis + 1;
    auto maxNumSlices = this->_cellsPerDimension[_dimsSortedByLength[0]] / minSliceThickness;
    return maxNumSlices > 0;
  }

  /**
   * Sets up the slice thicknesses to create as many slices as possible while respecting minSliceThickness.
   * @param minSliceThickness
   */
  virtual void initSliceThickness(unsigned long minSliceThickness) {
    auto numSlices = this->_cellsPerDimension[_dimsSortedByLength[0]] / minSliceThickness;
    _sliceThickness.clear();

    // abort if domain is too small -> cleared _sliceThickness array indicates non applicability
    if (numSlices < 1) return;

    _sliceThickness.insert(_sliceThickness.begin(), numSlices, minSliceThickness);
    auto rest = this->_cellsPerDimension[_dimsSortedByLength[0]] - _sliceThickness[0] * numSlices;
    // remaining slices to distribute the remaining layers on
    auto remSlices = std::min(rest, numSlices);
    for (size_t i = 0; i < remSlices; ++i) {
      _sliceThickness[i] += rest / (remSlices - i);
      rest -= rest / (remSlices - i);
    }
    if (_spaciallyForward) {
      // decreases last _sliceThickness by _overlapLongestAxis to account for the way we handle base cells
      _sliceThickness.back() -= _overlapLongestAxis;
    }
  }

  /**
   * Load Data Layouts and sets up slice thicknesses.
   */
  void initTraversal() override {
    loadDataLayout();
    // split domain across its longest dimension
    auto minSliceThickness = _overlapLongestAxis + 1;
    initSliceThickness(minSliceThickness);
  }

  /**
   * Write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
      /// @todo find a condition on when to use omp or when it is just overhead
      AUTOPAS_OPENMP(parallel for num_threads(autopas::autopas_get_preferred_num_threads()))
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

 protected:
  /**
   * Resets the cell structure of the traversal.
   */
  void init();

  /**
   * Load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  virtual void loadDataLayout() {
    if (this->_cells) {
      auto &cells = *(this->_cells);
      /// @todo find a condition on when to use omp or when it is just overhead
      AUTOPAS_OPENMP(parallel for num_threads(autopas::autopas_get_preferred_num_threads()))
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
  }

  /**
   * Update cell information based on VerletClusterLists
   * @param vcl Pointer to the VerletClusterLists object from which cell information is extracted.
   */
  void reinitForVCL(const VerletClusterLists<typename ParticleCell::ParticleType> *vcl) {
    // Reinitialize the sliced traversal with up to date tower information
    const auto towerSideLength = vcl->getTowerSideLength();
    this->_cellLength = {towerSideLength[0], towerSideLength[1], vcl->getBoxMax()[2] - vcl->getBoxMin()[2]};
    const auto towersPerDim = vcl->getTowersPerDimension();
    this->_cellsPerDimension = {towersPerDim[0], towersPerDim[1], 1};

    // reinitialize
    init();
  }

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;

  /**
   * Store ids of dimensions (0, 1, 2 = x, y, z) sorted by in decreasing order of number of cells in that dimension.
   */
  std::array<int, 3> _dimsSortedByLength;

  /**
   * Overlap of interacting cells along the longest axis.
   */
  unsigned long _overlapLongestAxis;

  /**
   * The number of cells per slice in the dimension that was sliced.
   */
  std::vector<unsigned long> _sliceThickness;

  /**
   * Whether the base step only covers neigboring cells tha are spacially forward (for example c08).
   */
  bool _spaciallyForward;

  /**
   * Cell length in CellBlock3D.
   */
  std::array<double, 3> _cellLength;

 private:
  /**
   * Interaction length (cutoff + skin).
   */
  double _interactionLength;

  /**
   * Data Layout Converter to be used with this traversal.
   */
  utils::DataLayoutConverter<Functor> _dataLayoutConverter;
};

template <class ParticleCell, class Functor>
void SlicedBasedTraversal<ParticleCell, Functor>::init() {
  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
    if (not _spaciallyForward) {
      // there is potentially overlap in both directions.
      _overlap[d] *= 2;
    }
  }

  // find longest dimension
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsSortedByLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsSortedByLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsSortedByLength[1] = 3 - (_dimsSortedByLength[0] + _dimsSortedByLength[2]);

  _overlapLongestAxis = _overlap[_dimsSortedByLength[0]];
}

}  // namespace autopas
