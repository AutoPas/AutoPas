/**
 * @file CSlicedBasedTraversal.h
 *
 * @date 24 May 2020
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the colored sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain in one slice (block) per thread along this dimension. Unlike the regular
 * sliced traversal, this version uses a 2-coloring to prevent race conditions, instead of
 * locking the starting layers. This could also be describes as a c02-traversal. This class
 * is however not derived from CBasedTraversal, as that would not allow varying slice thicknesses,
 * and would prevent us from selecting the dimension in which we cut the slices.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class CSlicedBasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the colored sliced traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit CSlicedBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                 const double interactionLength, const std::array<double, 3> &cellLength)
      : CellPairTraversal<ParticleCell>(dims),
        _overlap{},
        _dimsPerLength{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlapLongestAxis(0),
        _sliceThickness{},
        _dataLayoutConverter(pairwiseFunctor) {
    init(dims);
  }

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  bool isApplicable() const override {
    return not(dataLayout == DataLayoutOption::cuda) and
           this->_cellsPerDimension[_dimsPerLength[0]] >= _overlapLongestAxis;
  }

  /**
   * Load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      /// @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
    // split domain across its longest dimension

    auto minSliceThickness = _overlapLongestAxis;
    auto numSlices = static_cast<size_t>(this->_cellsPerDimension[_dimsPerLength[0]] / minSliceThickness);

    _sliceThickness.clear();

    // abort if domain is too small -> cleared _sliceThickness array indicates non applicability
    if (numSlices < 1) return;

    _sliceThickness.insert(_sliceThickness.begin(), numSlices, minSliceThickness);
    auto rest = this->_cellsPerDimension[_dimsPerLength[0]] - _sliceThickness[0] * numSlices;
    for (size_t i = 0; i < rest; ++i) ++_sliceThickness[i];
    // decreases last _sliceThickness by _overlapLongestAxis to account for the way we handle base cells
    _sliceThickness.back() -= _overlapLongestAxis;
  }

  /**
   * Write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      /// @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

 protected:
  /**
   * Resets the cell structure of the traversal.
   * @param dims
   */
  void init(const std::array<unsigned long, 3> &dims);

  /**
   * The main traversal of the colored sliced traversal.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   *
   * @tparam allCells Defines whether or not to iterate over all cells with the loop body given as argument. By default
   * (allCells=false) it will not iterate over all cells and instead skip the last few cells, because they will be
   * covered by the base step. If you plan to use the default base step of the traversal on this function, use
   * allCells=false, if you plan to just iterate over all cells, e.g., to iterate over verlet lists saved within the
   * cells, use allCells=true. For the sliced step if allCells is false, iteration will not occur over the last layer of
   * cells (for _overlap=1) (in x, y and z direction).
   */
  template <bool allCells = false, typename LoopBody>
  inline void cSlicedTraversal(LoopBody &&loopBody);

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  std::array<unsigned long, 3> _overlap;

 private:
  /**
   * Store ids of dimensions ordered by number of cells per dimensions.
   */
  std::array<int, 3> _dimsPerLength;

  /**
   * Interaction length (cutoff + skin).
   */
  double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  std::array<double, 3> _cellLength;

  /**
   * Overlap of interacting cells along the longest axis.
   */
  unsigned long _overlapLongestAxis;

  /**
   * The number of cells per slice in the dimension that was sliced.
   */
  std::vector<unsigned long> _sliceThickness;

  /**
   * Data Layout Converter to be used with this traversal.
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void CSlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::init(
    const std::array<unsigned long, 3> &dims) {
  for (unsigned int d = 0; d < 3; d++) {
    _overlap[d] = std::ceil(_interactionLength / _cellLength[d]);
  }

  // find longest dimension
  auto minMaxElem = std::minmax_element(this->_cellsPerDimension.begin(), this->_cellsPerDimension.end());
  _dimsPerLength[0] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.second);
  _dimsPerLength[2] = (int)std::distance(this->_cellsPerDimension.begin(), minMaxElem.first);
  _dimsPerLength[1] = 3 - (_dimsPerLength[0] + _dimsPerLength[2]);

  _overlapLongestAxis = _overlap[_dimsPerLength[0]];
}
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
template <bool allCells, typename LoopBody>
void CSlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::cSlicedTraversal(
    LoopBody &&loopBody) {
  using std::array;

  auto numSlices = _sliceThickness.size();
  // 0) check if applicable

  std::array<size_t, 2> overLapps23{_overlap[_dimsPerLength[1]], _overlap[_dimsPerLength[2]]};
  if (allCells) {
    overLapps23 = {0ul, 0ul};
    _sliceThickness.back() += _overlapLongestAxis;
  }

  for (size_t offset = 0; offset < 2; offset++) {
#ifdef AUTOPAS_OPENMP
// although every thread gets exactly one iteration (=slice) this is faster than a normal parallel region
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t slice = offset; slice < numSlices; slice += 2) {
      array<unsigned long, 3> myStartArray{0, 0, 0};
      for (size_t i = 0; i < slice; ++i) {
        myStartArray[_dimsPerLength[0]] += _sliceThickness[i];
      }

      const auto lastLayer = myStartArray[_dimsPerLength[0]] + _sliceThickness[slice];
      for (unsigned long dimSlice = myStartArray[_dimsPerLength[0]]; dimSlice < lastLayer; ++dimSlice) {
        for (unsigned long dimMedium = 0; dimMedium < this->_cellsPerDimension[_dimsPerLength[1]] - overLapps23[0];
             ++dimMedium) {
          for (unsigned long dimShort = 0; dimShort < this->_cellsPerDimension[_dimsPerLength[2]] - overLapps23[1];
               ++dimShort) {
            array<unsigned long, 3> idArray = {};
            idArray[_dimsPerLength[0]] = dimSlice;
            idArray[_dimsPerLength[1]] = dimMedium;
            idArray[_dimsPerLength[2]] = dimShort;
            loopBody(idArray[0], idArray[1], idArray[2]);
          }
        }
      }
    }
  }
  if (allCells) {
    _sliceThickness.back() -= _overlapLongestAxis;
  }
}

}  // namespace autopas
