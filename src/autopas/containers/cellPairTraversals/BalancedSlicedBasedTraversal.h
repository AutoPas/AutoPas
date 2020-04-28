/**
 * @file BalancedSlicedBasedTraversal.h
 *
 * @date 24 Apr 2020
 * @author fischerv
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/containers/cellPairTraversals/SlicedBasedTraversal.h"
#include "autopas/containers/loadEstimators/cellBasedHeuristics.h"

namespace autopas {

/**
 * This class provides a load balanced version of the base sliced traversal
 *
 * The domain is still cut into slices along the longest dimension, but the
 * slices are now chosen, so that the computational load for each slice is
 * roughly equal. Differen heuristics can be chosen to estimate this load.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class BalancedSlicedBasedTraversal
    : public SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> {
 public:
  /**
   * Constructor of the balanced sliced traversal.
   * @copydetails SlicedBasedTraversal::SlicedBasedTraversal()
   * @param heuristic The algorithm used for estimating the load
   */
  explicit BalancedSlicedBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                        const double interactionLength, const std::array<double, 3> &cellLength,
                                        const autopas::loadEstimators::CellBasedHeuristic heuristic)
      : SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                    interactionLength, cellLength),
        _heuristic(heuristic) {}

  /**
   * @copydoc SlicedBasedTraversal::initTraversal()
   * Calculates slice thickness according to estimates loads
   */
  void initTraversal() override {
    this->loadDataLayout();
    // make sure locks and thicknesses are empty
    this->_sliceThickness.clear();
    this->_locks.clear();

    // estimate loads along longest axis
    auto maxDimension = this->_dimsPerLength[0];
    auto maxDimensionLength = this->_cellsPerDimension[this->_dimsPerLength[0]];

    std::vector<unsigned long> loads;

    unsigned long fullLoad = 0;
    std::array<unsigned long, 3> lowerCorner = {0, 0, 0};
    std::array<unsigned long, 3> upperCorner = this->_cellsPerDimension;
    // upper corner is inclusive, so subtract 1 from each coordinate
    upperCorner[0]--;
    upperCorner[1]--;
    upperCorner[2]--;
    for (auto x = 0; x < maxDimensionLength; x++) {
      lowerCorner[maxDimension] = x;
      upperCorner[maxDimension] = x;
      auto load = loadEstimators::estimateCellBasedLoad<ParticleCell>(
          this->_heuristic, *(this->_cells), this->_cellsPerDimension, lowerCorner, upperCorner);
      fullLoad += load;
      loads.push_back(fullLoad);
    }

    auto numSlices = (size_t)autopas_get_max_threads();
    AutoPasLog(debug, "{} threads available.", numSlices);
    auto minSliceThickness = this->_overlapLongestAxis + 1;

    // using greedy algorithm to assign slice thicknesses. May lead to less slices being used.
    unsigned int totalThickness = 0;
    /* minimum load for the next slice. Ideally exactly this load is reached. If this is not possible
     * the remaining load is again divided upon the remaining slices.
     */
    auto min = fullLoad / numSlices;
    for (auto s = 0; s < numSlices; s++) {
      unsigned int thickness;
      if (s == numSlices - 1) {
        thickness = maxDimensionLength - totalThickness;
      } else {
        thickness = minSliceThickness;
        while (totalThickness + thickness + 1 < maxDimensionLength && loads[totalThickness + thickness - 1] < min) {
          thickness++;
        }
      }
      if (totalThickness + thickness > maxDimensionLength || thickness < minSliceThickness) {
        // if minSlicethickness can no longer be satisfied, add remaining space to last slice
        this->_sliceThickness[s - 1] += maxDimensionLength - totalThickness;
        AutoPasLog(debug, "Balanced Sliced traversal only using {} threads because of greedy algorithm.", s);
        numSlices = s;
        break;

      } else {
        totalThickness += thickness;
        this->_sliceThickness.push_back(thickness);
        if (s != numSlices - 1) {
          // add avg of remaining load over remaining threads to min
          min = loads[totalThickness - 1] + (fullLoad - loads[totalThickness - 1]) / (numSlices - s - 1);
        }
      }
    }

    // decreases last _sliceThickness by _overlapLongestAxis to account for the way we handle base cells
    this->_sliceThickness.back() -= this->_overlapLongestAxis;

    this->_locks.resize((numSlices - 1) * this->_overlapLongestAxis);
  }

 protected:
  /**
   * Algorithm to use for estimating load
   */
  loadEstimators::CellBasedHeuristic _heuristic;
};

}  // namespace autopas
