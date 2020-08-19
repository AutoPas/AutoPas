/**
 * @file BalancedSlicedBasedTraversal.h
 *
 * @date 24 Apr 2020
 * @author fischerv
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/containers/cellPairTraversals/BalancedTraversal.h"
#include "autopas/containers/cellPairTraversals/LockedSlicedBasedTraversal.h"
#include "autopas/utils/Timer.h"

namespace autopas {

/**
 * This class provides a load balanced version of the base sliced traversal.
 *
 * The domain is still cut into slices along the longest dimension, but the
 * slices are now chosen, so that the computational load for each slice is
 * roughly equal. Different heuristics can be chosen to estimate this load.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class BalancedSlicedBasedTraversal
    : public LockedSlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
      public BalancedTraversal {
 public:
  /**
   * Constructor of the balanced sliced traversal.
   * @copydetails SlicedBasedTraversal::SlicedBasedTraversal()
   */
  explicit BalancedSlicedBasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                        const double interactionLength, const std::array<double, 3> &cellLength)
      : LockedSlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(
            dims, pairwiseFunctor, interactionLength, cellLength) {
    // As we create exactly one slice per thread, dynamic scheduling makes little sense.
    this->_dynamic = false;
  }

  /**
   * Calculates slice thickness according to estimates loads
   * @param minSliceThickness
   */
  void initSliceThickness(unsigned long minSliceThickness) override {
    // make thicknesses are empty
    this->_sliceThickness.clear();

    // estimate loads along longest axis
    auto maxDimension = this->_dimsPerLength[0];
    auto maxDimensionLength = this->_cellsPerDimension[this->_dimsPerLength[0]];

    std::vector<unsigned long> loads;
    utils::Timer timer;
    timer.start();
    loads.resize(maxDimensionLength);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
    for (auto x = 0; x < maxDimensionLength; x++) {
      std::array<unsigned long, 3> lowerCorner = {0, 0, 0};
      std::array<unsigned long, 3> upperCorner = this->_cellsPerDimension;
      // upper corner is inclusive, so subtract 1 from each coordinate
      upperCorner[0]--;
      upperCorner[1]--;
      upperCorner[2]--;
      lowerCorner[maxDimension] = x;
      upperCorner[maxDimension] = x;
      if (not this->_loadEstimator) {
        autopas::utils::ExceptionHandler::exception(
            "AutoPas internal error: SlicedBalancedBasedTraversal's _loadEstimator is null.");
      }
      auto load = this->_loadEstimator(this->_cellsPerDimension, lowerCorner, upperCorner);
      loads[x] = load;
    }
    for (auto i = 1; i < loads.size(); i++) {
      loads[i] += loads[i - 1];
    }
    auto fullLoad = loads.back();
    auto loadEstimationTime = timer.stop();
    AutoPasLog(debug, "load estimation took {} nanoseconds", loadEstimationTime);

    auto numSlices = (size_t)autopas_get_max_threads();
    AutoPasLog(debug, "{} threads available.", numSlices);
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
        /// @todo reserve to numSlices before push_back
        this->_sliceThickness.push_back(thickness);
        if (s != numSlices - 1) {
          // add avg of remaining load over remaining threads to min
          min = loads[totalThickness - 1] + (fullLoad - loads[totalThickness - 1]) / (numSlices - s - 1);
        }
      }
    }

    // some analysis output that is only relevant when logger is set to debug
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      std::string thicknessStr;
      std::string loadStr;
      auto lastLoad = 0;
      totalThickness = 0;
      for (auto t : this->_sliceThickness) {
        thicknessStr += std::to_string(t) + ", ";
        totalThickness += t;
        loadStr += std::to_string(loads[totalThickness - 1] - lastLoad) + ", ";
        lastLoad = loads[totalThickness - 1];
      }

      /// @todo: use autopas::utils::ArrayUtils::to_string()
      AutoPasLog(debug, "Slice Thicknesses: [{}]", thicknessStr);
      AutoPasLog(debug, "Slice loads: [{}]", loadStr);
    }

    // decreases last _sliceThickness by _overlapLongestAxis to account for the way we handle base cells
    this->_sliceThickness.back() -= this->_overlapLongestAxis;
  }
};

}  // namespace autopas
