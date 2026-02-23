/**
 * @file SlicedC02BasedTraversal.h
 *
 * @date 24 May 2020
 * @author fischerv
 */

#pragma once

#include <algorithm>

#include "SlicedBasedTraversal.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the colored sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into as many slices as possible along this dimension. Unlike the regular
 * sliced traversal, this version uses a 2-coloring to prevent race conditions, instead of
 * locking the starting layers. This could also be describes as a c02-traversal. This class
 * is however not derived from ColorBasedTraversal, as that would not allow varying slice thicknesses,
 * and would prevent us from selecting the dimension in which we cut the slices.
 *
 * @tparam ParticleCell The type of cells.
 * @tparam Functor The functor that defines the interaction between particles.
 */
template <class ParticleCell, class Functor>
class SlicedC02BasedTraversal : public SlicedBasedTraversal<ParticleCell, Functor> {
 public:
  /**
   * Constructor of the colored sliced traversal.
   * @copydetails SlicedBasedTraversal::SlicedBasedTraversal()
   */
  explicit SlicedC02BasedTraversal(const std::array<unsigned long, 3> &dims, Functor *functor,
                                   const double interactionLength, const std::array<double, 3> &cellLength,
                                   DataLayoutOption dataLayout, bool useNewton3, bool spaciallyForward)
      : SlicedBasedTraversal<ParticleCell, Functor>(dims, functor, interactionLength, cellLength, dataLayout,
                                                    useNewton3, spaciallyForward) {}

  /**
   * The main traversal of the colored sliced traversal.
   * This provides the structure of the loops and its parallelization.
   *
   * @copydetails C01BasedTraversal::c01Traversal()
   *
   */
  template <typename LoopBody>
  inline void cSlicedTraversal(LoopBody &&loopBody);

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  [[nodiscard]] bool isApplicable() const override {
    return this->_cellsPerDimension[this->_dimsSortedByLength[0]] >= this->_overlapLongestAxis;
  }

  /**
   * Load Data Layouts and sets up slice thicknesses.
   */
  void initTraversal() override {
    this->loadDataLayout();
    // split domain across its longest dimension
    auto minSliceThickness = this->_overlapLongestAxis;
    this->initSliceThickness(minSliceThickness);
  }
};

template <class ParticleCell, class Functor>
template <typename LoopBody>
void SlicedC02BasedTraversal<ParticleCell, Functor>::cSlicedTraversal(LoopBody &&loopBody) {
  using std::array;

  auto numSlices = this->_sliceThickness.size();
  // check if applicable

  std::array<size_t, 2> overLapps23{this->_overlap[this->_dimsSortedByLength[1]],
                                    this->_overlap[this->_dimsSortedByLength[2]]};

  if (not this->_spaciallyForward) {
    overLapps23 = {0ul, 0ul};
  }

  for (size_t offset = 0; offset < 2; offset++) {
    // although every thread gets exactly one iteration (=slice) this is faster than a normal parallel region
    AUTOPAS_OPENMP(parallel for schedule(dynamic, 1))
    for (size_t slice = offset; slice < numSlices; slice += 2) {
      array<unsigned long, 3> myStartArray{0, 0, 0};
      for (size_t i = 0; i < slice; ++i) {
        myStartArray[this->_dimsSortedByLength[0]] += this->_sliceThickness[i];
      }

      const auto lastLayer = myStartArray[this->_dimsSortedByLength[0]] + this->_sliceThickness[slice];
      for (unsigned long dimSlice = myStartArray[this->_dimsSortedByLength[0]]; dimSlice < lastLayer; ++dimSlice) {
        for (unsigned long dimMedium = 0;
             dimMedium < this->_cellsPerDimension[this->_dimsSortedByLength[1]] - overLapps23[0]; ++dimMedium) {
          for (unsigned long dimShort = 0;
               dimShort < this->_cellsPerDimension[this->_dimsSortedByLength[2]] - overLapps23[1]; ++dimShort) {
            array<unsigned long, 3> idArray = {};
            idArray[this->_dimsSortedByLength[0]] = dimSlice;
            idArray[this->_dimsSortedByLength[1]] = dimMedium;
            idArray[this->_dimsSortedByLength[2]] = dimShort;
            loopBody(idArray[0], idArray[1], idArray[2]);
          }
        }
      }
    }
  }
}

}  // namespace autopas
