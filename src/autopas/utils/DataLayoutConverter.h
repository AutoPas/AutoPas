/**
 *@file DataLayoutConverter.h
 *@author jspahl
 *@date 2.4.19
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"

namespace autopas::utils {

/**
 * This converts cells to the target data Layout using the given functor
 *
 * @tparam Functor The functor that defines the interaction of two particles.
 */
template <class FunctorSoaWrapper>
class DataLayoutConverter {
 public:
  /**
   * Constructor
   * @tparam Functor Functor Type
   * @param functor responsible for the conversion
   * @param dataLayout The data layout to be used.
   */
  explicit DataLayoutConverter(FunctorSoaWrapper *functor, DataLayoutOption::Value dataLayout)
      : _functor(functor), _dataLayout(dataLayout) {}

  /**
   * loads the target dataLayout in a cell
   * @tparam ParticleCell Cell type
   * @param cell to load the data in
   */
  template <class ParticleCell>
  void loadDataLayout(ParticleCell &cell) {
    switch (_dataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        _functor->SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
        return;
      }
    }
  }

  /**
   * converts the dataLayout to aos
   * @tparam ParticleCell Cell type
   * @param cell to load the data in
   */
  template <class ParticleCell>
  void storeDataLayout(ParticleCell &cell) {
    switch (_dataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        _functor->SoAExtractor(cell, cell._particleSoABuffer, 0);
        return;
      }
    }
  }

 private:
  /**
   *  Functor to convert cells
   */
  FunctorSoaWrapper *_functor;

  const DataLayoutOption::Value _dataLayout;
};

}  // namespace autopas::utils
