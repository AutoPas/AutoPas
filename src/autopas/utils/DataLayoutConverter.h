/*
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
 * @tparam DataLayout the Layout to convert to
 */
template <class Functor, DataLayoutOption::Value dataLayout>
class DataLayoutConverter {
 public:
  /**
   * Constructor
   * @tparam Functor Functor Type
   * @param functor responsible for the conversion
   */
  explicit DataLayoutConverter(Functor *functor) : _functor(functor) {}

  /**
   * loads the target dataLayout in a cell
   * @tparam ParticleCell Cell type
   * @param cell to load the data in
   */
  template <class ParticleCell>
  void loadDataLayout(ParticleCell &cell) {
    switch (dataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        _functor->SoALoader(cell, cell._particleSoABuffer, 0);
        return;
      }
      case DataLayoutOption::cuda: {
        _functor->SoALoader(cell, cell._particleSoABuffer, 0);
        _functor->deviceSoALoader(cell._particleSoABuffer, cell._particleSoABufferDevice);
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
    switch (dataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        _functor->SoAExtractor(cell, cell._particleSoABuffer, 0);
        return;
      }
      case DataLayoutOption::cuda: {
        _functor->deviceSoAExtractor(cell._particleSoABuffer, cell._particleSoABufferDevice);
        _functor->SoAExtractor(cell, cell._particleSoABuffer, 0);
        return;
      }
    }
  }

 private:
  /**
   *  Functor to convert cells
   */
  Functor *_functor;
};

}  // namespace autopas::utils
