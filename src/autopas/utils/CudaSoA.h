/**
 * @file CudaSoAType.h
 * Provides a struct to easily generate the desired SoAType.
 * @author jspahl
 * @date 10.02.19
 */

#pragma once

#include <tuple>
#include <vector>

#include "autopas/utils/CudaSoAType.h"

namespace autopas {

/**
 * This class stores a soa on the GPU.
 * @tparam CudaSoAArraysType types of the soa vectors
 */
template <class CudaSoAArraysType>
class CudaSoA {
 public:
  /**
   * Default constructor.
   */
  CudaSoA() = default;

  /**
   * delete constructor.
   * @param soa CudaSoA to copy.
   */
  CudaSoA(const CudaSoA &soa) = default;

  /**
   * Get the device vector at the specific entry of the storage
   * @tparam soaAttribute the attribute for which the device vector should be returned
   * @return a reference to the vector for the specific attribute
   */
  template <size_t soaAttribute>
  auto &get() {
    return std::get<soaAttribute>(soaStorageTuple);
  }

  /**
   * @copydoc get()
   * @note const variant
   */
  template <size_t soaAttribute>
  const auto &get() const {
    return std::get<soaAttribute>(soaStorageTuple);
  }

 private:
  CudaSoAArraysType soaStorageTuple;
};
}  // namespace autopas
