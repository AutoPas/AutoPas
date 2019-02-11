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

template <class CudaSoAArraysType>
class CudaSoA {
 public:
  /**
   * @brief Default constructor.
   */
  CudaSoA() = default;

  /**
   * @brief delete constructor.
   * @param soa SoA to copy.
   */
  CudaSoA(const CudaSoA& soa) = default;

  template <size_t soaAttribute>
  auto& get() {
    return std::get<soaAttribute>(soaStorageTuple);
  }

  template <size_t soaAttribute>
  const auto& get() const {
    return std::get<soaAttribute>(soaStorageTuple);
  }

 private:
  CudaSoAArraysType soaStorageTuple;
};
}  // namespace autopas
