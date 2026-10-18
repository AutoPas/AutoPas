/**
* @file SoAMixingCache.h
 * @author muehlhaeusser
 * @date 18.10.2026
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>
#include <array>
#include "autopas/utils/AlignedAllocator.h"
#include <algorithm>
#include <tuple>


namespace mdLib {

  template <typename FloatType, size_t numParams, size_t maxStaticTypes = 8>
  class SoAMixingCache {
  public:
    SoAMixingCache() noexcept = default;

    SoAMixingCache(const size_t numTypes, const size_t numParticlesJ) noexcept
    {
      for (size_t paramIdx = 0; paramIdx < numParams; ++paramIdx) {
        _dynamicBuffer[paramIdx].resize(numTypes);
        for (size_t t = 0; t < numTypes; ++t) {
          _dynamicBuffer[paramIdx][t].resize(numParticlesJ);
        }
      }
        _heapPopulated = std::make_unique<bool[]>(numTypes);
        std::fill_n(_heapPopulated.get(), numTypes, false);
    }

    template <size_t... ParamIndices, typename TableType>
    auto getPreloadedArrays(
    const TableType &table, size_t typeI, const size_t *const __restrict__ typeptrJ, size_t numParticlesJ) noexcept {
      static_assert(sizeof...(ParamIndices) == numParams);
      bool *populatedArr = _heapPopulated.get();

      if (not populatedArr[typeI]) {
        populatedArr[typeI] = true;
        populateRow<0, ParamIndices...>(table, typeI, typeptrJ, numParticlesJ);
      }

      return [&]<size_t... Is>(std::index_sequence<Is...>) {
        return std::make_tuple(static_cast<const FloatType *__restrict__>
          (_dynamicBuffer[Is][typeI].data())...);
      }(std::make_index_sequence<numParams>{});
    }

  private:
    template <size_t SlotIdx, size_t FirstParam, size_t... RestParams, typename TableType>
  void populateRow(const TableType &table, size_t typeI,
                                          const size_t *const __restrict__ typeptrJ, size_t numParticlesJ) noexcept {
      FloatType *__restrict__ dst = _dynamicBuffer[SlotIdx][typeI].data();
      const FloatType *const __restrict__ srcRow = table.template getMixingRow<FirstParam>(typeI);

#pragma omp simd
      for (size_t j = 0; j < numParticlesJ; ++j) {
        dst[j] = srcRow[typeptrJ[j]];
      }

      if constexpr (sizeof...(RestParams) > 0) {
        populateRow<SlotIdx + 1, RestParams...>(table, typeI, typeptrJ, numParticlesJ);
      }
    }


    std::array<std::vector<std::vector<FloatType, autopas::AlignedAllocator<FloatType, 64>>>, numParams> _dynamicBuffer{};
    std::unique_ptr<bool[]> _heapPopulated{nullptr};
  };

} // namespace mdLib