/**
 * @file MixingTable.h
 * @author muehlhaeusser
 * @date 17.10.2026
 */

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"

namespace mdLib {

  template <size_t numParams, typename FloatType = double, size_t maxStaticTypes = 8>
  class MixingTable {
  public:
    static constexpr size_t staticCapacity = maxStaticTypes * maxStaticTypes;

    struct MixingBuffer {
      alignas(64) std::array<FloatType, staticCapacity> _staticBuffer{};
      std::vector<FloatType, autopas::AlignedAllocator<FloatType>> _dynamicBuffer{};
      const FloatType *__restrict__ _ptr{_staticBuffer.data()};
    };

    MixingTable() noexcept {
      rebindPointers();
    }

    MixingTable(const MixingTable &other) noexcept : _numTypes(other._numTypes), _mixingData(other._mixingData) {
      rebindPointers();
    }

    MixingTable(MixingTable &&other) noexcept : _numTypes(other._numTypes), _mixingData(std::move(other._mixingData)) {
      rebindPointers();
    }

    MixingTable &operator=(const MixingTable &other) noexcept {
      if (this != &other) {
        _numTypes = other._numTypes;
        _mixingData = other._mixingData;
        rebindPointers();
      }
      return *this;
    }

    MixingTable &operator=(MixingTable &&other) noexcept {
      if (this != &other) {
            _numTypes = other._numTypes;
            _mixingData = std::move(other._mixingData);
            rebindPointers();
      }
      return *this;
    }

    ~MixingTable() = default;

    template <typename MixingRulesType>
    void build(const size_t numTypes, MixingRulesType &&mixingRules) {
      resize(numTypes);
        for (size_t i = 0; i < numTypes; ++i) {
          for (size_t j = 0; j < numTypes; ++j) {
            const std::array<FloatType, numParams> mixedParams = mixingRules(i, j);
            const size_t index = i * numTypes + j;
            for (size_t channelIdx = 0; channelIdx < numParams; ++channelIdx) {
              writeEntry(channelIdx, index, mixedParams[channelIdx]);
            }
          }
        }
    }

    template <size_t paramIdx>
    inline FloatType get(const size_t i, const size_t j) const noexcept {
      return _mixingData[paramIdx]._ptr[i * _numTypes + j];
    }

    template <size_t paramIdx>
    inline FloatType get(const size_t i, const size_t j, const size_t k) const noexcept {
      return _mixingData[paramIdx]._ptr[i * _numTypes * _numTypes + j * _numTypes + k];
    }

    template <size_t paramIdx>
    inline const FloatType *getMixingRow(const size_t i) const noexcept {
      return _mixingData[paramIdx]._ptr + i * _numTypes;
    }

    template<size_t paramIdx>
        inline const FloatType *getMixingRow(const size_t i, const size_t j, const size_t k) const noexcept {
      return _mixingData[paramIdx]._ptr + i * _numTypes * _numTypes + j * _numTypes + k;
    }

  private:

    void rebindPointers() noexcept {
      const bool useStatic = (_numTypes <= maxStaticTypes);
      for (size_t i = 0; i < numParams; ++i) {
        _mixingData[i]._ptr = useStatic ? _mixingData[i]._staticBuffer.data() : _mixingData[i]._dynamicBuffer.data();
      }
    }

    void resize(const size_t numTypes) noexcept {
      _numTypes = numTypes;
      if (_numTypes > maxStaticTypes) {
        const size_t requiredCapacity = numTypes * numTypes;
        for (size_t i = 0; i < numParams; ++i) {
          _mixingData[i]._dynamicBuffer.resize(requiredCapacity);
        }
      }
      rebindPointers();
    }

    inline void writeEntry(const size_t channelIndex, const size_t index, const FloatType value) noexcept {
      if (_numTypes <= maxStaticTypes) {
        _mixingData[channelIndex]._staticBuffer[index] = value;
      } else {
        _mixingData[channelIndex]._dynamicBuffer[index] = value;
      }
    }

    size_t _numTypes{0};
    std::array<MixingBuffer, numParams> _mixingData{};
  };

} // namespace mdLib