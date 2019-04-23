/**
 * @file RandomGenerator.cpp
 * @author seckler
 * @date 22.05.18
 */

#include "RandomGenerator.h"

template <typename floatType>
floatType RandomGenerator::fRand(floatType fMin, floatType fMax) {
  floatType f = static_cast<floatType>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

template <typename floatType>
std::array<floatType, 3> RandomGenerator::randomPosition(const std::array<floatType, 3>& boxMin,
                                                         const std::array<floatType, 3>& boxMax) {
  std::array<floatType, 3> r{0, 0, 0};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand<floatType>(boxMin[d], boxMax[d]);
  }
  return r;
}

template float RandomGenerator::fRand<float>(float fMin, float fMax);
template double RandomGenerator::fRand<double>(double fMin, double fMax);

template std::array<float, 3> RandomGenerator::randomPosition<float>(const std::array<float, 3>& boxMin,
                                                                     const std::array<float, 3>& boxMax);
template std::array<double, 3> RandomGenerator::randomPosition<double>(const std::array<double, 3>& boxMin,
                                                                       const std::array<double, 3>& boxMax);
