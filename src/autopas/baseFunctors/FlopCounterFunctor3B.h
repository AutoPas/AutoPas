/**
 * @file FlopCounterFunctor3B.h
 *
 * @author muehlhaeusser
 * @date 23.08.2023
 */

#pragma once

#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * This class helps in getting the number of performed floating point operations. It is a functor that only calculates
 * the amount of floating point operations for a 3-body force functor that requires 3 distance checks between each of
 * the particles, where each distance must be smaller than the cutoff.
 *
 * @todo: we may want the possibility of doing this faster in cases where the number of flops per kernel call is
 * constant
 * @tparam Particle
 * @tparam ForceFunctorType
 */
template <class Particle, class ForceFunctorType>
class FlopCounterFunctor3B : public TriwiseFunctor<Particle, FlopCounterFunctor3B<Particle, ForceFunctorType>> {
 public:
  std::string getName() final { return "FlopCounterFunctor3B"; }

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * constructor of FlopCounterFunctor3B
   * @param forceFunctor 3-body force functor whose flops are counted
   * @param cutoffRadius the cutoff radius
   */
  explicit FlopCounterFunctor3B<Particle, ForceFunctorType>(ForceFunctorType forceFunctor, double cutoffRadius)
      : autopas::TriwiseFunctor<Particle, FlopCounterFunctor3B<Particle, ForceFunctorType>>(cutoffRadius),
        _forceFunctor(forceFunctor),
        _cutoffSquare(cutoffRadius * cutoffRadius),
        _distanceCalculations(0ul),
        _kernelCalls(0ul),
        _kernelFlops(0ul) {}

  void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy() or k.isDummy()) {
      return;
    }
    auto drij = j.getR() - i.getR();
    auto drjk = k.getR() - j.getR();
    auto drki = i.getR() - k.getR();

    double dr2ij = autopas::utils::ArrayMath::dot(drij, drij);
    double dr2jk = autopas::utils::ArrayMath::dot(drjk, drjk);
    double dr2ki = autopas::utils::ArrayMath::dot(drki, drki);

    _distanceCalculations.fetch_add(3, std::memory_order_relaxed);

    if (dr2ij <= _cutoffSquare and dr2jk <= _cutoffSquare and dr2ki <= _cutoffSquare) {
      _kernelCalls.fetch_add(1, std::memory_order_relaxed);
      _kernelFlops.fetch_add(
          _forceFunctor.getNumFlopsPerKernelCall(i.getTypeId(), j.getTypeId(), k.getTypeId(), newton3),
          std::memory_order_relaxed);
    }
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static std::array<typename Particle::AttributeNames, 3> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::posX, Particle::AttributeNames::posY, Particle::AttributeNames::posZ};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static std::array<typename Particle::AttributeNames, 3> getNeededAttr(std::false_type) {
    return getNeededAttr();
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 0>{/*Nothing*/};
  }

  /**
   * get the hit rate of the tri-wise interaction, i.e. the ratio of the number
   * of kernel calls compared to the number of distance checks (Checking all 3 distances)
   * @return the hit rate
   */
  double getHitRate() { return static_cast<double>(_kernelCalls) / static_cast<double>(_distanceCalculations / 3.0); }

  /**
   * get the total number of flops
   * @return
   */
  [[nodiscard]] size_t getFlops() const {
    const auto distFlops = numFlopsPerDistanceCalculation * _distanceCalculations;
    return distFlops + _kernelFlops;
  }

  /**
   * get the number of calculated distance operations
   * @return
   */
  [[nodiscard]] size_t getDistanceCalculations() const { return _distanceCalculations; }

  /**
   * get the number of kernel calls, i.e. the number of pairs of particles with
   * a distance not larger than the cutoff
   * @return
   */
  [[nodiscard]] size_t getKernelCalls() const { return _kernelCalls; }

  /**
   * Get the number of kernel flops, i.e. the flops done after the distance calculation and cutoff comparison.
   * @return
   */
  [[nodiscard]] size_t getKernelFlops() const { return _kernelFlops; }

  /**
   * number of flops for one distance calculation.
   * 3 sub + 3 square + 2 add
   */
  static constexpr size_t numFlopsPerDistanceCalculation = 8ul;

 private:
  double _cutoffSquare;
  std::atomic<size_t> _distanceCalculations, _kernelCalls, _kernelFlops;
  ForceFunctorType _forceFunctor;
};

}  // namespace autopas
