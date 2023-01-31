/**
 * @file FlopCounterFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * This class helps in getting the number of performed floating point
 * operations. It is a functor that only calculated the amount of floating point
 * operations.
 * @todo this class currently is limited to the following case:
 *  - constant cutoff radius
 *  - constant amount of floating point operations for one kernel call (distance < cutoff)
 * @todo: we may want the possibility of doing this faster in cases where the number of flops per kernel call is constant
 * @tparam Particle
 * @tparam ForceFunctorType
 */
template <class Particle, class ForceFunctorType>
class FlopCounterFunctor : public Functor<Particle, FlopCounterFunctor<Particle, ForceFunctorType>> {
 public:
  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * constructor of FlopCounterFunctor
   * @param forceFunctor force functor whose flops are counted
   * @param cutoffRadius the cutoff radius
   */
  explicit FlopCounterFunctor<Particle, ForceFunctorType>(ForceFunctorType forceFunctor, double cutoffRadius)
      : autopas::Functor<Particle, FlopCounterFunctor<Particle, ForceFunctorType>>(cutoffRadius),
        _forceFunctor(forceFunctor),
        _cutoffSquare(cutoffRadius * cutoffRadius),
        _distanceCalculations(0ul),
        _kernelCalls(0ul),
        _kernelFlops(0ul) {}

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const auto displacement = utils::ArrayMath::sub(i.getR(), j.getR());
    const auto distanceSquared = utils::ArrayMath::dot(displacement, displacement);
    _distanceCalculations.fetch_add(1, std::memory_order_relaxed);

    if (distanceSquared <= _cutoffSquare) {
      _kernelCalls.fetch_add(1, std::memory_order_relaxed);
      _kernelFlops.fetch_add(_forceFunctor.getNumFlopsPerKernelCall(i, j, newton3), std::memory_order_relaxed);
    }
  }

  /**
   * See Functor::SoAFunctorSingle()
   * @note: SoA variant of FlopCounterFunctor should not be used.
   */
  void SoAFunctorSingle(SoAView<typename Particle::SoAArraysType> soa, bool newton3) final {
    autopas::utils::ExceptionHandler::exception("Use FlopCounterFunctor only with the AoS variant!");
  }

  /**
   * See Functor::SoAFunctorPair()
   * @note: SoA variant of FlopCounterFunctor should not be used.
   */
  void SoAFunctorPair(SoAView<typename Particle::SoAArraysType> soa1, SoAView<typename Particle::SoAArraysType> soa2,
                      bool newton3) final {
    autopas::utils::ExceptionHandler::exception("Use FlopCounterFunctor only with the AoS variant!");
  }

  /**
   * See Functor::SoAFunctorVerlet()
   * @note: SoA variant of FlopCounterFunctor should not be used.
   */
  void SoAFunctorVerlet(SoAView<typename Particle::SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    autopas::utils::ExceptionHandler::exception("Use FlopCounterFunctor only with the AoS variant!");
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
   * get the hit rate of the pair-wise interaction, i.e. the ratio of the number
   * of kernel calls compared to the number of distance calculations
   * @return the hit rate
   */
  double getHitRate() { return static_cast<double>(_kernelCalls) / static_cast<double>(_distanceCalculations); }

  /**
   * get the total number of flops
   * @param numFlopsPerKernelCall
   * @return
   */
  [[nodiscard]] size_t getFlops(size_t numFlopsPerKernelCall) const {
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
