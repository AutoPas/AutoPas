/**
 * @file LJFunctorKokkos.h
 * @date 11.12.2025
 * @author Luis Gall
 */

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/SoAView.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules)
 * @tparam MemSpace TODO
 * @tparam Particle_T The type of the particle
 * @tparam applyShift Switch for the lj potential to be truncated shifted
 * @tparam useMixing Switch for the functor to be used with multiple particle types
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted
 * @tparam useNewton3 Switch for the functor to support newton3 on, off, or both. See FunctorN3Modes for possible values
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial)
 * @tparam countFLOPs counts FLOPs and hitrate
 * @tparam relevantForTuning Whether of not the auto-tuner should consider this functor
 */
template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorKokkos
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorKokkos<Particle_T, applyShift, useMixing, useNewton3,
                                                                  calculateGlobals, countFLOPs, relevantForTuning>> {
 public:
  /**
   * Structure of the SoAs defined by the particle
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Precision of SoA entries
   */
  using FloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  explicit LJFunctorKokkos(double cutoff, ParticlePropertiesLibrary<FloatPrecision, size_t> &)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorKokkos>(cutoff),
        _cutoffSquared{static_cast<FloatPrecision>(cutoff * cutoff)} {}

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    // No Op, TODO: make sure this is never uses (never!)

    std::cout << "Trying to call non-existing function" << std::endl;
  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    // No Op, TODO: make sure this is never used (also not in remainder traversal)

    if (soa.size() == 0) {
      return;
    }

    std::cout << "Trying to call non-existing function" << std::endl;
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {
    // No Op, TODO: make sure this is never used (also not in remainder traversal)

    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    std::cout << "Trying to call non-existing function" << std::endl;
  }

  KOKKOS_INLINE_FUNCTION
  void ForceKernelKokkos(const FloatPrecision &x1, const FloatPrecision &y1, const FloatPrecision &z1,
                       const autopas::utilsKokkos::KokkosStorage<Particle_T>& storage2, FloatPrecision &fxAcc, FloatPrecision &fyAcc,
                       FloatPrecision &fzAcc, FloatPrecision cutoffSquared, int i, int j) final {
    // const auto owned2 =
    //     soa2.template operator()<
    //         Particle_T::AttributeNames::ownershipState,
    //         true, false>(j);

    // if (owned2 != autopas::OwnershipState::dummy) {

    const auto x2 = storage2.template operator()<Particle_T::AttributeNames::posX, false>(j);

    const auto y2 = storage2.template operator()<Particle_T::AttributeNames::posY, false>(j);

    const auto z2 = storage2.template operator()<Particle_T::AttributeNames::posZ, false>(j);

    ljPair(x1, y1, z1, x2, y2, z2, cutoffSquared, fxAcc, fyAcc, fzAcc);

    // } // owned2 != dummy
  }

  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 9>{
        Particle_T::AttributeNames::id,
        Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,
        Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::forceX,
        Particle_T::AttributeNames::forceY,
        Particle_T::AttributeNames::forceZ,
        Particle_T::AttributeNames::typeId,
        Particle_T::AttributeNames::ownershipState,
    };
  }

  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
  }

  /* Interface required stuff */
  std::string getName() final { return "LJFunctorKokkos"; }

  bool isRelevantForTuning() final { return true; }

  bool allowsNewton3() final {
    // TODO
    return true;
  }

  bool allowsNonNewton3() final {
    // TODO
    return true;
  }

 private:
  KOKKOS_INLINE_FUNCTION
  static void ljPair(FloatPrecision x1, FloatPrecision y1, FloatPrecision z1, FloatPrecision x2, FloatPrecision y2,
                     FloatPrecision z2, FloatPrecision cutoffSquared, FloatPrecision &fx, FloatPrecision &fy,
                     FloatPrecision &fz) {
    const FloatPrecision drX = x1 - x2;
    const FloatPrecision drY = y1 - y2;
    const FloatPrecision drZ = z1 - z2;

    const FloatPrecision dr2 = drX * drX + drY * drY + drZ * drZ;

    if (dr2 <= cutoffSquared && dr2 > 0.) {
      const FloatPrecision sigmaSquared = 1.;
      const FloatPrecision epsilon24 = 24.;

      const FloatPrecision invDr2 = 1. / dr2;

      FloatPrecision lj6 = sigmaSquared * invDr2;
      lj6 = lj6 * lj6 * lj6;

      const FloatPrecision lj12 = lj6 * lj6;
      const FloatPrecision lj12m6 = lj12 - lj6;
      const FloatPrecision fac = epsilon24 * (lj12 + lj12m6) * invDr2;

      fx += fac * drX;
      fy += fac * drY;
      fz += fac * drZ;
    }
  }

 private:
  FloatPrecision _cutoffSquared;
};
}  // namespace mdLib