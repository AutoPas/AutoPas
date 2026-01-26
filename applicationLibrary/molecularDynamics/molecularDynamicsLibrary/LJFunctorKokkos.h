/**
 * @file FunctorKokkos.h
 * @date 11.12.2025
 * @author Luis Gall
 */

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/SoAView.h"

namespace mdLib {
template <class MemSpace, class Particle_T, bool applyShift = false, bool useMixing = false,
    autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
    bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorKokkos : public autopas::PairwiseFunctor<Particle_T, LJFunctorKokkos<MemSpace, Particle_T, applyShift, useMixing, useNewton3,
        calculateGlobals, countFLOPs, relevantForTuning>, MemSpace> {

public:
  using SoAArraysType = typename Particle_T::SoAArraysType;
  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  explicit LJFunctorKokkos(double cutoff, ParticlePropertiesLibrary<FloatPrecision, size_t> &)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorKokkos, MemSpace>(cutoff),
      _cutoffSquared{static_cast<FloatPrecision>(cutoff * cutoff)}
  {}

  /* Overrides for actual execution */
  void AoSFunctor(Particle_T& i, Particle_T& j, bool newton3) final {

  }

  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    // No-op as nothing should happen here
  }

  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {
    // No-op as nothing should happen here
  }

  void SoAFunctorSingleKokkos(const Particle_T::KokkosSoAArraysType& soa, bool newton3) final {

    const size_t N = soa.size();
    const FloatPrecision cutoffSquared = _cutoffSquared;

    Kokkos::parallel_for(Kokkos::RangePolicy<typename MemSpace::execution_space>(0, N), KOKKOS_LAMBDA(int i) {
        //int i = team.league_rank();

        const auto owned1 = soa.template operator()<Particle_T::AttributeNames::ownershipState, true, false>(i);

        if (owned1 != autopas::OwnershipState::dummy) {
            FloatPrecision fxAcc = 0.;
            FloatPrecision fyAcc = 0.;
            FloatPrecision fzAcc = 0.;

            const auto x1 = soa.template operator()<Particle_T::AttributeNames::posX, true, false>(i);
            const auto y1 = soa.template operator()<Particle_T::AttributeNames::posY, true, false>(i);
            const auto z1 = soa.template operator()<Particle_T::AttributeNames::posZ, true, false>(i);

            for (int j = 0; j < N; ++j) {
                if (i != j) {

                    const auto owned2 = soa.template operator()<Particle_T::AttributeNames::ownershipState, true, false>(j);

                    if (owned2 != autopas::OwnershipState::dummy) {
                        const FloatPrecision x2 = soa.template operator()<Particle_T::AttributeNames::posX, true, false>(j);
                        const FloatPrecision y2 = soa.template operator()<Particle_T::AttributeNames::posY, true, false>(j);
                        const FloatPrecision z2 = soa.template operator()<Particle_T::AttributeNames::posZ, true, false>(j);

                        const FloatPrecision drX = x1 - x2;
                        const FloatPrecision drY = y1 - y2;
                        const FloatPrecision drZ = z1 - z2;

                        const FloatPrecision drX2 = drX * drX;
                        const FloatPrecision drY2 = drY * drY;
                        const FloatPrecision drZ2 = drZ * drZ;

                        const FloatPrecision dr2 = drX2 + drY2 + drZ2;

                        if (dr2 <= cutoffSquared) {
                            // TODO: consider mixing based on type or some sort of parameter injection
                            const FloatPrecision sigmaSquared = 1.;
                            const FloatPrecision epsilon24 = 24.;

                            const FloatPrecision invDr2 = 1.0 / dr2;
                            FloatPrecision lj6 = sigmaSquared * invDr2;
                            lj6 = lj6 * lj6 * lj6;
                            const FloatPrecision lj12 = lj6 * lj6;
                            const FloatPrecision lj12m6 = lj12 - lj6;
                            const FloatPrecision fac = epsilon24 * (lj12 + lj12m6) * invDr2;

                            const FloatPrecision fX = fac * drX;
                            const FloatPrecision fY = fac * drY;
                            const FloatPrecision fZ = fac * drZ;

                            // TODO: consider newton3 if enabled
                            fxAcc += fX;
                            fxAcc += fY;
                            fzAcc += fZ;
                        }
                    }
                }
            }

            const FloatPrecision oldFx = soa.template operator()<Particle_T::AttributeNames::forceX, true, false>(i);
            const FloatPrecision oldFy = soa.template operator()<Particle_T::AttributeNames::forceY, true, false>(i);
            const FloatPrecision oldFz = soa.template operator()<Particle_T::AttributeNames::forceZ, true, false>(i);

            const FloatPrecision newFx = oldFx + fxAcc;
            const FloatPrecision newFy = oldFy + fyAcc;
            const FloatPrecision newFz = oldFz + fzAcc;

            soa.template operator()<Particle_T::AttributeNames::forceX, true, false>(i) = newFx;
            soa.template operator()<Particle_T::AttributeNames::forceY, true, false>(i) = newFy;
            soa.template operator()<Particle_T::AttributeNames::forceZ, true, false>(i) = newFz;
        }
    });
  }

  void SoAFunctorPairKokkos(const Particle_T::KokkosSoAArraysType& soa1, const Particle_T::KokkosSoAArraysType& soa2, bool newton3) final {
    const size_t N = soa1.size();
    const size_t M = soa2.size();

    const FloatPrecision cutoffSquared = _cutoffSquared;

    Kokkos::parallel_for(Kokkos::RangePolicy<typename MemSpace::execution_space>(0, N), KOKKOS_LAMBDA(int i) {

        const auto owned1 = soa1.template operator()<Particle_T::AttributeNames::ownershipState, true, false>(i);

        if (owned1 != autopas::OwnershipState::dummy) {
            FloatPrecision fxAcc = 0.;
            FloatPrecision fyAcc = 0.;
            FloatPrecision fzAcc = 0.;

            const auto x1 = soa1.template operator()<Particle_T::AttributeNames::posX, true, false>(i);
            const auto y1 = soa1.template operator()<Particle_T::AttributeNames::posY, true, false>(i);
            const auto z1 = soa1.template operator()<Particle_T::AttributeNames::posZ, true, false>(i);

            for (int j = 0; j < M; ++j) {
                const auto owned2 = soa2.template operator()<Particle_T::AttributeNames::ownershipState, true, false>(j);

                if (owned2 != autopas::OwnershipState::dummy) {
                    const FloatPrecision x2 = soa2.template operator()<Particle_T::AttributeNames::posX, true, false>(j);
                    const FloatPrecision y2 = soa2.template operator()<Particle_T::AttributeNames::posY, true, false>(j);
                    const FloatPrecision z2 = soa2.template operator()<Particle_T::AttributeNames::posZ, true, false>(j);

                    const FloatPrecision drX = x1 - x2;
                    const FloatPrecision drY = y1 - y2;
                    const FloatPrecision drZ = z1 - z2;

                    const FloatPrecision drX2 = drX * drX;
                    const FloatPrecision drY2 = drY * drY;
                    const FloatPrecision drZ2 = drZ * drZ;

                    const FloatPrecision dr2 = drX2 + drY2 + drZ2;

                    if (dr2 <= cutoffSquared) {
                        // TODO: consider mixing based on type or some sort of parameter injection
                        const FloatPrecision sigmaSquared = 1.;
                        const FloatPrecision epsilon24 = 24.;

                        const FloatPrecision invDr2 = 1.0 / dr2;
                        FloatPrecision lj6 = sigmaSquared * invDr2;
                        lj6 = lj6 * lj6 * lj6;
                        const FloatPrecision lj12 = lj6 * lj6;
                        const FloatPrecision lj12m6 = lj12 - lj6;
                        const FloatPrecision fac = epsilon24 * (lj12 + lj12m6) * invDr2;

                        const FloatPrecision fX = fac * drX;
                        const FloatPrecision fY = fac * drY;
                        const FloatPrecision fZ = fac * drZ;

                        // TODO: consider newton3 if enabled
                        fxAcc += fX;
                        fyAcc += fY;
                        fzAcc += fZ;
                    }
                }
            }


            const FloatPrecision oldFx = soa1.template operator()<Particle_T::AttributeNames::forceX, true, false>(i);
            const FloatPrecision oldFy = soa1.template operator()<Particle_T::AttributeNames::forceY, true, false>(i);
            const FloatPrecision oldFz = soa1.template operator()<Particle_T::AttributeNames::forceZ, true, false>(i);

            const FloatPrecision newFx = oldFx + fxAcc;
            const FloatPrecision newFy = oldFy + fyAcc;
            const FloatPrecision newFz = oldFz + fzAcc;

            soa1.template operator()<Particle_T::AttributeNames::forceX, true, false>(i) = newFx;
            soa1.template operator()<Particle_T::AttributeNames::forceY, true, false>(i) = newFy;
            soa1.template operator()<Particle_T::AttributeNames::forceZ, true, false>(i) = newFz;
        }
    });
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
      Particle_T::AttributeNames::id,
      Particle_T::AttributeNames::posX,
      Particle_T::AttributeNames::posY,
      Particle_T::AttributeNames::posZ,
      Particle_T::AttributeNames::typeId,
      Particle_T::AttributeNames::ownershipState};
  }

  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
      Particle_T::AttributeNames::forceX,
      Particle_T::AttributeNames::forceY,
      Particle_T::AttributeNames::forceZ
  };
  }

  /* Interface required stuff */
  std::string getName() final {
    return "LJFunctorKokkos";
  }

  bool isRelevantForTuning() final {
    return true;
  }

  bool allowsNewton3() final {
    // TODO
    return true;
  }

  bool allowsNonNewton3() final {
    // TODO
    return true;
  }

private:

  FloatPrecision _cutoffSquared;
};
}