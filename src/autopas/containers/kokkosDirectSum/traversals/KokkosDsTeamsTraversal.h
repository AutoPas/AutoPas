/**
 * @file KokkosDsTeamsTraversal.h
 * @date 31. October 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/containers/kokkosDirectSum/traversals/DSKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/baseFunctors/KokkosFunctor.h"

namespace autopas {

template <class Functor, class Particle_T, class DeviceSpace>
struct TeamTraversalFunctor {
    using FloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;
    using MemberType = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>::member_type;

    utilsKokkos::KokkosStorage<Particle_T> _storageA;
    utilsKokkos::KokkosStorage<Particle_T> _storageB;
    Functor* _func;
    FloatPrecision _cutoffSquared;
    size_t _M;

    struct ReductionResult {
      FloatPrecision virialSum;
      FloatPrecision uPotSum;
    };

    KOKKOS_INLINE_FUNCTION
    void init(ReductionResult& result) const {
      result.virialSum = 0.;
      result.uPotSum = 0.;
    }

    KOKKOS_INLINE_FUNCTION
    void join(ReductionResult& result, const ReductionResult& src) const {
      result.virialSum += src.virialSum;
      result.uPotSum += src.uPotSum;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType& teamHandle, ReductionResult& localResult) const {
        const int i = teamHandle.league_rank();

        FloatPrecision fxAcc = 0.;
        FloatPrecision fyAcc = 0.;
        FloatPrecision fzAcc = 0.;

        const auto x1 = _storageA.template operator()<Particle_T::AttributeNames::posX, false>(i);
        const auto y1 = _storageA.template operator()<Particle_T::AttributeNames::posY, false>(i);
        const auto z1 = _storageA.template operator()<Particle_T::AttributeNames::posZ, false>(i);

        FloatPrecision virialSum = 0.;
        FloatPrecision uPotSum = 0;

        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamHandle, _M), [&](int j,
            FloatPrecision& localFxAcc,
            FloatPrecision& localFyAcc,
            FloatPrecision& localFzAcc,
            FloatPrecision& localVirialSum,
            FloatPrecision& localUPotSum) {
                _func->ForceKernelKokkos(x1, y1, z1, _storageB, localFxAcc, localFyAcc, localFzAcc, localVirialSum, localUPotSum, _cutoffSquared, i, j);
            }, fxAcc, fyAcc, fzAcc, virialSum, uPotSum);

        Kokkos::single(Kokkos::PerTeam(teamHandle), [&]() {
          localResult.virialSum += virialSum;
          localResult.uPotSum += uPotSum;

          _storageA.template operator()<Particle_T::AttributeNames::forceX, false>(i) += fxAcc;
          _storageA.template operator()<Particle_T::AttributeNames::forceY, false>(i) += fyAcc;
          _storageA.template operator()<Particle_T::AttributeNames::forceZ, false>(i) += fzAcc;
        });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType& teamHandle) const {
        const int i = teamHandle.league_rank();

        FloatPrecision fxAcc = 0.;
        FloatPrecision fyAcc = 0.;
        FloatPrecision fzAcc = 0.;

        const auto x1 = _storageA.template operator()<Particle_T::AttributeNames::posX, false>(i);
        const auto y1 = _storageA.template operator()<Particle_T::AttributeNames::posY, false>(i);
        const auto z1 = _storageA.template operator()<Particle_T::AttributeNames::posZ, false>(i);

        FloatPrecision virialSum = 0.;
        FloatPrecision uPotSum = 0;

        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamHandle, _M), [&](int j,
            FloatPrecision& localFxAcc,
            FloatPrecision& localFyAcc,
            FloatPrecision& localFzAcc) {
                _func->ForceKernelKokkos(x1, y1, z1, _storageB, localFxAcc, localFyAcc, localFzAcc, virialSum, uPotSum, _cutoffSquared, i, j);
            }, fxAcc, fyAcc, fzAcc);

        Kokkos::single(Kokkos::PerTeam(teamHandle), [&]() {
          const int index = teamHandle.league_rank();
          const FloatPrecision oldFx = _storageA.template operator()<Particle_T::AttributeNames::forceX, false>(index);
          const FloatPrecision oldFy = _storageA.template operator()<Particle_T::AttributeNames::forceY, false>(index);
          const FloatPrecision oldFz = _storageA.template operator()<Particle_T::AttributeNames::forceZ, false>(index);
        });
    }
};

/**
 * This class defines the traversal typically used by the KokkosDirectSumContainer
 *
 * @tparam Functor
 * @tparam Particle_T
 */
template <class Functor, class Particle_T>
class KokkosDsTeamsTraversal : public DSKokkosTraversalInterface<Particle_T> {

public:
  /**
   * Constructor for the KokkkosDSNaiveParallelTraversal
   * @param functor the functor that defines the interaction of particles
   * @param dataLayout The data layout wth which this traversal should be initialized
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not
   * @param teamSize Size of the Kokkos Teams
   */
explicit KokkosDsTeamsTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3, size_t teamSize)
        : DSKokkosTraversalInterface<Particle_T>(dataLayout, useNewton3), _functor{functor}, _teamSize(teamSize) {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::ds_kokkos_teams; }

    [[nodiscard]] bool isApplicable() const final { 
        // TODO
        return true;
    }

  void initTraversal() final {
  const auto I = std::make_index_sequence<Functor::getNeededAttr().size()>{};

  syncNeeded<typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace::execution_space>(DSKokkosTraversalInterface<Particle_T>::_ownedParticles, I);
  syncNeeded<typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace::execution_space>(DSKokkosTraversalInterface<Particle_T>::_haloParticles, I);
}

  void endTraversal() final {
  constexpr auto J = std::make_index_sequence<Functor::getComputedAttr().size()>{};
  modifyComputed<typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace::execution_space>(DSKokkosTraversalInterface<Particle_T>::_ownedParticles, J);
}

protected:
  void performTraversal(const utilsKokkos::KokkosStorage<Particle_T>& storageA, const utilsKokkos::KokkosStorage<Particle_T>& storageB) final {
    const size_t N = storageA.size();
    const size_t M = storageB.size();

    if (N == 0 || M == 0) {
      return;
    }

    auto func = _functor;
    typename DSKokkosTraversalInterface<Particle_T>::FloatPrecision cutoffSquared = func->getCutoff() * func->getCutoff();

    TeamTraversalFunctor<Functor, Particle_T, typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace> functor {storageA, storageB, func, cutoffSquared, M};
    auto teamPolicy = Kokkos::TeamPolicy<typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace::execution_space>(N, _teamSize, Kokkos::AUTO);

    constexpr bool calculateGlobals = Functor::globalCalculationRequested();

    if constexpr (calculateGlobals) {
      typename TeamTraversalFunctor<Functor, Particle_T, typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace>::ReductionResult globalResult {};
      Kokkos::parallel_reduce("autopas::KokkosDsTeamsTraversal_Globals", teamPolicy, functor, globalResult);

      AutoPasLog(INFO, "Final potential energy {}", globalResult.uPotSum / 12.);
      AutoPasLog(INFO, "Final virial           {}", globalResult.virialSum / 2.);

      auto kokkosFunc = dynamic_cast<KokkosFunctor*>(func);

      kokkosFunc->setPotentialEnergy(static_cast<double>(globalResult.uPotSum / 12.));
      kokkosFunc->setVirial(static_cast<double>(globalResult.virialSum / 2.));
    } else {
      Kokkos::parallel_for("autopas::KokkosDsTeamsTraversal_NoGlobals", teamPolicy, functor);
    }
  }

private:

  template <typename ExecSpace, std::size_t... I>
  void syncNeeded(auto& particles, std::index_sequence<I...>) {
    (particles.template sync<ExecSpace, Functor::getNeededAttr()[I]-1>(), ...);
  }

  template <typename ExecSpace, std::size_t... I>
  void modifyComputed(auto& particles, std::index_sequence<I...>) {
    (particles.template modify<ExecSpace, Functor::getComputedAttr()[I]-1>(), ...);
  }

  Functor* _functor;

  const size_t _teamSize {0};
};

}