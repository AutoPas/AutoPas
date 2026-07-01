/**
 * @file KokkosDsFlatTraversal.h
 * @date 31. October 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/containers/kokkosDirectSum/traversals/DSKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/baseFunctors/KokkosFunctor.h"

#include <Kokkos_Core.hpp>

namespace autopas {

template <class Functor, class Particle_T>
struct FlatTraversalReductionFunctor {
  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  utilsKokkos::KokkosStorage<Particle_T> _storageA;
  utilsKokkos::KokkosStorage<Particle_T> _storageB;
  Functor* _func;
  FloatPrecision _cutoffSquared;
  size_t M;

  struct ReductionResult {
    FloatPrecision virialSum;
    FloatPrecision uPotSum;
  };

  using value_type = ReductionResult;

  KOKKOS_INLINE_FUNCTION
  void init(value_type& result) const {
    result.virialSum = 0.;
    result.uPotSum = 0.;
  }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& result, const value_type& src) const {
    result.virialSum += src.virialSum;
    result.uPotSum += src.uPotSum;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type& localResult) const {
    FloatPrecision fxAcc = 0.;
    FloatPrecision fyAcc = 0.;
    FloatPrecision fzAcc = 0.;

    const auto x1 = _storageA.template operator()<Particle_T::AttributeNames::posX, false>(i);
    const auto y1 = _storageA.template operator()<Particle_T::AttributeNames::posY, false>(i);
    const auto z1 = _storageA.template operator()<Particle_T::AttributeNames::posZ, false>(i);

    FloatPrecision virialSum = 0.;
    FloatPrecision uPotSum = 0;

    for (int j = 0; j < M; ++j) {
      _func->ForceKernelKokkos(x1, y1, z1, _storageB, fxAcc, fyAcc, fzAcc, virialSum, uPotSum, _cutoffSquared, i, j);
    }

    _storageA.template operator()<Particle_T::AttributeNames::forceX, false>(i) += fxAcc;
    _storageA.template operator()<Particle_T::AttributeNames::forceY, false>(i) += fyAcc;
    _storageA.template operator()<Particle_T::AttributeNames::forceZ, false>(i) += fzAcc;

    localResult.virialSum += virialSum;
    localResult.uPotSum += uPotSum;
  }
};

template <class Functor, class Particle_T>
struct FlatTraversalFunctor {
  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  utilsKokkos::KokkosStorage<Particle_T> _storageA;
  utilsKokkos::KokkosStorage<Particle_T> _storageB;
  Functor* _func;
  FloatPrecision _cutoffSquared;
  size_t M;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    FloatPrecision fxAcc = 0.;
    FloatPrecision fyAcc = 0.;
    FloatPrecision fzAcc = 0.;

    FloatPrecision virialSum = 0.;
    FloatPrecision uPotSum = 0;

    const auto x1 = _storageA.template operator()<Particle_T::AttributeNames::posX, false>(i);
    const auto y1 = _storageA.template operator()<Particle_T::AttributeNames::posY, false>(i);
    const auto z1 = _storageA.template operator()<Particle_T::AttributeNames::posZ, false>(i);

    for (int j = 0; j < M; ++j) {
      _func->ForceKernelKokkos(x1, y1, z1, _storageB, fxAcc, fyAcc, fzAcc, virialSum, uPotSum, _cutoffSquared, i, j);
    }

    _storageA.template operator()<Particle_T::AttributeNames::forceX, false>(i) += fxAcc;
    _storageA.template operator()<Particle_T::AttributeNames::forceY, false>(i) += fyAcc;
    _storageA.template operator()<Particle_T::AttributeNames::forceZ, false>(i) += fzAcc;
  }
};

/**
 * This class defines the traversal typically used by the KokkosDirectSumContainer
 *
 * @tparam Functor
 * @tparam Particle_T
 */
template <class Functor, class Particle_T>
class KokkosDsFlatTraversal : public DSKokkosTraversalInterface<Particle_T> {

public:
  /**
   * Constructor for the KokkkosDSNaiveParallelTraversal
   * @param functor the functor that defines the interaction of particles
   * @param dataLayout The data layout wth which this traversal should be initialized
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not
   */
explicit KokkosDsFlatTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
        : DSKokkosTraversalInterface<Particle_T>(dataLayout, useNewton3), _functor{functor} {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::ds_kokkos_flat; }

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

    auto rangePolicy = Kokkos::RangePolicy<typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace::execution_space>(0, N);

    constexpr bool calculateGlobals = Functor::globalCalculationRequested();
    if constexpr (calculateGlobals) {
      typename FlatTraversalReductionFunctor<Functor, Particle_T>::value_type globalResult {};
      FlatTraversalReductionFunctor<Functor, Particle_T> reductionFunctor {storageA, storageB, func, cutoffSquared, M};
      Kokkos::parallel_reduce("autopas::KokkosDsFlatTraversal_Globals", rangePolicy, reductionFunctor, globalResult);

      AutoPasLog(INFO, "Final potential energy {}", static_cast<double>(globalResult.uPotSum) / 12.);
      AutoPasLog(INFO, "Final virial           {}", static_cast<double>(globalResult.virialSum) * 0.5);

      auto kokkosFunc = dynamic_cast<KokkosFunctor*>(func);

      kokkosFunc->setPotentialEnergy(static_cast<double>(globalResult.uPotSum) / 12.);
      kokkosFunc->setVirial(static_cast<double>(globalResult.virialSum) * 0.5);

    } else {
      FlatTraversalFunctor<Functor, Particle_T> functor {storageA, storageB, func, cutoffSquared, M};
      Kokkos::parallel_for("autopas::KokkosDsFlatTraversal_NoGlobals", rangePolicy, functor);
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
};

}