/**
 * @file KokkosDsFlatTraversal.h
 * @date 31. October 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/kokkosDirectSum/traversals/DSKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

#include <Kokkos_Core.hpp>

namespace autopas {

template <class Functor, class Particle_T>
struct SoAFlatTraversalFunctor {
  using KokkosSoAArraysType = typename Particle_T::KokkosSoAArraysType;
  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  KokkosSoAArraysType _soa1;
  KokkosSoAArraysType _soa2;
  Functor* _func;
  FloatPrecision _cutoffSquared;
  size_t M;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    FloatPrecision fxAcc = 0.;
    FloatPrecision fyAcc = 0.;
    FloatPrecision fzAcc = 0.;

    const auto x1 = _soa1.template operator()<Particle_T::AttributeNames::posX, false>(i);
    const auto y1 = _soa1.template operator()<Particle_T::AttributeNames::posY, false>(i);
    const auto z1 = _soa1.template operator()<Particle_T::AttributeNames::posZ, false>(i);

    for (int j = 0; j < M; ++j) {
      _func->SoAKernelKokkos(x1, y1, z1, _soa2, fxAcc, fyAcc, fzAcc, _cutoffSquared, i, j);
    }

    _soa1.template operator()<Particle_T::AttributeNames::forceX, false>(i) += fxAcc;
    _soa1.template operator()<Particle_T::AttributeNames::forceY, false>(i) += fyAcc;
    _soa1.template operator()<Particle_T::AttributeNames::forceZ, false>(i) += fzAcc;
  }
};

/**
 * This class defines the traversal typically used by the KokkosDirectSumContainer
 *
 * @tparam Functor
 * @tparam Particle_T
 */
template <class Functor, class Particle_T>
class KokkosDsFlatTraversal : public TraversalInterface, public DSKokkosTraversalInterface<Particle_T> {

public:
  /**
   * Constructor for the KokkkosDSNaiveParallelTraversal
   * @param functor the functor that defines the interaction of particles
   * @param dataLayout The data layout wth which this traversal should be initialized
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not
   */
explicit KokkosDsFlatTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
        : TraversalInterface(dataLayout, useNewton3), DSKokkosTraversalInterface<Particle_T>(), _functor{functor} {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::ds_kokkos_flat; }

    [[nodiscard]] bool isApplicable() const final { 
        // TODO
        return true;
    }

    void initTraversal() final {
    }

    void traverseParticles() final {
      const bool newton3 = _useNewton3;
      const auto func = _functor;

      size_t N = DSKokkosTraversalInterface<Particle_T>::_ownedParticles.size();

      // TODO: this should only be executed on the CPU
      if (_dataLayout == DataLayoutOption::aos) {
#ifdef KOKKOS_ENABLE_CUDA
        return; // TODO: log error / exception
#else
        Kokkos::parallel_for("traverseParticlesAoS", Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, N), KOKKOS_LAMBDA(int i)  {
          for (int j = (newton3 ? i+1 : 0); j < N; ++j) {
            if (newton3 or i != j) {
              func->AoSFunctorKokkos(
                DSKokkosTraversalInterface<Particle_T>::_ownedParticles.getAoS().getParticle(i),
                DSKokkosTraversalInterface<Particle_T>::_ownedParticles.getAoS().getParticle(j),
                newton3);
            }
          }
        });

        // TODO: consider halo particles
        // Maybe even execute halo traversal simultaneously on the host with some sort of result merge mechanism in the end
#endif
      }
      else if (_dataLayout == DataLayoutOption::soa) {

        auto& ownedSoA = DSKokkosTraversalInterface<Particle_T>::_ownedParticles.getSoA();
        auto& haloSoA = DSKokkosTraversalInterface<Particle_T>::_haloParticles.getSoA();

        const auto I = std::make_index_sequence<Functor::getNeededAttr().size()>{};

        syncNeeded<DeviceSpace::execution_space>(ownedSoA, I);
        syncNeeded<DeviceSpace::execution_space>(haloSoA, I);

        performSoATraversal(ownedSoA, ownedSoA);
        performSoATraversal(ownedSoA, haloSoA);

        constexpr auto J = std::make_index_sequence<Functor::getComputedAttr().size()>{};
        modifyComputed<DeviceSpace::execution_space>(ownedSoA, J);
      }
    }

    void endTraversal() final {
    }

private:

  template <typename ExecSpace, std::size_t... I>
  void syncNeeded(auto& particles, std::index_sequence<I...>) {
    (particles.template sync<ExecSpace, Functor::getNeededAttr()[I]-1>(), ...);
  }

  template <typename ExecSpace, std::size_t... I>
  void modifyComputed(auto& particles, std::index_sequence<I...>) {
    (particles.template markModified<ExecSpace, Functor::getComputedAttr()[I]-1>(), ...);
  }

  void performSoATraversal(const Particle_T::KokkosSoAArraysType& soa1, const Particle_T::KokkosSoAArraysType& soa2) {
    const size_t N = soa1.size();
    const size_t M = soa2.size();

    if (N == 0 || M == 0) {
      return;
    }

    auto func = _functor;
    FloatPrecision cutoffSquared = func->getCutoff() * func->getCutoff();

    SoAFlatTraversalFunctor<Functor, Particle_T> functor {soa1, soa2, func, cutoffSquared, M};

    auto rangePolicy = Kokkos::RangePolicy<typename DeviceSpace::execution_space>(0, N);
    Kokkos::Timer timer {};
    Kokkos::parallel_for("autopas::FlatTraversalSoA", rangePolicy, functor);
    std::cout << timer.seconds() << std::endl;
  }


#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif

  using HostSpace = Kokkos::HostSpace;

  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  Functor *_functor;
};

}