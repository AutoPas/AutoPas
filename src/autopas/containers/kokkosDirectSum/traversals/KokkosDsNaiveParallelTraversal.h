/**
 * @file KokkosDsNaiveParallelTraversal.h
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
class KokkosDsNaiveParallelTraversal : public TraversalInterface, public DSKokkosTraversalInterface<Particle_T> {

public:

    explicit KokkosDsNaiveParallelTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3, size_t teamSize, size_t chunkSize)
        : TraversalInterface(dataLayout, useNewton3), DSKokkosTraversalInterface<Particle_T>(), _functor{functor}, _teamSize(teamSize), _chunkSize(chunkSize) {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::kokkos_ds_naive_parallel; }

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

    //auto teamPolicy = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>(N, Kokkos::AUTO, Kokkos::AUTO);
    //using MemberType = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>::member_type;

    auto rangePolicy = Kokkos::RangePolicy<typename DeviceSpace::execution_space>(0, N, Kokkos::ChunkSize(_chunkSize));
    // Kokkos::parallel_for("traversal", teamPolicy, KOKKOS_LAMBDA(const MemberType& teamHandle) {
    // const int i = teamHandle.league_rank();

    Kokkos::parallel_for("traversal", rangePolicy, KOKKOS_LAMBDA(const int i) {
      FloatPrecision fxAcc = 0.;
      FloatPrecision fyAcc = 0.;
      FloatPrecision fzAcc = 0.;

      const auto x1 = soa1.template operator()<Particle_T::AttributeNames::posX, true, false>(i);
      const auto y1 = soa1.template operator()<Particle_T::AttributeNames::posY, true, false>(i);
      const auto z1 = soa1.template operator()<Particle_T::AttributeNames::posZ, true, false>(i);

      /*
      // Kokkos::parallel_reduce(Kokkos::TeamVectorRange(teamHandle, M), [&](int j,
          FloatPrecision& localFxAcc,
          FloatPrecision& localFyAcc,
          FloatPrecision& localFzAcc) {
          */
      for (int j = 0; j < M; ++j) {
        func->SoAKernelKokkos(x1, y1, z1, soa2, fxAcc, fyAcc, fzAcc, cutoffSquared, i, j);
      }
        //}, fxAcc, fyAcc, fzAcc);

        //Kokkos::single(Kokkos::PerTeam(teamHandle), [&]() {
          // int index = teamHandle.league_rank();
          soa1.template operator()<Particle_T::AttributeNames::forceX, true, false>(i) += fxAcc;
          soa1.template operator()<Particle_T::AttributeNames::forceY, true, false>(i) += fyAcc;
          soa1.template operator()<Particle_T::AttributeNames::forceZ, true, false>(i) += fzAcc;
        // });
    });
  }


#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif

  using HostSpace = Kokkos::HostSpace;

  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  Functor *_functor;

  const size_t _teamSize {0};

  const size_t _chunkSize {0};
};

}