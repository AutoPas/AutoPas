/**
 * @file KokkosDsNaiveParallelTraversal.h
 * @date 31. October 2025
 * @author Luis Gall
 */

#pragma once

#include <Kokkos_Core.hpp>

#include "../../../../../cmake-build-debug/_deps/kokkos-src/containers/src/Kokkos_ScatterView.hpp"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/kokkosDirectSum/traversals/DSKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

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

    const size_t chunkSize = _chunkSize;
    const size_t numChunks = N / chunkSize;

    auto teamPolicy = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>(numChunks+1, _teamSize, Kokkos::AUTO);

    using ScratchView = Kokkos::View<FloatPrecision*, typename DeviceSpace::execution_space::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Atomic>>;
    int size = 3 * ScratchView::shmem_size(chunkSize);
    teamPolicy.set_scratch_size(0, Kokkos::PerTeam(size));

    using MemberType = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>::member_type;
    Kokkos::parallel_for("traversal", teamPolicy, KOKKOS_LAMBDA(const MemberType& teamHandle) {
      const int k = teamHandle.league_rank();

      size_t offset = k* chunkSize;
      size_t rest = N - offset;
      size_t upper = std::min(rest, chunkSize);

      ScratchView fX (teamHandle.team_scratch((0)), upper);
      ScratchView fY (teamHandle.team_scratch((0)), upper);
      ScratchView fZ (teamHandle.team_scratch((0)), upper);

      Kokkos::parallel_for(Kokkos::TeamThreadMDRange(teamHandle, upper, M), [&](int i, int j) {
        const auto x1 = soa1.template operator()<Particle_T::AttributeNames::posX, true, false>(i+offset);
        const auto y1 = soa1.template operator()<Particle_T::AttributeNames::posY, true, false>(i+offset);
        const auto z1 = soa1.template operator()<Particle_T::AttributeNames::posZ, true, false>(i+offset);

        FloatPrecision localFxAcc = 0;
        FloatPrecision localFyAcc = 0;
        FloatPrecision localFzAcc = 0;
        func->SoAKernelKokkos(x1, y1, z1, soa2, localFxAcc, localFyAcc, localFzAcc, cutoffSquared, i, j);

        fX(i) += localFxAcc;
        fY(i) += localFyAcc;
        fZ(i) += localFzAcc;
      });

      teamHandle.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamHandle, chunkSize), [&](int i) {
        int index = i + offset;
        const FloatPrecision oldFx = soa1.template operator()<Particle_T::AttributeNames::forceX, true, false>(index);
        const FloatPrecision oldFy = soa1.template operator()<Particle_T::AttributeNames::forceY, true, false>(index);
        const FloatPrecision oldFz = soa1.template operator()<Particle_T::AttributeNames::forceZ, true, false>(index);

        const FloatPrecision newFx = oldFx + fX(i);
        const FloatPrecision newFy = oldFy + fY(i);
        const FloatPrecision newFz = oldFz + fZ(i);

        soa1.template operator()<Particle_T::AttributeNames::forceX, true, false>(index) = newFx;
        soa1.template operator()<Particle_T::AttributeNames::forceY, true, false>(index) = newFy;
        soa1.template operator()<Particle_T::AttributeNames::forceZ, true, false>(index) = newFz;
      });
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