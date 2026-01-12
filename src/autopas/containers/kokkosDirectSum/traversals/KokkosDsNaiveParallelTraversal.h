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

    explicit KokkosDsNaiveParallelTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
        : TraversalInterface(dataLayout, useNewton3), DSKokkosTraversalInterface<Particle_T>(), _functor{functor} {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::kokkos_ds_naive_parallel; }

    [[nodiscard]] bool isApplicable() const final { 
        // TODO
        return true;
    }

    void initTraversal() final {
    }

    void traverseParticles() final {
      // TODO: be aware of data layout
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
        auto& haloSoA = DSKokkosTraversalInterface<Particle_T>::_haloParticles.getSoA().getSoA();

        constexpr auto tupleSize = ownedSoA.tupleSize();
        constexpr auto I = std::make_index_sequence<tupleSize>();

        ownedSoA.template sync<DeviceSpace::execution_space>(I);
        haloSoA.template sync<DeviceSpace::execution_space>(I);

        func->SoAFunctorSingleKokkos(ownedSoA, newton3);
        func->SoAFunctorPairKokkos(ownedSoA, haloSoA, newton3);

        // TODO: only modify changed attributes (Functor will have to return which attributes he did change)
        ownedSoA.template markModified<DeviceSpace::execution_space>(I);
        // TODO: consider halo particles
      }
    }

    void endTraversal() final {
    }

private:


#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif

  using HostSpace = Kokkos::HostSpace;

  Functor *_functor;

};

}