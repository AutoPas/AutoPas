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
      auto& owned = DSKokkosTraversalInterface<Particle_T>::_ownedParticlesAoS;
      const bool newton3 = _useNewton3;

      const auto func = _functor;
      const int N = owned.size();

      // TODO: think about MD Range policy to allow parallelism for both loops
      Kokkos::parallel_for("traverseParticles", Kokkos::RangePolicy<DeviceSpace::execution_space>(0, N), KOKKOS_LAMBDA(int i) {
        for (int j = (newton3 ? i+1 : 0); j < N; ++j) {
          if (newton3 or i != j) {
            func->AoSFunctor(owned(i), owned(j), newton3);
          }
        }
      });

      // TODO: consider halo particles
      // Maybe even execute halo traversal simultaneously on the host with some sort of data storage mechanism
    }

    void endTraversal() final {
    }

private:

  Functor *_functor;

};

}