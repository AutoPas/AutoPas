/**
 * @file VerletListsKokkosTraversalFlat.h
 * @date 28.05.2026
 * @author Franziska Duhr
 * @note scaffolding taken from DSKokkosTraversalFlat.h
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/verletListsKokkos/traversals/VerletListsKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

#include <Kokkos_Core.hpp>

namespace autopas {

/**
 * This class defines the traversal typically used by the VerletListsKokkos container
 *
 * @tparam Functor
 * @tparam Particle_T
 */
template <class Functor, class Particle_T>
class VerletListsKokkosTraversalFlat : public TraversalInterface, public VerletListsKokkosTraversalInterface<Particle_T> {
    public:
  /**
   * Constructor for the VerletListsKokkosTraversalFlat.
   * @param functor the functor that defines the interaction of particles
   * @param dataLayout The data layout wth which this traversal should be initialized
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not
   */
explicit VerletListsKokkosTraversalFlat(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
        : TraversalInterface(dataLayout, useNewton3), VerletListsKokkosTraversalInterface<Particle_T>(), _functor{functor} {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::vl_kokkos_traversal_flat; }

    [[nodiscard]] bool isApplicable() const final { 
        // TODO
        return true;
    }

    void initTraversal() final {
    }

    void traverseParticles() final {
      const bool newton3 = _useNewton3;
      const auto func = _functor;

      size_t N = VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.size();

      // TODO: this should only be executed on the CPU
      if (_dataLayout == DataLayoutOption::aos) {
#ifdef KOKKOS_ENABLE_CUDA
        return; // TODO: log error / exception
#elif defined(KOKKOS_ENABLE_HIP)
        return; // TODO: log error / exception
#else
        const auto offsets = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListOffsets;
        const auto entries = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListEntries;

        // Full neighbor list: each pair appears under both partners, so the gather-only
        // update of particle i is correct without newton3. Force newton3 off here.
        Kokkos::parallel_for("traverseParticlesAoS", Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, N), KOKKOS_LAMBDA(int i)  {
          for (size_t k = offsets(i); k < offsets(i + 1); ++k) {
            const size_t j = entries(k);
            func->AoSFunctorKokkos(
              VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.getAoS().getParticle(i),
              VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.getAoS().getParticle(j),
              false);
          }
        });

        // TODO: consider halo particles
        // Maybe even execute halo traversal simultaneously on the host with some sort of result merge mechanism in the end
#endif
      }
      else if (_dataLayout == DataLayoutOption::soa) {
        auto& ownedSoA = VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.getSoA();
        auto& haloSoA = VerletListsKokkosTraversalInterface<Particle_T>::_haloParticles.getSoA();

        const auto I = std::make_index_sequence<Functor::getNeededAttr().size()>{};

        syncNeeded<DeviceSpace::execution_space>(ownedSoA, I);
        syncNeeded<DeviceSpace::execution_space>(haloSoA, I);

        performSoATraversal(ownedSoA);
        // TODO: owned-halo interactions need their own (owned-halo) neighbor list before they
        // can be traversed here.

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

  // Owned-owned traversal over the (full) neighbor list: particle i only interacts with the
  // entries stored for it, which is the whole point of the Verlet list traversal.
  void performSoATraversal(const Particle_T::KokkosSoAArraysType& soa) {
    const size_t N = soa.size();

    if (N == 0) {
      return;
    }

    FloatPrecision cutoffSquared = _functor->getCutoff() * _functor->getCutoff();
    const auto soaDevice = soa.deviceView();
    const auto offsets = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListOffsets;
    const auto entries = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListEntries;

    auto rangePolicy = Kokkos::RangePolicy<typename DeviceSpace::execution_space>(0, N);
    Kokkos::parallel_for("vl_kokkos_flat", rangePolicy, KOKKOS_LAMBDA(const int i) {
      FloatPrecision fxAcc = 0.;
      FloatPrecision fyAcc = 0.;
      FloatPrecision fzAcc = 0.;

      const auto x1 = soaDevice.template operator()<Particle_T::AttributeNames::posX, true>(i);
      const auto y1 = soaDevice.template operator()<Particle_T::AttributeNames::posY, true>(i);
      const auto z1 = soaDevice.template operator()<Particle_T::AttributeNames::posZ, true>(i);

      for (size_t k = offsets(i); k < offsets(i + 1); ++k) {
        const size_t j = entries(k);
        Functor::SoAKernelKokkosStatic(x1, y1, z1, soaDevice, fxAcc, fyAcc, fzAcc, cutoffSquared, i, j);
      }

      soaDevice.template operator()<Particle_T::AttributeNames::forceX, true>(i) += fxAcc;
      soaDevice.template operator()<Particle_T::AttributeNames::forceY, true>(i) += fyAcc;
      soaDevice.template operator()<Particle_T::AttributeNames::forceZ, true>(i) += fzAcc;
    });
  }


#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ENABLE_HIP)
  using DeviceSpace = Kokkos::HIPSpace; 
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif

  using HostSpace = Kokkos::HostSpace;

  using FloatPrecision = Particle_T::ParticleSoAFloatPrecision;

  Functor *_functor;
};
} 
