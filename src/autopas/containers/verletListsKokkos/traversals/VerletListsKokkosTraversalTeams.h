/**
 * @file VerletListsKokkosTraversalTeams.h
 * @date 01.06.2026
 * @author Franziska Duhr
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/verletListsKokkos/traversals/VerletListsKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/logging/Logger.h"

#include <Kokkos_Core.hpp>

namespace autopas {

/**
 * Traversal for VerletListsKokkos that uses Kokkos::TeamPolicy for intra-particle parallelism.
 *
 * @tparam Functor
 * @tparam Particle_T
 */
template <class Functor, class Particle_T>
class VerletListsKokkosTraversalTeams : public TraversalInterface, public VerletListsKokkosTraversalInterface<Particle_T> {

public:
  /**
   * Constructor for the VerletListsKokkosTraversalTeams
   * @param functor the functor that defines the interaction of particles
   * @param dataLayout The data layout wth which this traversal should be initialized
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not
   * @param teamSize Size of the Kokkos Teams
   */
explicit VerletListsKokkosTraversalTeams(Functor *functor, DataLayoutOption dataLayout, bool useNewton3, size_t teamSize)
        : TraversalInterface(dataLayout, useNewton3), VerletListsKokkosTraversalInterface<Particle_T>(), _functor{functor}, _teamSize(teamSize) {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::vl_kokkos_traversal_teams; }

    [[nodiscard]] bool isApplicable() const final { 
        // TODO
        return true;
    }

    void initTraversal() final {
    }

  void traverseParticles() final {
      [[maybe_unused]] const bool newton3 = _useNewton3;  // kept for a future newton3 implementation
      const auto func = _functor;

      size_t N = VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.size();

      spdlog::debug("VerletListsKokkosTraversalTeams::traverseParticles: dataLayout={}, ownedParticles={}",
                    _dataLayout.to_string(), N);

      // TODO: this should only be executed on the CPU
      if (_dataLayout == DataLayoutOption::aos) {
    #ifdef KOKKOS_ENABLE_CUDA
        spdlog::warn(
            "VerletListsKokkosTraversalTeams: AoS data layout is not supported on the CUDA backend - no forces "
            "were computed. Use SoA.");
        return;
    #elif defined(KOKKOS_ENABLE_HIP)
        spdlog::warn(
            "VerletListsKokkosTraversalTeams: AoS data layout is not supported on the HIP backend - no forces "
            "were computed. Use SoA.");
        return;
    #else
        const auto offsets = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListOffsets;
        const auto entries = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListEntries;
        const size_t stride = VerletListsKokkosTraversalInterface<Particle_T>::_neighborListStride;
        const auto haloOffsets = VerletListsKokkosTraversalInterface<Particle_T>::_haloNeighborListOffsets;
        const auto haloEntries = VerletListsKokkosTraversalInterface<Particle_T>::_haloNeighborListEntries;
        const size_t haloStride = VerletListsKokkosTraversalInterface<Particle_T>::_haloNeighborListStride;

        // Full neighbor list: each pair appears under both partners, so the gather-only
        // update of particle i is correct without newton3
        Kokkos::parallel_for("traverseParticlesAoS", Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(0, N), KOKKOS_LAMBDA(int i)  {
          auto& ownedI = VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.getAoS().getParticle(i);
          for (size_t k = (stride ? static_cast<size_t>(i) * stride : offsets(i)),
                      kEnd = (stride ? offsets(i) : offsets(i + 1)); k < kEnd; ++k) {
            const size_t j = entries(k);
            func->AoSFunctorKokkos(
              ownedI,
              VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.getAoS().getParticle(j),
              false);
          }
          for (size_t k = (haloStride ? static_cast<size_t>(i) * haloStride : haloOffsets(i)),
                      kEnd = (haloStride ? haloOffsets(i) : haloOffsets(i + 1)); k < kEnd; ++k) {
            const size_t j = haloEntries(k);
            func->AoSFunctorKokkos(
              ownedI,
              VerletListsKokkosTraversalInterface<Particle_T>::_haloParticles.getAoS().getParticle(j),
              false);
          }
        });
    #endif
      }
      else if (_dataLayout == DataLayoutOption::soa) {
        auto& ownedSoA = VerletListsKokkosTraversalInterface<Particle_T>::_ownedParticles.getSoA();
        auto& haloSoA = VerletListsKokkosTraversalInterface<Particle_T>::_haloParticles.getSoA();

        const auto I = std::make_index_sequence<Functor::getNeededAttr().size()>{};

        syncNeeded<DeviceSpace::execution_space>(ownedSoA, I);
        syncNeeded<DeviceSpace::execution_space>(haloSoA, I);

        performSoATraversal(ownedSoA, ownedSoA,
                            VerletListsKokkosTraversalInterface<Particle_T>::_neighborListOffsets,
                            VerletListsKokkosTraversalInterface<Particle_T>::_neighborListEntries,
                            VerletListsKokkosTraversalInterface<Particle_T>::_neighborListStride);
        performSoATraversal(ownedSoA, haloSoA,
                            VerletListsKokkosTraversalInterface<Particle_T>::_haloNeighborListOffsets,
                            VerletListsKokkosTraversalInterface<Particle_T>::_haloNeighborListEntries,
                            VerletListsKokkosTraversalInterface<Particle_T>::_haloNeighborListStride);

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

  // Traverse a neighbor list: each owned particle i (in soa1) gathers force from the partner
  // particles (in soa2) stored for it in the list. One Kokkos team handles one particle i, with
  // the team threads sharing the work over i's neighbor list.
  void performSoATraversal(const Particle_T::KokkosSoAArraysType& soa1, const Particle_T::KokkosSoAArraysType& soa2,
                           const Kokkos::View<size_t*>& offsets, const Kokkos::View<size_t*>& entries, size_t stride) {
    const size_t N = soa1.size();

    spdlog::debug("VerletListsKokkosTraversalTeams::performSoATraversal: soa1.size()={}, soa2.size()={}", N,
                  soa2.size());

    if (N == 0 || soa2.size() == 0) {
      spdlog::debug(
          "VerletListsKokkosTraversalTeams::performSoATraversal: skipping kernel launch (soa1.size()={}, "
          "soa2.size()={})",
          N, soa2.size());
      return;
    }

    FloatPrecision cutoffSquared = _functor->getCutoff() * _functor->getCutoff();
    const auto soa1Device = soa1.deviceView();
    const auto soa2Device = soa2.deviceView();

    auto teamPolicy = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>(N, _teamSize, Kokkos::AUTO);
    using MemberType = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>::member_type;
    Kokkos::parallel_for("traversal", teamPolicy, KOKKOS_LAMBDA(const MemberType& teamHandle) {
      const int i = teamHandle.league_rank();

      const auto x1 = soa1Device.template operator()<Particle_T::AttributeNames::posX, true>(i);
      const auto y1 = soa1Device.template operator()<Particle_T::AttributeNames::posY, true>(i);
      const auto z1 = soa1Device.template operator()<Particle_T::AttributeNames::posZ, true>(i);

      const size_t listBegin = (stride ? static_cast<size_t>(i) * stride : offsets(i));
      const size_t listEnd = (stride ? offsets(i) : offsets(i + 1));
      const size_t numNeighbors = listEnd - listBegin;

      // The team cooperatively reduces particle i's neighbor list into a single force contribution.
      FloatPrecision fxAcc = 0.;
      FloatPrecision fyAcc = 0.;
      FloatPrecision fzAcc = 0.;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(teamHandle, numNeighbors),
          [&](const size_t k, FloatPrecision& fx, FloatPrecision& fy, FloatPrecision& fz) {
            const size_t j = entries(listBegin + k);
            Functor::SoAKernelKokkosStatic(x1, y1, z1, soa2Device, fx, fy, fz, cutoffSquared, i, j);
          },
          fxAcc, fyAcc, fzAcc);

      Kokkos::single(Kokkos::PerTeam(teamHandle), [&]() {
        soa1Device.template operator()<Particle_T::AttributeNames::forceX, true>(i) += fxAcc;
        soa1Device.template operator()<Particle_T::AttributeNames::forceY, true>(i) += fyAcc;
        soa1Device.template operator()<Particle_T::AttributeNames::forceZ, true>(i) += fzAcc;
      });
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

  const size_t _teamSize {0};
};

}