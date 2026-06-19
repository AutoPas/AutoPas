/**
 * @file KokkosDsChunksTraversal.h
 * @date 31. October 2025
 * @author Luis Gall
 */

#pragma once

#include "autopas/containers/kokkosDirectSum/traversals/DSKokkosTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

#include <Kokkos_Core.hpp>

namespace autopas {

template <class Functor, class Particle_T, class DeviceSpace>
struct ChunksTraversalFunctor {
  using FloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;
  using MemberType = Kokkos::TeamPolicy<typename DeviceSpace::execution_space>::member_type;

  utilsKokkos::KokkosStorage<Particle_T> _storageA;
  utilsKokkos::KokkosStorage<Particle_T> _storageB;
  Functor* _func;
  FloatPrecision _cutoffSquared;
  size_t _M;
  size_t _N;
  size_t _chunkSize;

  KOKKOS_INLINE_FUNCTION
  void operator()(const MemberType& teamHandle) const {
    const int k = teamHandle.league_rank();

    size_t offset = k * _chunkSize;
    size_t rest = _N - offset;
    size_t upper = rest < _chunkSize ? rest : _chunkSize;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamHandle, upper), [&](int i) {
        FloatPrecision fxAcc = 0.;
        FloatPrecision fyAcc = 0.;
        FloatPrecision fzAcc = 0.;

        const auto x1 = _storageA.template operator()<Particle_T::AttributeNames::posX, false>(i + offset);
        const auto y1 = _storageA.template operator()<Particle_T::AttributeNames::posY, false>(i + offset);
        const auto z1 = _storageA.template operator()<Particle_T::AttributeNames::posZ, false>(i + offset);

        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamHandle, _M), [&](int j,
            FloatPrecision& localFxAcc,
            FloatPrecision& localFyAcc,
            FloatPrecision& localFzAcc) {
                _func->ForceKernelKokkos(x1, y1, z1, _storageB, localFxAcc, localFyAcc, localFzAcc, _cutoffSquared, i, j);
            }, fxAcc, fyAcc, fzAcc);

        _storageA.template operator()<Particle_T::AttributeNames::forceX, false>(i + offset) += fxAcc;
        _storageA.template operator()<Particle_T::AttributeNames::forceY, false>(i + offset) += fyAcc;
        _storageA.template operator()<Particle_T::AttributeNames::forceZ, false>(i + offset) += fzAcc;
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
class KokkosDsChunksTraversal : public DSKokkosTraversalInterface<Particle_T> {

public:
  /**
   * Constructor for the KokkkosDSNaiveParallelTraversal
   * @param functor the functor that defines the interaction of particles
   * @param dataLayout The data layout wth which this traversal should be initialized
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not
   * @param teamSize Size of the Kokkos Teams that work on one chunk
   * @param chunkSize Size of the Chunks for the Kokkos Teams, one team works one chunkSize consecutive i-particles
   */
explicit KokkosDsChunksTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3, size_t teamSize, size_t chunkSize)
        : DSKokkosTraversalInterface<Particle_T>(dataLayout, useNewton3), _functor{functor}, _teamSize(teamSize), _chunkSize(chunkSize) {}

    [[nodiscard]] TraversalOption getTraversalType() const final { return TraversalOption::ds_kokkos_chunks; }

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

    const size_t chunkSize = _chunkSize;
    const size_t numChunks = N / chunkSize;

    ChunksTraversalFunctor<Functor, Particle_T, typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace> functor {storageA, storageB, func, cutoffSquared, M, N, chunkSize};

    auto teamPolicy = Kokkos::TeamPolicy<typename DSKokkosTraversalInterface<Particle_T>::DeviceSpace::execution_space>(numChunks+1, _teamSize, Kokkos::AUTO);
    Kokkos::parallel_for("autopas::ChunksTraversalSoA", teamPolicy, functor);
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

  const size_t _chunkSize {0};
};

}