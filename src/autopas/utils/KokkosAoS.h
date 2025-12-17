/**
 * @file KokkosAoS.h
 * @author Luis Gall
 * @date 27.11.2025
 */

#pragma once

#include <Kokkos_Core.hpp>
#include "Kokkos_DualView.hpp"

namespace autopas::utils {

  template <class Particle_T>
  class KokkosAoS {

  public:

    // TODO: guarantee that MemSpace is not Cuda Space
    explicit KokkosAoS()
      : KokkosAoS(0) {}

    explicit KokkosAoS(size_t numParticles)
    {
      resize(numParticles);
    }

    /* Get/Set/Allocation */
    void resize(size_t numParticles) {
      if (numParticles == 0) {
        return;
      }
      Kokkos::realloc(view, numParticles);
    }

    void addParticle(size_t index, const Particle_T& p) {
      view(index) = p;
    }

    template <size_t attribute, bool>
    constexpr auto& operator() (int i) const {
      return view(i).template operator()<static_cast<Particle_T::AttributeNames>(attribute)>();
    }

    template <size_t attribute, bool>
    auto get(size_t index) {
      return view(index).template get<static_cast<Particle_T::AttributeNames>(attribute)>();
    }

    template <size_t attribute, bool>
    const auto get(size_t index) const {
      return view(index).template get<static_cast<Particle_T::AttributeNames>(attribute)>();
    }

    auto& getView() {
      return view;
    }

    template <size_t attribute, bool, typename Type>
    void set(Type value, size_t index) {
      view(index).template set<static_cast<Particle_T::AttributeNames>(attribute)>(value);
    }

    /* Meta Data */
    size_t size() const {
      return view.extent(0);
    }

    Particle_T& getParticle (size_t index) {
      return view(index);
    }

    const Particle_T& getParticle (size_t index) const {
      return view(index);
    }

  private:

    Kokkos::View<Particle_T*, Kokkos::HostSpace> view;
  };

}