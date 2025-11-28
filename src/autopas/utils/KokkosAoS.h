/**
 * @file KokkosAoS.h
 * @author Luis Gall
 * @date 27.11.2025
 */

#pragma once

#include <Kokkos_Core.hpp>

namespace autopas::utils {

  template <class MemSpace, class Particle_T>
  class KokkosAoS {

  public:

    explicit KokkosAoS()
      : KokkosAoS(0) {}

    explicit KokkosAoS(size_t numParticles)
    {
      resize(numParticles);
    }

    /* Get/Set/Allocation */
    void resize(size_t numParticles) {
      Kokkos::realloc(view, numParticles);
    }

    template <size_t attribute>
    auto get(size_t index) {
      return view(index).template get<static_cast<Particle_T::AttributeNames>(attribute)>();
    }

    auto& getView() {
      return view;
    }

    template <size_t attribute, typename Type>
    void set(Type value, size_t index) {
      view(index).template set<static_cast<Particle_T::AttributeNames>(attribute)>(value);
    }

    /* Meta Data */
    size_t size() const {
      return view.extent(0);
    }

    Particle_T& operator() (size_t index) {
      return view(index);
    }

    const Particle_T& operator() (size_t index) const {
      return view(index);
    }

  private:

    Kokkos::View<Particle_T*, MemSpace> view;
  };

}