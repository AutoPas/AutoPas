/**
 * @file KokkosAoS.h
 * @author Luis Gall
 * @date 27.11.2025
 */

#pragma once

#include <Kokkos_Core.hpp>
#include "Kokkos_DualView.hpp"

namespace autopas::utilsKokkos {

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
    void realloc(size_t numParticles) {
      Kokkos::realloc(view, numParticles);
    }

    void resize(size_t numParticles) {
      if (numParticles == 0) {
        return;
      } else if (view.extent(0) == 0) {
        realloc(numParticles);
      } else {
        Kokkos::resize(view, numParticles);
      }
    }

    template <typename Target>
    void sync() {
      view.template sync<Target>();
    }

    template <typename Target>
    void modify() {
      view.template modify<Target>();
    }

    template <bool useHostView>
    void addParticle(size_t index, const Particle_T& p) {
      if constexpr (useHostView) {
        view.view_host()(index) = p;
      } else {
        view.view_device()(index) = p;
      }
    }

    template <size_t attribute, bool useHostView>
    constexpr auto& operator() (size_t index) const {
      if constexpr(useHostView) {
        return view.view_host()(index).template operator()<static_cast<Particle_T::AttributeNames>(attribute)>();
      } else {
        return view.view_device()(index).template operator()<static_cast<Particle_T::AttributeNames>(attribute)>();
      }
    }

    auto& getView() {
      return view;
    }

    /* Meta Data */
    size_t size() const {
      return view.extent(0);
    }

    template <bool useHostView>
    Particle_T& getParticle (size_t index) {
      return useHostView? view.view_host()(index) : view.view_device()(index);
    }

    template <bool useHostView>
    const Particle_T& getParticle (size_t index) const {
      return useHostView? view.view_host()(index) : view.view_device()(index);
    }

  private:

#ifdef KOKKOS_ENABLE_CUDA
    using DeviceSpace = Kokkos::CudaSpace;
#else
    using DeviceSpace = Kokkos::HostSpace;
#endif

    // TODO: think about converting this to a Kokkos::DualView and allow to store AoS particles on the GPU
    Kokkos::DualView<Particle_T*, DeviceSpace::device_type> view {};
  };

}