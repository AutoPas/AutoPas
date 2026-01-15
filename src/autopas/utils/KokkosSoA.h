/**
 * @file KokkosSoA.h
 * @author Luis Gall
 * @date 13.11.2025
 */

#pragma once
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <tuple>

namespace autopas::utils {

  template <typename ... Types>
  class KokkosSoA {

  public:

    explicit KokkosSoA()
      : KokkosSoA(0) {}

    explicit KokkosSoA(size_t N) {
      resize(N);
    }

    KokkosSoA(const KokkosSoA<Types...>& other) {
      views = other.views;
      std::cout << "Explicit copy of KokkosSoA" << std::endl;
    }

    KOKKOS_FUNCTION
    KokkosSoA(size_t N, const std::string& label) : views{Kokkos::DualView<Types>(label, N)...} {}

    /* Get/Set/Allocation */
    void resize(size_t numParticles) {
      if (numParticles == 0) {
        return;
      }
      constexpr auto tupleSize = std::tuple_size<decltype(views)>::value;
      constexpr auto I = std::make_index_sequence<tupleSize>{};

      resizeImpl(numParticles, I);
    }

    template <size_t attribute, bool offset, bool host = false>
    KOKKOS_INLINE_FUNCTION
    constexpr auto& operator() (int i) const {
      if constexpr (host) {
        return std::get<attribute - (offset ? 1 : 0)>(views).view_host()(i);
      }
      else {
        return std::get<attribute - (offset ? 1 : 0)>(views).view_device()(i);
      }
    }

    template <class Particle_T>
    void addParticle(size_t position, const Particle_T& p) {
      constexpr auto tupleSize = std::tuple_size<decltype(views)>::value;
      constexpr auto I = std::make_index_sequence<tupleSize>{};

      addParticleImpl(position, p, I);
    }

    template <size_t attribute>
    KOKKOS_INLINE_FUNCTION
    constexpr auto& getView() const {
      return std::get<attribute>(views);
    }

    /* Meta Data */
    KOKKOS_INLINE_FUNCTION
    size_t size() const {
      return std::get<0>(views).extent(0);
    }

    constexpr static size_t tupleSize() {
      return std::tuple_size<decltype(views)>::value;
    }

    /* Data copies */
    /*
    template <class SrcSoA, std::size_t ... I>
    void copyFrom(SrcSoA src, std::index_sequence<I...>) {
      (Kokkos::deep_copy(std::get<I>(views), src.template getView<I>()), ...);
    }
    */

    template <typename Target, std::size_t... I>
    void markModified(std::index_sequence<I...>) {
      (std::get<I>(views).template modify<Target>(), ...);
    }

    template <typename Target, std::size_t... I>
    void sync(std::index_sequence<I...>) {
      (std::get<I>(views).template sync<Target>(), ...);
    }

    void operator= (KokkosSoA<Types...> &other) {
      views = other.views;
    }

    void operator= (const KokkosSoA<Types...> &other) {
      views = other.views;
    }

  private:
    template <std::size_t... I>
    void resizeImpl(size_t numParticles, std::index_sequence<I...>) {
      (std::get<I>(views).resize(numParticles), ...);
    }

    template <class Particle_T, std::size_t... I>
    void addParticleImpl(size_t position, const Particle_T& p, std::index_sequence<I...>) {
      ((operator()<I, false, true>(position) = p.template get<static_cast<Particle_T::AttributeNames>(I+1)>()), ...);
    }

#ifdef KOKKOS_ENABLE_CUDA
    std::tuple<Kokkos::DualView<Types, Kokkos::CudaSpace::device_type>...> views {};
#else
    std::tuple<Kokkos::DualView<Types, Kokkos::HostSpace::device_type>...> views {};
#endif
  };

}