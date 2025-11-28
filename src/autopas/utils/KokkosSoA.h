/**
 * @file KokkosSoA.h
 * @author Luis Gall
 * @date 13.11.2025
 */

#pragma once
#include <tuple>

#include <Kokkos_Core.hpp>

namespace autopas::utils {

  template <class MemSpace, typename ... Types>
  class KokkosSoA {

  public:

    explicit KokkosSoA()
      : KokkosSoA(0) {}

    explicit KokkosSoA(size_t numParticles)
    {
      resize(numParticles);
    }

    /* Get/Set/Allocation */
    void resize(size_t numParticles) {
      constexpr auto tupleSize = std::tuple_size<decltype(views)>::value;
      constexpr auto I = std::make_index_sequence<tupleSize>{};

      resizeImpl(numParticles, I);
    }

    template <class Particle_T>
    void addParticle(size_t position, const Particle_T& p) {
      constexpr auto tupleSize = std::tuple_size<decltype(views)>::value;
      constexpr auto I = std::make_index_sequence<tupleSize>{};

      addParticleImpl(position, p, I);
    }

    template <size_t attribute>
    constexpr std::tuple_element<attribute, std::tuple<Kokkos::View<Types, MemSpace>...>>::type::value_type get(size_t index) {
      return std::get<attribute>(views)(index);
    }

    template <size_t attribute>
        constexpr std::tuple_element<attribute, std::tuple<Kokkos::View<Types, MemSpace>...>>::type& getView() {
      return std::get<attribute>(views);
    }

    template <size_t attribute>
    void set(std::tuple_element<attribute, std::tuple<Kokkos::View<Types, MemSpace>...>>::type::value_type value, size_t index) {
      (std::get<attribute>(views))(index) = value;
    }

    /* Meta Data */
    size_t size() const {
      return std::get<0>(views).extent(0);
    }

    constexpr static size_t tupleSize() {
      return std::tuple_size<decltype(views)>::value;
    }

    /* Data copies */
    template <class SrcSoA, std::size_t ... I>
    void copyFrom(SrcSoA src, std::index_sequence<I...>) {
      (Kokkos::deep_copy(std::get<I>(views), src.template getView<I>()), ...);
    }

  private:
    template <std::size_t... I>
    void resizeImpl(size_t numParticles, std::index_sequence<I...>) {
      (Kokkos::realloc(std::get<I>(views), numParticles), ...);
    }

    template <class Particle_T, std::size_t... I>
    void addParticleImpl(size_t position, const Particle_T& p, std::index_sequence<I...>) {
      (set<I>(p.template get<static_cast<Particle_T::AttributeNames>(I+1)>(), position), ...);
    }

    std::tuple<Kokkos::View<Types, MemSpace>...> views {};
  };

}