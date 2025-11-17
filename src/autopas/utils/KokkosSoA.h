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
      resizeSoA(numParticles);
    }

    void resizeSoA(size_t numParticles) {
      constexpr auto tupleSize = std::tuple_size<decltype(views)>::value;
      constexpr auto I = std::make_index_sequence<tupleSize>{} ;

      resizeSoAImpl(numParticles, I);
    }

    template <size_t attribute>
    constexpr typename std::tuple_element<attribute, std::tuple<Kokkos::View<Types, MemSpace>...>>::type::value_type get(size_t index) {
      return std::get<attribute>(views)(index);
    }

    template <size_t attribute>
    void set(typename std::tuple_element<attribute, std::tuple<Kokkos::View<Types, MemSpace>...>>::type::value_type value, size_t index) {
      (std::get<attribute>(views))(index) = value;
    }

  private:
    template <std::size_t... I>
    void resizeSoAImpl(size_t numParticles, std::index_sequence<I...>) {
      (Kokkos::realloc(std::get<I>(views), numParticles), ...);
    }

    std::tuple<Kokkos::View<Types, MemSpace>...> views {};
  };

}