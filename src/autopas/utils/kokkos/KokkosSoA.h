/**
 * @file KokkosSoA.h
 * @author Luis Gall
 * @date 13.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>
#include <type_traits>

#include "Kokkos_DualView.hpp"
#include "autopas/utils/kokkos/KokkosSpaces.h"

namespace autopas::utils::kokkos {

/**
 * Structure of arrays with one DualView per attribute.
 *
 * All public methods address columns by attribute name. The columns store the attributes [firstStored, firstStored +
 * sizeof...(AttributeTypes)) in the order of AttributeNames_T, so attributes before firstStored (e.g. ptr) are not
 * stored. The conversion between attribute names and column indices happens only in this class.
 *
 * @tparam AttributeNames_T Enum that names the attributes, usually Particle_T::AttributeNames.
 * @tparam firstStored Attribute that is stored in the first column.
 * @tparam AttributeTypes Data types of the columns, e.g. float *.
 */
template <typename AttributeNames_T, AttributeNames_T firstStored, typename... AttributeTypes>
class KokkosSoA {
  static_assert(std::is_enum_v<AttributeNames_T>, "KokkosSoA columns have to be named by an enum.");

 public:
  using Type = std::tuple<Kokkos::DualView<AttributeTypes, DeviceExecSpace::device_type>...>;

  template <AttributeNames_T Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto &operator()(size_t index) const {
    if constexpr (useDeviceView) {
      return std::get<column<Attribute>()>(_particleData).view_device()(index);
    } else {
      return std::get<column<Attribute>()>(_particleData).view_host()(index);
    }
  }

  void resize(size_t size) { resizeImpl(size, Columns{}); }

  template <class Particle_T>
  void addParticle(size_t index, const Particle_T &particle) {
    static_assert(std::is_same_v<typename Particle_T::AttributeNames, AttributeNames_T>,
                  "The particle's attribute names have to match the ones of the SoA.");
    addParticleImpl<Particle_T>(index, particle, Columns{});
  }

  template <AttributeNames_T Attribute>
  const auto &getView() const {
    return std::get<column<Attribute>()>(_particleData);
  }

  template <typename Space, AttributeNames_T Attribute>
  void modify() {
    modifyColumn<Space, column<Attribute>()>();
  }

  template <typename Space>
  void modifyAll() {
    modifyAllImpl<Space>(Columns{});
  }

  template <typename Space, AttributeNames_T Attribute>
  void sync() {
    syncColumn<Space, column<Attribute>()>();
  }

  template <typename Space>
  void syncAll() {
    syncAllImpl<Space>(Columns{});
  }

 private:
  using Columns = std::index_sequence_for<AttributeTypes...>;

  /**
   * The only place that converts an attribute name to a column index.
   * @tparam Attribute
   * @return Index of the column that stores Attribute.
   */
  template <AttributeNames_T Attribute>
  static constexpr size_t column() {
    static_assert(Attribute >= firstStored, "Attribute is not stored in this SoA (e.g. ptr).");
    constexpr size_t col = static_cast<size_t>(Attribute) - static_cast<size_t>(firstStored);
    static_assert(col < sizeof...(AttributeTypes), "Attribute is out of range for this SoA.");
    return col;
  }

  /**
   * The only place that converts a column index to an attribute name.
   * @tparam Col
   * @return Attribute stored in column Col.
   */
  template <size_t Col>
  static constexpr AttributeNames_T attributeOf() {
    return static_cast<AttributeNames_T>(Col + static_cast<size_t>(firstStored));
  }

  template <typename Space, size_t Col>
  void modifyColumn() {
    std::get<Col>(_particleData).template modify<Space>();
  }

  template <typename Space, size_t Col>
  void syncColumn() {
    std::get<Col>(_particleData).template sync<Space>();
  }

  template <size_t... Col>
  void resizeImpl(size_t size, std::index_sequence<Col...>) {
    (std::get<Col>(_particleData).resize(Kokkos::WithoutInitializing, size), ...);
  }

  template <typename Space, size_t... Col>
  void syncAllImpl(std::index_sequence<Col...>) {
    (syncColumn<Space, Col>(), ...);
  }

  template <typename Space, size_t... Col>
  void modifyAllImpl(std::index_sequence<Col...>) {
    (modifyColumn<Space, Col>(), ...);
  }

  template <class Particle_T, size_t... Col>
  void addParticleImpl(size_t index, const Particle_T &particle, std::index_sequence<Col...>) {
    ((std::get<Col>(_particleData).view_host()(index) = particle.template get<attributeOf<Col>()>()), ...);
  }

  Type _particleData{};
};

}  // namespace autopas::utils::kokkos

#endif
