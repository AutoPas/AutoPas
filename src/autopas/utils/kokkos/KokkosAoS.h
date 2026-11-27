/**
 * @file KokkosAoS.h
 * @author Luis Gall
 * @date 27.11.2025
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>
#include <type_traits>

#include "Kokkos_DualView.hpp"
#include "autopas/utils/kokkos/KokkosSpaces.h"

namespace autopas::utils::kokkos {

/**
 * Kokkos version of the Array of Structure Particle Vector
 *
 * @note This class does not manage the Kokkos::DualView sync state itself (unlike KokkosList). Callers have to uphold
 * this protocol:
 * 1. sync<Space>() before the first write in a space: after resize(), the device is marked modified, and
 *    modify<Space>() aborts if the other space is modified but not synced.
 * 2. modify<Space>() after every direct write (incl. addParticle()) and BEFORE any resize(): DualView::resize()
 *    discards the side that is not flagged as newest, so unflagged writes are silently lost.
 * Rule of thumb: sync → write → modify, never write → resize.
 *
 * @tparam Particle_T The particle type to store
 */
template <typename Particle_T>
class KokkosAoS {
  static_assert(std::is_trivially_copyable_v<Particle_T>,
                "KokkosAoS requires trivially copyable particles: elements are allocated WithoutInitializing and "
                "deep-copied between memory spaces, so e.g. virtual functions are not allowed.");

 public:
  /**
   * Operator to access the data of the requested particle attribute at the given position.
   * @tparam Attribute Attribute name
   * @tparam useDeviceView Flag to choose whether Kokkos should access the device-side memory
   * @param index particle index in the list
   * @return value of the particle attribute in the list at position index
   */
  template <typename Particle_T::AttributeNames Attribute, bool useDeviceView = true>
  KOKKOS_INLINE_FUNCTION auto &operator()(size_t index) const {
    if constexpr (useDeviceView) {
      return _particles.view_device()(index).template operator()<Attribute>();
    } else {
      return _particles.view_host()(index).template operator()<Attribute>();
    }
  }

  /**
   * Resizes the list to the given size
   * @param size new size of the list
   */
  void resize(size_t size) { _particles.resize(Kokkos::WithoutInitializing, size); }

  /**
   * Adds the given particle at the given index
   * @param index position where the new particle should be added
   * @param particle particle oject to add to the list
   */
  void addParticle(size_t index, const Particle_T &particle) { _particles.view_host()(index) = particle; }

  /**
   * Getter for the Kokkos::DualView<Particle>
   * @return the list of particles
   */
  const auto &getView() const { return _particles; }

  /**
   * Marks the given execution space as modified. Operation is delegated to Kokkos::DualView
   * @tparam Space Execution Space to mark as modified
   */
  template <typename Space>
  void modify() {
    _particles.template modify<Space>();
  }

  /**
   * Syncs the content to the given execution space. Operation is delegated to Kokkos::DualView
   * @tparam Space Execution Space to sync to
   */
  template <typename Space>
  void sync() {
    _particles.template sync<Space>();
  }

 private:
  /**
   * List of AoS particle data
   */
  Kokkos::DualView<Particle_T *, DeviceExecSpace::device_type> _particles{};
};

}  // namespace autopas::utils::kokkos

#endif