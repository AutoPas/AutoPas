/**
 * @file OwnershipState.h
 *
 * @author seckler
 * @date 26.05.20
 */

#pragma once

namespace autopas {
/**
 * Enum that specifies the state of ownership.
 * @note this type has int64_t as an underlying type to be compatible with LJFunctorAVX. For that a size equal to the
 * precision of the particles is required.
 */
enum class OwnershipState : int64_t {
  /// Dummy or deleted state, a particle with this state is not an actual particle!
  /// @note LJFunctorAVX requires that the Dummy state should always be the integer zero and the state with the lowest
  /// value.
  dummy = 0,
  /// Owned state, a particle with this state is an actual particle and owned by the current AutoPas object!
  owned,
  /// Halo state, a particle with this state is an actual particle, but not owned by the current AutoPas object!
  halo
};

}  // namespace autopas