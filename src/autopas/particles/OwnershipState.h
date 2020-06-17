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
 * @note This type has unsigned char as an underlying type to be compatible with LJFunctorAVX. You may only change this
 * type, when you adapt LJFunctorAVX accordingly.
 */
enum class OwnershipState : unsigned char {
  /// Dummy or deleted state, a particle with this state is not an actual particle!
  dummy = 0,
  /// Owned state, a particle with this state is an actual particle and owned by the current AutoPas object!
  owned = 1,
  /// Halo state, a particle with this state is an actual particle, but not owned by the current AutoPas object!
  halo = 2
};

}  // namespace autopas