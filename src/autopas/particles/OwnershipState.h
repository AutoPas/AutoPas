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
 */
enum class OwnershipState : char {
  /// Dummy or deleted state, a particle with this state is not an actual particle!
  dummy,
  /// Owned state, a particle with this state is an actual particle and owned by the current AutoPas object!
  owned,
  /// Halo state, a particle with this state is an actual particle, but not owned by the current AutoPas object!
  halo
};

}  // namespace autopas