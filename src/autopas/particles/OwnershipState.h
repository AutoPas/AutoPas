/**
 * @file OwnershipState.h
 *
 * @author seckler
 * @date 26.05.20
 */

#pragma once

#include <bitset>
#include <iostream>

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
  /* TODO: check whether it might be possible to relax the prev constraint and have it as the highest value 0b100 as this would align to the IteratorBehaviors */
  /* TODO: Yes, LJFunctorAVX wants dummy to be 0 BUT it does not need to be the lowest value (at least I did not see why) */
  /* TODO: maybe change first value to an unused entry such that dummy can be 0b0100 */
  dummy = 0b0000,  // 0
  /// Owned state, a particle with this state is an actual particle and owned by the current AutoPas object!
  owned = 0b0001,  // 1
  /// Halo state, a particle with this state is an actual particle, but not owned by the current AutoPas object!
  halo = 0b0010,  // 2
};

/**
 * Bitwise AND operator for OwnershipState
 * @param a first operand
 * @param b second operand
 * @return a & b
 */
constexpr OwnershipState operator&(const OwnershipState a, const OwnershipState b) {
  return static_cast<OwnershipState>(static_cast<int64_t>(a) & static_cast<int64_t>(b));
}

/**
 * Bitwise OR operator for OwnershipState
 * @param a first operand
 * @param b second operand
 * @return a | b
 */
constexpr OwnershipState operator|(const OwnershipState a, const OwnershipState b) {
  return static_cast<OwnershipState>(static_cast<int64_t>(a) | static_cast<int64_t>(b));
}

/**
 * Returns the int64_t value of a given OwnershipState
 *
 * @param a OwnershipState
 * @return const int64_t value of a given OwnershipState
 */
constexpr int64_t toInt64(const OwnershipState a) { return static_cast<int64_t>(a); }

/**
 * Insertion operator for OwnershipState.
 * This function enables passing ownershipState to an ostream via `<<`.
 * @param os
 * @param ownershipState
 * @return os
 */
[[maybe_unused]] static std::ostream &operator<<(std::ostream &os, const OwnershipState &ownershipState) {
  switch (ownershipState) {
    case OwnershipState::dummy:
      os << "dummy";
      break;
    case OwnershipState::owned:
      os << "owned";
      break;
    case OwnershipState::halo:
      os << "halo";
      break;
    default:
      os << "unknown state: 0b" << std::bitset<4>(static_cast<int64_t>(ownershipState));
      break;
  }
  return os;
}
}  // namespace autopas