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
 * This function is needed to pass an ownership state to spdlog/fmt.
 * @param state the OwnershipState to format
 * @return string representation of the ownership state
 */
inline std::string format_as(const OwnershipState &state) {
  switch (state) {
    case OwnershipState::dummy:
      return "dummy";
    case OwnershipState::owned:
      return "owned";
    case OwnershipState::halo:
      return "halo";
    default:
      return "unknown state: 0b" + std::bitset<4>(static_cast<int64_t>(state)).to_string();
  }
}

/**
 * Insertion operator for OwnershipState.
 * This function enables passing ownershipState to an ostream via `<<`.
 * @param os
 * @param ownershipState
 * @return os
 */
inline std::ostream &operator<<(std::ostream &os, const OwnershipState &ownershipState) {
  // Reuse format_as to prevent duplicated logic
  return os << format_as(ownershipState);
}
}  // namespace autopas