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
#if AUTOPAS_PRECISION_MODE == SPSP
using OwnershipType = int32_t;
#elif AUTOPAS_PRECISION_MODE == DPDP
using OwnershipType = int64_t;
#else
using OwnershipType = int32_t;
#endif
/**
 * Enum that specifies the state of ownership.
 * @note this type has either int32_t or int64_t depending on the precision selected as an underlying type to be
 * compatible with LJFunctorAVX. For that a size equal to the precision of the particles is required.
 */
enum class OwnershipState : OwnershipType {
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
  return static_cast<OwnershipState>(static_cast<OwnershipType>(a) & static_cast<OwnershipType>(b));
}

/**
 * Bitwise OR operator for OwnershipState
 * @param a first operand
 * @param b second operand
 * @return a | b
 */
constexpr OwnershipState operator|(const OwnershipState a, const OwnershipState b) {
  return static_cast<OwnershipState>(static_cast<OwnershipType>(a) | static_cast<OwnershipType>(b));
}

/**
 * Returns the OwnershipType value of a given OwnershipState
 *
 * @param a OwnershipState
 * @return const OwnershipType value of a given OwnershipState
 */
constexpr OwnershipType toInt64(const OwnershipState a) { return static_cast<OwnershipType>(a); }

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