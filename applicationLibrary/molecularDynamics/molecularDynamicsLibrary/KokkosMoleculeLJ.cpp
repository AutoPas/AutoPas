/**
 * @file KokkosMoleculeLJ.cpp
 * @date 23 Jan 2026 copied from MoleculeLJ.h
 * @author Luis Gall
 */

#include "KokkosMoleculeLJ.h"

namespace mdLib {

size_t KokkosMoleculeLJ::getTypeId() const { return _typeId; }
void KokkosMoleculeLJ::setTypeId(size_t typeId) { _typeId = typeId; }

std::string KokkosMoleculeLJ::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  // clang-format off
  text << "MoleculeLJ"
     << "\nID                 : " << _id
    /*
     << "\nPosition           : " << _r
     << "\nVelocity           : " << _v
     << "\nForce              : " << _f
     << "\nOld Force          : " << _oldF
     << "\nType ID            : " << _typeId
     */
     << "\nOwnershipState     : " << _ownershipState;
  // clang-format on
  return text.str();
}
}  // namespace mdLib
