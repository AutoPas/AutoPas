/**
 * @file MoleculeLJ.cpp
 * @date Originally 17/01/2018. Code pulled from header on 15/06/2023.
 * @author tchipevn
 */

#include "MoleculeLJ.h"

namespace mdLib {
MoleculeLJ::MoleculeLJ(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long moleculeId,
                       unsigned long typeId)
    : MoleculeLJPrecisionBase(), _typeId(typeId) {
  std::array<MoleculeLJ::AoSCalcType, 3> posCast{};
  std::array<MoleculeLJ::AoSCalcType, 3> vCast{};

  for (std::size_t d = 0; d < 3; ++d) {
    posCast[d] = static_cast<MoleculeLJ::AoSCalcType>(pos[d]);
    vCast[d] = static_cast<MoleculeLJ::AoSCalcType>(v[d]);
  }

  this->setR(posCast);
  this->setV(vCast);
  this->setID(moleculeId);
}

const std::array<MoleculeLJ::ParticleSoAFloatPrecision, 3> &MoleculeLJ::getOldF() const { return _oldF; }
void MoleculeLJ::setOldF(const std::array<MoleculeLJ::ParticleSoAFloatPrecision, 3> &oldForce) { _oldF = oldForce; }

size_t MoleculeLJ::getTypeId() const { return _typeId; }
void MoleculeLJ::setTypeId(size_t typeId) { _typeId = typeId; }

std::string MoleculeLJ::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  // clang-format off
  text << "MoleculeLJ"
     << "\nID                 : " << _id
     << "\nPosition           : " << _r
     << "\nVelocity           : " << _v
     << "\nForce              : " << _f
     << "\nOld Force          : " << _oldF
     << "\nType ID            : " << _typeId
     << "\nOwnershipState     : " << _ownershipState;
  // clang-format on
  return text.str();
}
}  // namespace mdLib
