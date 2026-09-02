/**
 * @file MoleculeLJ.cpp
 * @date Originally 17/01/2018. Code pulled from header on 15/06/2023.
 * @author tchipevn
 */

#include "MoleculeLJ.h"

namespace mdLib {
MoleculeLJ::MoleculeLJ(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long moleculeId,
                       unsigned long typeId)
    : autopas::ParticleBaseFP64(pos, v, moleculeId), _typeId(typeId) {}

const std::array<double, 3> &MoleculeLJ::getOldF() const { return _oldF; }
void MoleculeLJ::setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

size_t MoleculeLJ::getTypeId() const { return _typeId; }
void MoleculeLJ::setTypeId(size_t typeId) { _typeId = typeId; }

double MoleculeLJ::getSqrtEpsilon() const { return _sqrtEpsilon; }
void MoleculeLJ::setSqrtEpsilon(double sqrtEpsilon) { _sqrtEpsilon = sqrtEpsilon; }

double MoleculeLJ::getHalfSigma() const { return _halfSigma; }
void MoleculeLJ::setHalfSigma(double halfSigma) { _halfSigma = halfSigma; }

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
     << "\nSqrtEpsilon        : " << _sqrtEpsilon
     << "\nHalfSigma          : " << _halfSigma
     << "\nOwnershipState     : " << _ownershipState;
  // clang-format on
  return text.str();
}
}  // namespace mdLib
