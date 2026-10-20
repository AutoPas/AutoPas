/**
 * @file MultisiteMoleculeLJ.cpp
 * @date 14/02/2022. Code pulled from header on 15/06/2023
 * @author S. Newcome
 */

#include "MultisiteMoleculeLJ.h"
#include "autopas/utils/Quaternion.h"

namespace mdLib {
MultisiteMoleculeLJ::MultisiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, std::array<double, 4> q,
                                         std::array<double, 3> angularVel, unsigned long moleculeId,
                                         unsigned long typeId, const static ParticlePropertiesLibrary<> &PPL)
    : MoleculeLJ(r, v, moleculeId, typeId), _q(q), _angularVel(angularVel), _torque({0., 0., 0.}) {
  _PPL = PPL;
}

const std::array<double, 4> &MultisiteMoleculeLJ::getQuaternion() const { return _q; }
void MultisiteMoleculeLJ::setQuaternion(const std::array<double, 4> &q) { _q = q; }

const std::array<double, 3> &MultisiteMoleculeLJ::getAngularVel() const { return _angularVel; }
void MultisiteMoleculeLJ::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
void MultisiteMoleculeLJ::addAngularVel(const std::array<double, 3> &angularVel) {
  _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
}

const std::array<double, 3> &MultisiteMoleculeLJ::getTorque() const { return _torque; }
void MultisiteMoleculeLJ::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
void MultisiteMoleculeLJ::addTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::add(_torque, torque);
}
void MultisiteMoleculeLJ::subTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::sub(_torque, torque);
}

const double MultisiteMoleculeLJ::getSitePosX(size_t i) const {
  const auto unrotatedUntranslatedSitePos = this->PPL.getSitePositions(_typeId)[i];
  const auto rotatedUntranslatedSitePosX = autopas::utils::quaternion::rotatePositionX(_q, unrotatedUntranslatedSitePos);
  return rotatedUntranslatedSitePosX + _r[0];
}

const double MultisiteMoleculeLJ::getSitePosY(size_t i) const {
  const auto unrotatedUntranslatedSitePos = PPL.getSitePositions(_typeId)[i];
  const auto rotatedUntranslatedSitePosY = autopas::utils::quaternion::rotatePositionY(_q, unrotatedUntranslatedSitePos);
  return rotatedUntranslatedSitePosY + _r[1];
}

const double MultisiteMoleculeLJ::getSitePosZ(size_t i) const {
  const auto unrotatedUntranslatedSitePos = PPL.getSitePositions(_typeId)[i];
  const auto rotatedUntranslatedSitePosZ = autopas::utils::quaternion::rotatePositionZ(_q, unrotatedUntranslatedSitePos);
  return rotatedUntranslatedSitePosZ + _r[2];
}

const size_t &MultisiteMoleculeLJ::getSiteTypeId(size_t i) const {
  return PPL.getSiteTypes(_typeId)[i];
}

void MultisiteMoleculeLJ::addSiteForceX(const size_t i, const double &siteForceX) {
  _f[0] += siteForceX;
  _torque[0] += autopas::utils::ArrayMath::
}

std::string MultisiteMoleculeLJ::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  // clang-format off
      text << "MultisiteMoleculeLJ"
         << "\nID                 : " << _id
         << "\nPosition           : " << _r
         << "\nVelocity           : " << _v
         << "\nForce              : " << _f
         << "\nOld Force          : " << _oldF
         << "\nQuaternion         : " << _q
         << "\nRotational Velocity: " << _angularVel
         << "\nTorque             : " << _torque
         << "\nType ID            : " << _typeId
         << "\nOwnershipState     : " << _ownershipState;
  // clang-format on
  return text.str();
}
}  // namespace mdLib