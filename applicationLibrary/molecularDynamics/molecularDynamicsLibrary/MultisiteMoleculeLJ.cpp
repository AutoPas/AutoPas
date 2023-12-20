/**
 * @file MultisiteMoleculeLJ.cpp
 * @date 14/02/2022. Code pulled from header on 15/06/2023
 * @author S. Newcome
 */

#include "MultisiteMoleculeLJ.h"

namespace mdLib {
MultisiteMoleculeLJ::MultisiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, std::array<double, 4> q,
                                         std::array<double, 3> angularVel, unsigned long moleculeId,
                                         unsigned long typeId)
    : mdLib::MoleculeLJ(r, v, moleculeId, typeId), _q(q), _angularVel(angularVel), _torque({0., 0., 0.}) {}

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

#if defined(MD_FLEXIBLE_TORQUE_AFTER_FORCE) and MD_FLEXIBLE_MODE==MULTISITE
void MultisiteMoleculeLJ::setForcesOnSitesX(std::vector<double> &&forcesOnSitesX);
_forcesOnSitesX = forcesOnSitesX;
}

void MultisiteMoleculeLJ::setForcesOnSitesY(std::vector<double> &&forcesOnSitesY);
_forcesOnSitesY = forcesOnSitesY;
}

void MultisiteMoleculeLJ::setForcesOnSitesZ(std::vector<double> &&forcesOnSitesZ);
_forcesOnSitesZ = forcesOnSitesZ;
}

const std::array<double, 3> MultisiteMoleculeLJ::getForceOnSite(size_t siteIndex) {
  return {_forcesOnSitesX[siteIndex], _forcesOnSitesY[siteIndex], _forcesOnSitesZ[siteIndex]};
}

void MultisiteMoleculeLJ::addToForceOnSite(size_t siteIndex, std::array<double, 3> &force) {
  _forcesOnSitesX[siteIndex] += force[0];
  _forcesOnSitesY[siteIndex] += force[1];
  _forcesOnSitesZ[siteIndex] += force[2];
}

void MultisiteMoleculeLJ::subToForceOnSite(size_t siteIndex, std::array<double, 3> &force) {
  _forcesOnSitesX[siteIndex] -= force[0];
  _forcesOnSitesY[siteIndex] -= force[1];
  _forcesOnSitesZ[siteIndex] -= force[2];
}

#endif
}  // namespace mdLib