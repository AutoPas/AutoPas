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

const std::array<double, 4> &MultisiteMoleculeLJ::getQ() const { return _q; }
void MultisiteMoleculeLJ::setQ(const std::array<double, 4> &q) { _q = q; }

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
}  // namespace mdLib