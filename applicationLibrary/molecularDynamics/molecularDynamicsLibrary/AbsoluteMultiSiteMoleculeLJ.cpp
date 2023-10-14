/**
* @file AbsoluteMultiSiteMoleculeLJ.h
* @date 10/10/2023
* @author Johannes Riemenschneider
 */

#include "AbsoluteMultiSiteMoleculeLJ.h"

namespace mdLib {
AbsoluteMultiSiteMoleculeLJ::AbsoluteMultiSiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, std::array<double, 4> q,
                                        std::array<double, 3> angularVel, unsigned long moleculeId,
                                        unsigned long typeId)
   : mdLib::MoleculeLJ(r, v, moleculeId, typeId), _q(q), _angularVel(angularVel), _torque({0., 0., 0.}) {}

const std::array<double, 4> &AbsoluteMultiSiteMoleculeLJ::getQuaternion() const { return _q; }
void AbsoluteMultiSiteMoleculeLJ::setQuaternion(const std::array<double, 4> &q) { _q = q; }

const std::array<double, 3> &AbsoluteMultiSiteMoleculeLJ::getAngularVel() const { return _angularVel; }
void AbsoluteMultiSiteMoleculeLJ::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
void AbsoluteMultiSiteMoleculeLJ::addAngularVel(const std::array<double, 3> &angularVel) {
 _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
}

const std::array<double, 3> &AbsoluteMultiSiteMoleculeLJ::getTorque() const { return _torque; }
void AbsoluteMultiSiteMoleculeLJ::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
void AbsoluteMultiSiteMoleculeLJ::addTorque(const std::array<double, 3> &torque) {
 _torque = autopas::utils::ArrayMath::add(_torque, torque);
}
void AbsoluteMultiSiteMoleculeLJ::subTorque(const std::array<double, 3> &torque) {
 _torque = autopas::utils::ArrayMath::sub(_torque, torque);
}

std::string AbsoluteMultiSiteMoleculeLJ::toString() const {
 using autopas::utils::ArrayUtils::operator<<;
 std::ostringstream text;
 // clang-format off
     text << "AbsoluteMultiSiteMoleculeLJ"
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