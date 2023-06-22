/**
* @file MultisiteMoleculeLJ.cpp
* @date 14/02/2022. Code pulled from header on 15/06/2023
* @author S. Newcome
*/

#include "MultisiteMoleculeLJ.h"

namespace mdLib {
  MultisiteMoleculeLJ::MultisiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, std::array<double, 4> q,
                                         std::array<double, 3> angularVel, unsigned long moleculeId, unsigned long typeId)
    : mdLib::MoleculeLJ(r, v, moleculeId, typeId), _q(q), _angularVel(angularVel), _torque({0., 0., 0.}) {}

  template <MultisiteMoleculeLJ::AttributeNames attribute, std::enable_if_t<attribute != MultisiteMoleculeLJ::AttributeNames::ptr, bool>>
  constexpr typename std::tuple_element<attribute, MultisiteMoleculeLJ::SoAArraysType>::type::value_type MultisiteMoleculeLJ::get() const {
    if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::velocityX) {
      return getV()[0];
    } else if constexpr (attribute == AttributeNames::velocityY) {
      return getV()[1];
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      return getV()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      return getOldF()[0];
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      return getOldF()[1];
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      return getOldF()[2];
    } else if constexpr (attribute == AttributeNames::quaternion0) {
      return getQ()[0];
    } else if constexpr (attribute == AttributeNames::quaternion1) {
      return getQ()[1];
    } else if constexpr (attribute == AttributeNames::quaternion2) {
      return getQ()[2];
    } else if constexpr (attribute == AttributeNames::quaternion3) {
      return getQ()[3];
    } else if constexpr (attribute == AttributeNames::angularVelX) {
      return getAngularVel()[0];
    } else if constexpr (attribute == AttributeNames::angularVelY) {
      return getAngularVel()[1];
    } else if constexpr (attribute == AttributeNames::angularVelZ) {
      return getAngularVel()[2];
    } else if constexpr (attribute == AttributeNames::torqueX) {
      return getTorque()[0];
    } else if constexpr (attribute == AttributeNames::torqueY) {
      return getTorque()[1];
    } else if constexpr (attribute == AttributeNames::torqueZ) {
      return getTorque()[2];
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("MultisiteMoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  template <MultisiteMoleculeLJ::AttributeNames attribute>
  constexpr void MultisiteMoleculeLJ::set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::id) {
      setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::velocityX) {
      _v[0] = value;
    } else if constexpr (attribute == AttributeNames::velocityY) {
      _v[1] = value;
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      _v[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      _oldF[0] = value;
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      _oldF[1] = value;
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      _oldF[2] = value;
    } else if constexpr (attribute == AttributeNames::quaternion0) {
      _q[0] = value;
    } else if constexpr (attribute == AttributeNames::quaternion1) {
      _q[1] = value;
    } else if constexpr (attribute == AttributeNames::quaternion2) {
      _q[2] = value;
    } else if constexpr (attribute == AttributeNames::quaternion3) {
      _q[3] = value;
    } else if constexpr (attribute == AttributeNames::angularVelX) {
      _angularVel[0] = value;
    } else if constexpr (attribute == AttributeNames::angularVelY) {
      _angularVel[1] = value;
    } else if constexpr (attribute == AttributeNames::angularVelZ) {
      _angularVel[2] = value;
    } else if constexpr (attribute == AttributeNames::torqueX) {
      _torque[0] = value;
    } else if constexpr (attribute == AttributeNames::torqueY) {
      _torque[1] = value;
    } else if constexpr (attribute == AttributeNames::torqueZ) {
      _torque[2] = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      _typeId = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("MultisiteMoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  const std::array<double, 4> &MultisiteMoleculeLJ::getQ() const { return _q; }
  void MultisiteMoleculeLJ::setQ(const std::array<double, 4> &q) { _q = q; }

  const std::array<double, 3> &MultisiteMoleculeLJ::getAngularVel() const { return _angularVel; }
  void MultisiteMoleculeLJ::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
  void MultisiteMoleculeLJ::addAngularVel(const std::array<double, 3> &angularVel) {
    _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
  }

  const std::array<double, 3> &MultisiteMoleculeLJ::getTorque() const { return _torque; }
  void MultisiteMoleculeLJ::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
  void MultisiteMoleculeLJ::addTorque(const std::array<double, 3> &torque) { _torque = autopas::utils::ArrayMath::add(_torque, torque); }
  void MultisiteMoleculeLJ::subTorque(const std::array<double, 3> &torque) { _torque = autopas::utils::ArrayMath::sub(_torque, torque); }

  std::string MultisiteMoleculeLJ::toString() const  {
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
}