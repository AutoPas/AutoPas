/**
* @file GranularDEM.cpp
* @author Joon Kim
* @date 27/03/2025
 */

#include "GranularDEM.h"

namespace demLib {
GranularDEM::GranularDEM(const std::array<double, 3> &pos, const std::array<double, 3> &v,
                         std::array<double, 3> angularVel, unsigned long particleId, unsigned long typeId, double temperature)
    : autopas::ParticleBaseFP64(pos, v, particleId), _angularVel(angularVel), _torque({0., 0., 0.}), _typeId(typeId), _temperature(temperature), _heatFlux(0.) {}

const std::array<double, 3> &GranularDEM::getOldF() const { return _oldF; }
void GranularDEM::setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

size_t GranularDEM::getTypeId() const { return _typeId; }
void GranularDEM::setTypeId(size_t typeId) { _typeId = typeId; }

const std::array<double, 3> &GranularDEM::getAngularVel() const { return _angularVel; }
void GranularDEM::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
void GranularDEM::addAngularVel(const std::array<double, 3> &angularVel) {
  _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
}

const std::array<double, 3> &GranularDEM::getTorque() const { return _torque; }
void GranularDEM::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
void GranularDEM::addTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::add(_torque, torque);
}
void GranularDEM::subTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::sub(_torque, torque);
}

double GranularDEM::getTemperature() const { return _temperature; }
void GranularDEM::setTemperature(double temperature) { _temperature = temperature; }
void GranularDEM::addTemperature(double temperature) { _temperature += temperature; }

double GranularDEM::getHeatFlux() const { return _heatFlux; }
void GranularDEM::setHeatFlux(double heatFlux) { _heatFlux = heatFlux; }
void GranularDEM::addHeatFlux(double heatFlux) { _heatFlux += heatFlux; }
void GranularDEM::subHeatFlux(double heatFlux) { _heatFlux -= heatFlux; }

std::string GranularDEM::toString() const {
  using autopas::utils::ArrayUtils::operator<<;
  std::ostringstream text;
  // clang-format off
      text << "GranularDEM"
         << "\nID                 : " << _id
         << "\nPosition           : " << _r
         << "\nVelocity           : " << _v
         << "\nForce              : " << _f
         << "\nOld Force          : " << _oldF
         << "\nAngular Velocity   : " << _angularVel
         << "\nTorque             : " << _torque
         << "\nType ID            : " << _typeId
         << "\nOwnershipState     : " << _ownershipState
         << "\nTemperature        : " << _temperature
         << "\nHeat Flux          : " << _heatFlux;
  // clang-format on
  return text.str();
}
}  // namespace demLib