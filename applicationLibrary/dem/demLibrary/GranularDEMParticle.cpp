/**
 * @file GranularDEM.cpp
 * @author Joon Kim
 * @date 27/03/2025
 */

#include "GranularDEMParticle.h"

namespace demLib {
GranularDEMParticle::GranularDEMParticle(const std::array<double, 3> &pos, const std::array<double, 3> &v,
                                         std::array<double, 3> angularVel, unsigned long particleId,
                                         unsigned long typeId, double temperature)
    : autopas::ParticleBaseFP64(pos, v, particleId),
      _angularVel(angularVel),
      _torque({0., 0., 0.}),
      _typeId(typeId),
      _temperature(temperature),
      _heatFlux(0.) {}

const std::array<double, 3> &GranularDEMParticle::getOldF() const { return _oldF; }
void GranularDEMParticle::setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

size_t GranularDEMParticle::getTypeId() const { return _typeId; }
void GranularDEMParticle::setTypeId(size_t typeId) { _typeId = typeId; }

const std::array<double, 3> &GranularDEMParticle::getAngularVel() const { return _angularVel; }
void GranularDEMParticle::setAngularVel(const std::array<double, 3> &angularVel) { _angularVel = angularVel; }
void GranularDEMParticle::addAngularVel(const std::array<double, 3> &angularVel) {
  _angularVel = autopas::utils::ArrayMath::add(_angularVel, angularVel);
}

const std::array<double, 3> &GranularDEMParticle::getTorque() const { return _torque; }
void GranularDEMParticle::setTorque(const std::array<double, 3> &torque) { _torque = torque; }
void GranularDEMParticle::addTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::add(_torque, torque);
}
void GranularDEMParticle::subTorque(const std::array<double, 3> &torque) {
  _torque = autopas::utils::ArrayMath::sub(_torque, torque);
}

double GranularDEMParticle::getTemperature() const { return _temperature; }
void GranularDEMParticle::setTemperature(double temperature) { _temperature = temperature; }
void GranularDEMParticle::addTemperature(double temperature) { _temperature += temperature; }

double GranularDEMParticle::getHeatFlux() const { return _heatFlux; }
void GranularDEMParticle::setHeatFlux(double heatFlux) { _heatFlux = heatFlux; }
void GranularDEMParticle::addHeatFlux(double heatFlux) { _heatFlux += heatFlux; }
void GranularDEMParticle::subHeatFlux(double heatFlux) { _heatFlux -= heatFlux; }

std::string GranularDEMParticle::toString() const {
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