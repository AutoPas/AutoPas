/**
 * @file SPHParticleBase.h
 * @author seckler
 * @date 19.01.18
 */

#pragma once

#include <cstring>
#include <vector>

#include "autopas/particles/Particle.h"

namespace sphLib {
/**
 * Basic SPHParticleBase class.
 */
template <typename CalcType, typename AccuType, typename idType>
class SPHParticleBase : public autopas::ParticleBase<CalcType, AccuType, idType> {
 public:
  /**
   * Default constructor of SPHParticleBase.
   * Will initialize all values to some basic defaults.
   */
  SPHParticleBase()
      : autopas::ParticleBase<CalcType, AccuType, idType>(),
        _density(0.),
        _pressure(0.),
        _mass(0.),
        _smth(0.),
        _snds(0.),
        // temporaries / helpers
        _v_sig_max(0.),
        _acc{0., 0., 0.},
        _energy_dot(0.),
        _energy(0.),
        _dt(0.),
        _vel_half{0., 0., 0.},
        _eng_half(0.) {}
  /**
   * Constructor of the SPHParticleBase class.
   * @param r position of the particle
   * @param v velocity of the particle
   * @param id id of the particle. This id should be unique
   */
  SPHParticleBase(const std::array<CalcType, 3> &r, const std::array<CalcType, 3> &v, idType id)
      : autopas::ParticleBase<CalcType, AccuType, idType>(r, v, id),
        _density(0.),
        _pressure(0.),
        _mass(0.),
        _smth(0.),
        _snds(0.),
        // temporaries / helpers
        _v_sig_max(0.),
        _acc{0., 0., 0.},
        _energy_dot(0.),
        _energy(0.),
        _dt(0.),
        _vel_half{0., 0., 0.},
        _eng_half(0.) {}

  /**
   * Constructor of the SPHParticleBase class.
   * @param r position of the particle
   * @param v velocity of the particle
   * @param id id of the particle. This id should be unique
   * @param mass mass of the particle
   * @param smth smoothing length of the particle
   * @param snds speed of sound (SouND Speed)
   */
  SPHParticleBase(const std::array<CalcType, 3> &r, const std::array<CalcType, 3> &v, unsigned long id, double mass,
                  double smth, double snds)
      : autopas::ParticleBase<CalcType, AccuType, idType>(r, v, id),
        _density(0.),
        _pressure(0.),
        _mass(mass),
        _smth(smth),
        _snds(snds),
        // temporaries / helpers
        _v_sig_max(0.),
        _acc{0., 0., 0.},
        _energy_dot(0.),
        _energy(0.),
        _dt(0.),
        _vel_half{0., 0., 0.},
        _eng_half(0.) {}

  /**
   * Destructor of the SPHParticleBase
   */
  ~SPHParticleBase() override = default;

  /**
   * Getter for the Density
   * @return the current density of the particle
   */
  double getDensity() const { return _density; }

  /**
   * Adds the given density to the current density
   * @param density density to be added
   */
  void addDensity(double density) { _density += density; }

  /**
   * Setter for Density
   * @param density The value of the density to be set as the particle's density
   */
  void setDensity(double density) { _density = density; }

  /**
   * Getter for Pressure
   * @return current pressure of the particle
   */
  double getPressure() const { return _pressure; }

  /**
   * Calculates the pressure within the particle from the energy and density of
   * the particle and updates the pressure and sound of speed
   */
  void calcPressure() {
    const double hcr = 1.4;
    _pressure = (hcr - 1.0) * _density * _energy;
    _snds = sqrt(hcr * _pressure / _density);
  };

  /**
   * Setter for the pressure
   * @param pressure pressure value to be set
   */
  void setPressure(double pressure) { _pressure = pressure; }

  /**
   * Getter for the mass of the particle
   * @return mass of particle
   */
  CalcType getMass() const { return _mass; }

  /**
   * Setter for the mass of the particle
   * @param mass mass to be set
   */
  void setMass(CalcType mass) { _mass = mass; }

  /**
   * Getter for the smoothing length of the particle
   * @return the smoothing length of the particle
   */
  double getSmoothingLength() const { return _smth; }

  /**
   * Setter for the smoothing length
   * @param smth smoothing lenth to be set
   */
  void setSmoothingLength(double smth) { _smth = smth; }

  /**
   * Getter for the speed of sound of the particle
   * @return speed of sound of the particle
   */
  double getSoundSpeed() const { return _snds; }

  /**
   * Setter for the speed of sound of the particle
   * @param snds speed of sound of the particle
   */
  void setSoundSpeed(double snds) { _snds = snds; }

  /**
   * Getter for the current maximally allowed signal velocity of the particle
   * @return the maximally allowed signal velocity of the particle
   */
  double getVSigMax() const { return _v_sig_max; }

  /**
   * Checks if the given signal velocity is higher than the current (local) one
   * and updates the local one if it is.
   * @param v_sig given signal velocity
   */
  void checkAndSetVSigMax(double v_sig) { _v_sig_max = std::max(v_sig, _v_sig_max); }

  /**
   * Setter for the maximally allowed signal velocity
   * @param v_sig_max the maximally allowed signal velocity
   */
  void setVSigMax(double v_sig_max) { _v_sig_max = v_sig_max; }

  /**
   * Getter for the acceleration of the particle
   * @return the acceleration of the particle
   */
  const std::array<double, 3> &getAcceleration() const { return _acc; }

  /**
   * Adds the given acceleration on the local acceleration.
   * Used to sum up different acceleration values.
   * @param acc Acceleration to be added
   */
  void addAcceleration(const std::array<double, 3> &acc) {
    using namespace autopas::utils::ArrayMath::literals;
    _acc += acc;
  };

  /**
   * Substracts the given acceleration from the local acceleration.
   * Used to sum up different negative acceleration values.
   * @param acc Acceleration to be substracted
   */
  void subAcceleration(const std::array<double, 3> &acc) {
    using namespace autopas::utils::ArrayMath::literals;
    _acc -= acc;
  };

  /**
   * Setter for the acceleration
   * @param acc Acceleration to be set
   */
  void setAcceleration(const std::array<double, 3> &acc) { _acc = acc; }

  /**
   * Getter for the time derivative of the energy of the particle
   * @return the time derivative of the energy of the particle
   */
  double getEngDot() const { return _energy_dot; }

  /**
   * Adds the given value to the current value of the time derivative of the
   * energy
   * @param eng_dot
   */
  void addEngDot(double eng_dot) { _energy_dot += eng_dot; }

  /**
   * Setter for the time derivative of the energy
   * @param eng_dot
   */
  void setEngDot(double eng_dot) { _energy_dot = eng_dot; }

  /**
   * Getter for the energy of the particle
   * @return the energy of the particle
   */
  double getEnergy() const { return _energy; }

  /**
   * Setter for the energy of the particle
   * @param energy the energy of the particle
   */
  void setEnergy(double energy) { _energy = energy; }

  /**
   * Adds the given energy to the energy of the particle
   * @param energy the energy to be added
   */
  void addEnergy(double energy) { _energy += energy; }

  /**
   * Getter for the maximally allowed time step for this particle
   * @return the maximally allowed time step for this particle
   */
  double getDt() const { return _dt; }

  /**
   * Set the maximally allowed time step for this particle
   * @param dt the maximally allowed time step for this particle
   */
  void setDt(double dt) { _dt = dt; }

  /**
   * Calculate the maximally allowed time step for the particle based on the
   * smoothing length and the signal velocity of the particle
   */
  void calcDt() {
    const double C_CFL = 0.3;
    _dt = C_CFL * 2.0 * _smth / _v_sig_max;
  }

  /**
   * Getter for velocity at half-time step (leapfrog)
   * @return
   */
  const std::array<CalcType, 3> &getVel_half() const { return _vel_half; }

  /**
   * Setter for velocity at half-time step (leapfrog)
   * @param vel_half
   */
  void setVel_half(const std::array<CalcType, 3> &vel_half) { SPHParticleBase::_vel_half = vel_half; }

  /**
   * Getter for energy at half-time step (leapfrog)
   * @return
   */
  double getEng_half() const { return _eng_half; }

  /**
   * Setter for energy at half-time step (leapfrog)
   * @param eng_half
   */
  void setEng_half(double eng_half) { SPHParticleBase::_eng_half = eng_half; }

  /**
   * function to serialize an SPHParticleBase
   * @return serialized vector of bytes (char)
   */
  std::vector<double> serialize() const {
    std::vector<double> stream;
    for (int i = 0; i < 3; i++) {
      stream.push_back(this->getR()[i]);
    }
    for (int i = 0; i < 3; i++) {
      stream.push_back(this->getV()[i]);
    }

    for (int i = 0; i < 3; i++) {
      // stream.push_back(this->getF()[i]);  // not actually needed
    }
    auto id = this->getID();
    double id_dbl;
    memcpy(&id_dbl, &id, sizeof(double));
    static_assert(sizeof(id) == sizeof(double), "sizes should be the same, otherwise the above will not work");

    stream.push_back(id_dbl);
    stream.push_back(_density);
    stream.push_back(_pressure);
    stream.push_back(_mass);
    stream.push_back(_smth);
    stream.push_back(_snds);

    for (int i = 0; i < 3; i++) {
      stream.push_back(this->getAcceleration()[i]);
    }
    stream.push_back(_energy_dot);
    stream.push_back(_energy);
    // stream.push_back(_dt); // not needed
    for (int i = 0; i < 3; i++) {
      stream.push_back(this->getVel_half()[i]);
    }
    stream.push_back(_eng_half);
    return stream;
  }

  /**
   * funtion to deserialize an SPHParticleBase
   * @param stream
   * @param index start index within the stream, will be increased while
   * deserializing to mark already processed data.
   * @return
   */
  static SPHParticleBase deserialize(const double *stream, size_t &index) {
    std::array<double, 3> r = {stream[index], stream[index + 1], stream[index + 2]};
    index += 3;
    std::array<double, 3> v = {stream[index], stream[index + 1], stream[index + 2]};
    index += 3;
    // std::array<double,3> F = {stream[index], stream[index+1],
    // stream[index+2]};  // not needed
    // index+=3  // not needed
    double id_dbl = stream[index++];
    unsigned long id;
    memcpy(&id, &id_dbl, sizeof(double));
    static_assert(sizeof(id) == sizeof(double), "sizes should be the same, otherwise the above will not work");

    double density = stream[index++];
    double pressure = stream[index++];

    double mass = stream[index++];
    double smth = stream[index++];
    double snds = stream[index++];

    std::array<double, 3> ac = {stream[index], stream[index + 1], stream[index + 2]};
    index += 3;
    double energy_dot = stream[index++];
    double energy = stream[index++];
    // double dt = stream[index++];  // not needed
    std::array<double, 3> vel_half = {stream[index], stream[index + 1], stream[index + 2]};
    index += 3;
    double eng_half = stream[index++];

    SPHParticleBase p = SPHParticleBase(r, v, id, mass, smth, snds);
    p.setDensity(density);
    p.setPressure(pressure);
    p.setAcceleration(ac);
    p.setEngDot(energy_dot);
    p.setEnergy(energy);
    p.setVel_half(vel_half);
    p.setEng_half(eng_half);
    return p;
  }

  /**
   * Attribute names for the soa arrays
   */
  enum AttributeNames : int {
    ptr,
    mass,
    posX,
    posY,
    posZ,
    smth,
    density,
    velX,
    velY,
    velZ,
    soundSpeed,
    pressure,
    vsigmax,
    accX,
    accY,
    accZ,
    engDot,
    ownershipState
  };

  /**
   * SoA arrays type, cf. AttributeNames
   */
  using SoAArraysType = typename autopas::utils::SoAType<SPHParticleBase<CalcType, AccuType, idType> *,
                                                         CalcType,  // mass
                                                         CalcType,  // posX
                                                         CalcType,  // posY
                                                         CalcType,  // posZ
                                                         CalcType,  // smth
                                                         CalcType,  // density
                                                         CalcType,  // velX
                                                         CalcType,  // velY
                                                         CalcType,  // velZ
                                                         CalcType,  // soundSpeed
                                                         CalcType,  // pressure
                                                         CalcType,  // vsigmax
                                                         AccuType,  // accX
                                                         AccuType,  // accY
                                                         AccuType,  // accZ
                                                         AccuType,  // engDot
                                                         autopas::OwnershipState>::Type;

  /**
   * Non-const getter for the pointer of this object.
   * @tparam attribute Attribute name.
   * @return this.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::mass) {
      return this->getMass();
    } else if constexpr (attribute == AttributeNames::posX) {
      return this->getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return this->getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return this->getR()[2];
    } else if constexpr (attribute == AttributeNames::smth) {
      return this->getSmoothingLength();
    } else if constexpr (attribute == AttributeNames::density) {
      return this->getDensity();
    } else if constexpr (attribute == AttributeNames::velX) {
      return this->getV()[0];
    } else if constexpr (attribute == AttributeNames::velY) {
      return this->getV()[1];
    } else if constexpr (attribute == AttributeNames::velZ) {
      return this->getV()[2];
    } else if constexpr (attribute == AttributeNames::soundSpeed) {
      return this->getSoundSpeed();
    } else if constexpr (attribute == AttributeNames::pressure) {
      return this->getPressure();
    } else if constexpr (attribute == AttributeNames::vsigmax) {
      return this->getVSigMax();
    } else if constexpr (attribute == AttributeNames::accX) {
      return this->getAcceleration()[0];
    } else if constexpr (attribute == AttributeNames::accY) {
      return this->getAcceleration()[1];
    } else if constexpr (attribute == AttributeNames::accZ) {
      return this->getAcceleration()[2];
    } else if constexpr (attribute == AttributeNames::engDot) {
      return this->getEngDot();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("SPHParticleBase::get: unknown attribute");
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::mass) {
      this->setMass(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      this->_r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      this->_r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      this->_r[2] = value;
    } else if constexpr (attribute == AttributeNames::smth) {
      this->setSmoothingLength(value);
    } else if constexpr (attribute == AttributeNames::density) {
      this->setDensity(value);
    } else if constexpr (attribute == AttributeNames::velX) {
      this->_v[0] = value;
    } else if constexpr (attribute == AttributeNames::velY) {
      this->_v[1] = value;
    } else if constexpr (attribute == AttributeNames::velZ) {
      this->_v[2] = value;
    } else if constexpr (attribute == AttributeNames::soundSpeed) {
      this->setSoundSpeed(value);
    } else if constexpr (attribute == AttributeNames::pressure) {
      this->setPressure(value);
    } else if constexpr (attribute == AttributeNames::vsigmax) {
      this->setVSigMax(value);
    } else if constexpr (attribute == AttributeNames::accX) {
      this->_acc[0] = value;
    } else if constexpr (attribute == AttributeNames::accY) {
      this->_acc[1] = value;
    } else if constexpr (attribute == AttributeNames::accZ) {
      this->_acc[2] = value;
    } else if constexpr (attribute == AttributeNames::engDot) {
      this->setEngDot(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("SPHParticleBase::set: unknown attribute");
    }
  }

 private:
  CalcType _density;   // density
  CalcType _pressure;  // pressure
  CalcType _mass;      // mass
  CalcType _smth;      // smoothing length
  CalcType _snds;      // speed of sound

  // temporaries / helpers
  CalcType _v_sig_max;

  // integrator
  std::array<AccuType, 3> _acc;  // acceleration
  AccuType _energy_dot;          // time derivative of the energy
  AccuType _energy;              // energy
  AccuType _dt;                  // local timestep allowed by this particle

  // integrator
  std::array<CalcType, 3> _vel_half;  // velocity at half time-step
  CalcType _eng_half;                 // energy at half time-step
};
}  // namespace sphLib