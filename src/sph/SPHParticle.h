//
// Created by seckler on 19.01.18.
//

#ifndef AUTOPAS_SPHPARTICLE_H
#define AUTOPAS_SPHPARTICLE_H

#include "particles/Particle.h"

namespace autopas {
namespace sph {
/**
 * Basic SPHParticle class.
 */
class SPHParticle : public autopas::Particle {
 public:
  /**
   * Default constructor of SPHParticle.
   * Will initialize all values to some basic defaults.
   */
  SPHParticle()
      : autopas::Particle(),
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
        _dt(0.) {}
  /**
  * Constructor of the SPHParticle class.
  * @param r position of the particle
  * @param v velocity of the particle
  * @param id id of the particle. This id should be unique
  */
  SPHParticle(std::array<double, 3> r, std::array<double, 3> v,
              unsigned long id)
      : autopas::Particle(r, v, id),
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
        _dt(0.) {}

  /**
   * Constructor of the SPHParticle class.
   * @param r position of the particle
   * @param v velocity of the particle
   * @param id id of the particle. This id should be unique
   * @param mass mass of the particle
   * @param smth smoothing length of the particle
   * @param snds speed of sound (SouND Speed)
   */
  SPHParticle(std::array<double, 3> r, std::array<double, 3> v,
              unsigned long id, double mass, double smth, double snds)
      : autopas::Particle(r, v, id),
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
        _dt(0.) {}

  /**
   * Destructor of the SPHParticle
   */
  ~SPHParticle() override = default;

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
  void calcPressure();

  /**
   * Setter for the pressure
   * @param pressure pressure value to be set
   */
  void setPressure(double pressure) { _pressure = pressure; }

  /**
   * Getter for the mass of the particle
   * @return mass of particle
   */
  double getMass() const { return _mass; }

  /**
   * Setter for the mass of the particle
   * @param mass mass to be set
   */
  void setMass(double mass) { _mass = mass; }

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
  void checkAndSetVSigMax(double v_sig) {
    _v_sig_max = std::max(v_sig, _v_sig_max);
  }

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
  void addAcceleration(const std::array<double, 3> &acc);

  /**
   * Substracts the given acceleration from the local acceleration.
   * Used to sum up different negative acceleration values.
   * @param acc Acceleration to be substracted
   */
  void subAcceleration(const std::array<double, 3> &acc);

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

 private:
  double _density;
  double _pressure;
  double _mass;
  double _smth;
  double _snds;

  // temporaries / helpers
  double _v_sig_max;
  std::array<double, 3> _acc;
  double _energy_dot;
  double _energy;
  double _dt;
};
}  // namespace sph
}  // namespace autopas
#endif  // AUTOPAS_SPHPARTICLE_H
