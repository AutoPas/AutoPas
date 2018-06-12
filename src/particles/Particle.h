/**
 * @file Particle.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once
#include <array>
#include <tuple>
#include "utils/SoAType.h"
#include "utils/ArrayMath.h"
#include "utils/SoAStorage.h"

namespace autopas {

/**
 * basic Particle class
 * This class can be used to build your own Particle class. However, you are
 * free to not use it as well.
 */
class Particle {
 public:
  Particle() : _r({0.0, 0.0, 0.0}), _v({0., 0., 0.}), _f({0.0, 0.0, 0.0}), _id(0) {}

  /**
   * Constructor of the Particle class
   * @param r position of the particle
   * @param v velocity of the particle
   * @param id id of the particle
   */
  Particle(std::array<double, 3> r, std::array<double, 3> v, unsigned long id)
      : _r(r), _v(v), _f({0.0, 0.0, 0.0}), _id(id) {}

  virtual ~Particle() = default;

  /**
   * get the force acting on the particle
   * @return force
   */
  const std::array<double, 3> &getF() const { return _f; }

  /**
   * Set the force acting on the particle
   * @param f force
   */
  void setF(const std::array<double, 3> &f) { _f = f; }

  /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
   */
  void addF(const std::array<double, 3> &f) { _f = ArrayMath::add(_f, f); }

  /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
   */
  void subF(const std::array<double, 3> &f) { _f = ArrayMath::sub(_f, f); }

  /**
   * Get the id of the particle
   * @return id
   */
  unsigned long getID() const { return _id; }

  /**
   * Set the id of the particle
   * @param id id
   */
  void setID(unsigned long id) { _id = id; }

  /**
   * Get the position of the particle
   * @return current position
   */
  const std::array<double, 3> &getR() const { return _r; }

  /**
   * Set the position of the particle
   * @param r new position
   */
  void setR(const std::array<double, 3> &r) { _r = r; }

  /**
   * Add a distance vector to the position of the particle
   * @param r vector to be added
   */
  void addR(const std::array<double, 3> &r) { _r = ArrayMath::add(_r, r); }

  /**
   * Checks whether the particle is within a cuboidal box specified by rmin and
   * rmax
   * @param rmin lower corner of the box
   * @param rmax higher corner of the box
   * @return true if the particle is in the box, false otherwise
   */
  bool inBox(const std::array<double, 3> &rmin, const std::array<double, 3> rmax) {
    bool in = true;
    for (int d = 0; d < 3; ++d) {
      in &= (_r[d] >= rmin[d] and _r[d] < rmax[d]);
    }
    return in;
  }

  /**
   * Get the velocity of the particle
   * @return current velocity
   */
  const std::array<double, 3> &getV() const { return _v; }

  /**
   * Set the velocity of the particle
   * @param v new velocity
   */
  void setV(const std::array<double, 3> &v) { _v = v; }

  /**
   * Add a vector to the current velocity of the particle
   * @param v vector to be added
   */
  void addV(const std::array<double, 3> &v) { _v = ArrayMath::add(_v, v); }

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { id, posX, posY, posZ, forceX, forceY, forceZ };

  /**
   * the type for the soa storage
   */
  typedef autopas::utils::SoAType<size_t, double, double, double, double, double, double>::Type SoAArraysType;

 private:
  std::array<double, 3> _r;
  std::array<double, 3> _v;
  std::array<double, 3> _f;
  unsigned long _id;
};

}  // namespace autopas
