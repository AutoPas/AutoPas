/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>
#include "autopas/particles/Particle.h"

namespace autopas {

/**
 * lennard jones molecule class
 */
class MoleculeLJ : public Particle {
 public:
  MoleculeLJ() = default;

  /**
   * constructor of a lennard jones molecule
   * @param r position of the molecule
   * @param v velocity of the molecule
   * @param id id of the molecule
   */
  explicit MoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : Particle(r, v, id) {}

  ~MoleculeLJ() override = default;

  /**
   * get epsilon (characteristic energy of the lj potential)
   * @return epsilon
   */
  static double getEpsilon() { return EPSILON; }

  /**
   * set epsilon (characteristic energy of the lj potential)
   * @param epsilon
   */
  static void setEpsilon(double epsilon) { EPSILON = epsilon; }

  /**
   * get sigma (characteristic length of the lj potential)
   * @return sigma
   */
  static double getSigma() { return SIGMA; }

  /**
   * set sigma (characteristic length of the lj potential)
   * @param sigma
   */
  static void setSigma(double sigma) { SIGMA = sigma; }

  /**get mass
   * @return MASS
   */
  static double getMass() { return MASS; }

  /**set mass
   * @param mass
   */
  static void setMass(double mass) { MASS = mass; }

  /**
   * the type for the soa storage
   */
  //  typedef autopas::utils::SoAType<size_t, double, double, double, double, double, double>::Type SoAArraysType;

  /**get OldForce
   * @return OLDF
   */
  std::array<double, 3> getOldf() const { return OLDF; }

  /**set OldForce
   * @param oldf
   */
  void setOldf(const std::array<double, 3> &oldf) { OLDF = oldf; }

  /**get TypeId
   * @return _typeId
   * */
  size_t getTypeId() const { return _typeId; }
  /**set _TypeId of Particle
   * @param typeId
   * */
  void setTypeId(size_t typeId) { _typeId = typeId; }

 private:
  static double EPSILON, SIGMA, MASS;

  /**
   * Particle type id.
   */
  size_t _typeId = 0;

 private:
  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> OLDF = {0., 0., 0.};
};

}  // namespace autopas