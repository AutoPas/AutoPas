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

  virtual ~MoleculeLJ() = default;

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
   * @return Mass
   */
  static double getMass() { return MASS;}

    /**set mass
   * @param Mass
   */
  static void setMass(double mass) {MASS = mass;}

    /**
     * the type for the soa storage
     */
  typedef autopas::utils::SoAType<size_t, double, double, double, double, double, double>::Type SoAArraysType;

    /**get OldForce
    * @return OldForce
    */
  static double getOldf() {return OLDF;}

    /**set OldForce
    * @param OldForce
    */
  static void setOldf(double oldf) { OLDF = oldf;}

private:
  static double EPSILON, SIGMA, MASS, OLDF;
};

}  // namespace autopas