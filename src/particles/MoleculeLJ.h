/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "Particle.h"

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

 private:
  static double EPSILON, SIGMA;
};

}  // namespace autopas