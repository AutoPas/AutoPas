/**
 * @file ParticleClassLibrary.h
 * @author N. Fottner
 * @date 7/4/19
 */

#pragma once
#include <cmath>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"

/**This Class is used to map the Particles to their Physical Properties used in the
 * Force calculations(epsilon,sigma) and in the timeDiscretization (mass)
 * The mapping key is the type of the Particles
 * for the lennard jones force Calculation this class provides the preprocessed values
 * of the mixing rules for epsilon and sigma values
 * */
class ParticlePropertiesLibrary {
 public:
  ParticlePropertiesLibrary() = default;

  ParticlePropertiesLibrary(double &epsilon, double &sigma, double mass);

  ParticlePropertiesLibrary(const ParticlePropertiesLibrary &pcl);

  ParticlePropertiesLibrary &operator=(const ParticlePropertiesLibrary &plc);

  /**adds the Properties of a particle type to the class members
   * calculates and adds the values of the mixing rules between the added particle type and the existing ones to the
   * appropriate members
   * @param typeID
   * @param epsilon
   * @param sigma
   * @param mass
   * */
  void addType(unsigned long typeID, double epsilon, double sigma, double mass);

  ~ParticlePropertiesLibrary() = default;
  /**Getter for Particle Epsilon*24
   * @param Particle
   * @return Epsilon*24
   */
  double get24Epsilon(unsigned long i);
  /**Getter for Particle Square Sigma
   * @param Particle
   * @return Sigma²
   */
  double getSigmaSquare(unsigned long i);

  /**Getter for Particle Mass
   * @param Particle
   * @return Sigma
   */
  double getMass(unsigned long i);

  /**Returns (Epsilon*24) of the MixingRule of 2 Particles precalculated in computedMixing24Epsilon
   * @param Particles; i and j
   * @return 24*(Epsilon of both)
   * */
  double mixing24Epsilon(unsigned long i, unsigned long j);
  /**Returns Sigma Square of the MixingRule of 2 Particles precalculated in computedMixingSigmaSquare
   * @param Particles; i and j
   * @return (Sigma of both)²
   * */
  double mixingSigmaSquare(unsigned long i, unsigned long j);

 private:
  std::map<unsigned long, double> Epsilon;
  std::map<unsigned long, double> Sigma;
  std::map<unsigned long, double> Mass;
  std::map<std::pair<unsigned long, unsigned long>, double> computedMixing24Epsilon;
  std::map<std::pair<unsigned long, unsigned long>, double> computedMixingSigmaSquare;
};