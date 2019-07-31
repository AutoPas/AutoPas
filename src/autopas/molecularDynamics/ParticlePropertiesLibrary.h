/**
 * @file ParticlePropertiesLibrary.h
 * @author N. Fottner
 * @date 7/4/19
 */

#pragma once
#include <cmath>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"

//@todo soon: add template parameter to support 32 or 64 bit values-> did resolve into undefined references

/**This Class is used to map the Particles to their Physical Properties used in the
 * Force calculations(epsilon,sigma) and in the timeDiscretization (mass)
 * The mapping key is the type of the Particles
 * for the lennard jones force Calculation this class provides the preprocessed values
 * of the mixing rules for epsilon and sigma values
 * */
class ParticlePropertiesLibrary {
 public:
  ParticlePropertiesLibrary() = default;
  /**Constructor if there is only one Type used
   * @param epsilon
   * @param sigma
   * @param mass
   * */
  ParticlePropertiesLibrary(double &epsilon, double &sigma, double mass);
  /**Copy Constructor
   * @param ParticlePropertiesLibrary
   * */
  ParticlePropertiesLibrary(const ParticlePropertiesLibrary &pcl);
  /**Copy assignment operator
   * @param ParticlePropertiesLibrary
   * */
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
   * @param i (typeId of Particle)
   * @return Epsilon*24
   */
  double get24Epsilon(unsigned long i);
  /**Getter for Particle Square Sigma
   * @param i (typeId of Particle)
   * @return Sigma²
   */
  double getSigmaSquare(unsigned long i);

  /**Getter for Particle Mass
   * @param i (typeId of Particle)
   * @return Sigma
   */
  double getMass(unsigned long i);

  /**Returns (Epsilon*24) of the MixingRule of 2 Particles precalculated in computedMixing24Epsilon
   * @param  i and j (typeId index)
   * @return 24*(epsilonMixingRule)
   * */
  double mixing24Epsilon(unsigned long i, unsigned long j);
  /**Returns Sigma Square of the MixingRule of 2 Particles precalculated in computedMixingSigmaSquare
   * @param i and j (typeId index)
   * @return (sigmaMixingRule)²
   * */
  double mixingSigmaSquare(unsigned long i, unsigned long j);

 private:
  std::map<unsigned long, double> Epsilon;
  std::map<unsigned long, double> Sigma;
  std::map<unsigned long, double> Mass;
  std::map<std::pair<unsigned long, unsigned long>, double> computedMixing24Epsilon;
  std::map<std::pair<unsigned long, unsigned long>, double> computedMixingSigmaSquare;
};
