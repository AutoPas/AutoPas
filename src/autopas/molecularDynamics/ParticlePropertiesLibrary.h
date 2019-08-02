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
   * @param pcl
   * */
  ParticlePropertiesLibrary(const ParticlePropertiesLibrary &pcl);
  /**Copy assignment operator
   * @param plc
   * @return ParticlePropertiesLibrary
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
   * @param  i (typeId index)
   * @param  j (typeId index)
   * @return 24*(epsilonMixingRule)
   * */
  inline double mixing24Epsilon(unsigned long i, unsigned long j) const {
    auto key = std::make_pair((i < j) ? i : j, (j > i) ? j : i);  // key in preprocessed maps: (i,j) with i<j
    return _computedMixing24Epsilon.at(key);
  }
  /**Returns Sigma Square of the MixingRule of 2 Particles precalculated in computedMixingSigmaSquare
   * @param i (typeId index)
   * @param j (typeId index)
   * @return (sigmaMixingRule)²
   * */
  inline double mixingSigmaSquare(unsigned long i, unsigned long j) const {
    auto key = std::make_pair((i < j) ? i : j, (j > i) ? j : i);  // key in preprocessed maps: (i,j) with i<j
    return _computedMixingSigmaSquare.at(key);
  }

 private:
  std::map<unsigned long, double> _epsilons;
  std::map<unsigned long, double> _sigmas;
  std::map<unsigned long, double> _masses;
  std::map<std::pair<unsigned long, unsigned long>, double> _computedMixing24Epsilon;
  std::map<std::pair<unsigned long, unsigned long>, double> _computedMixingSigmaSquare;
};
