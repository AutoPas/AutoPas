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

/**
 * This class stores the (physical) properties of particle types.
 *
 * It also provides mixed values for (force) calculations between known types.
 */
class ParticlePropertiesLibrary {
 public:
  /**
   * Default constructor.
   */
  ParticlePropertiesLibrary() = default;

  /**
   * Copy Constructor.
   * @param particlePropertiesLibrary
   */
  ParticlePropertiesLibrary(const ParticlePropertiesLibrary &particlePropertiesLibrary);

  /**
   * Copy assignment operator.
   * @param particlePropertiesLibrary
   * @return
   */
  ParticlePropertiesLibrary &operator=(const ParticlePropertiesLibrary &particlePropertiesLibrary);

  /**
   * Adds the properties of a particle type to the library.
   *
   * This function also precomputes all possible mixed values with already known particle types.
   * If the type id already exists the values will be overwritten.
   * @param typeID
   * @param epsilon
   * @param sigma
   * @param mass
   */
  void addType(unsigned long typeID, double epsilon, double sigma, double mass);

  ~ParticlePropertiesLibrary() = default;

  /**
   * Getter for the particle's epsilon*24.
   * @param i typeId of the particle.
   * @return 24*epsilon_i
   */
  double get24Epsilon(unsigned long i);

  /**
   * Getter for the particle's squared sigma.
   * @param i typeId of the particle.
   * @return sigma_i²
   */
  double getSigmaSquare(unsigned long i);

  /**
   * Getter for the particle's mass.
   * @param i typeId of the particle.
   * @return mass_i
   */
  double getMass(unsigned long i);

  /**
   * Returns the precomputed mixed epsilon24.
   * @param  i typeId index of particle one.
   * @param  j typeId index of particle two.
   * @return 24*epsilon_ij
   * */
  inline double mixing24Epsilon(unsigned long i, unsigned long j) const {
    auto key = std::make_pair((i < j) ? i : j, (j > i) ? j : i);  // key in preprocessed maps: (i,j) with i<j
    return _computedMixing24Epsilon.at(key);
  }
  /**
   * Returns precomputed mixed squared sigma.
   * @param i typeId index of particle one.
   * @param j typeId index of particle two.
   * @return sigma_ij²
   */
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
