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

/**
 * This class stores the (physical) properties of particle types.
 *
 * It also provides mixed values for (force) calculations between known types.
 */
template <typename floatType = double, typename intType = unsigned long>
class ParticlePropertiesLibrary {
 public:
  /**
   * Type for floating point numbers.
   */
  using ParticlePropertiesLibraryFloatType = floatType;

  /**
   * Type for integer numbers.
   */
  using ParticlePropertiesLibraryIntType = intType;

  /**
   * Default constructor.
   */
  ParticlePropertiesLibrary() = default;

  /**
   * Copy Constructor.
   * @param particlePropertiesLibrary
   */
  ParticlePropertiesLibrary(const ParticlePropertiesLibrary &particlePropertiesLibrary) = default;

  /**
   * Copy assignment operator.
   * @param particlePropertiesLibrary
   * @return
   */
  ParticlePropertiesLibrary &operator=(const ParticlePropertiesLibrary &particlePropertiesLibrary) = default;

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
  void addType(intType typeID, floatType epsilon, floatType sigma, floatType mass);

  ~ParticlePropertiesLibrary() = default;

  /**
   * Getter for the particle's epsilon*24.
   * @param i typeId of the particle.
   * @return 24*epsilon_i
   */
  floatType get24Epsilon(intType i) const;

  /**
   * Getter for the particle's squared sigma.
   * @param i typeId of the particle.
   * @return sigma_i²
   */
  floatType getSigmaSquare(intType i) const;

  /**
   * Getter for the particle's mass.
   * @param i typeId of the particle.
   * @return mass_i
   */
  floatType getMass(intType i) const;

  /**
   * Returns the precomputed mixed epsilon24.
   * @param  i typeId index of particle one.
   * @param  j typeId index of particle two.
   * @return 24*epsilon_ij
   */
  inline floatType mixing24Epsilon(intType i, intType j) const {
    auto key = std::make_pair((i < j) ? i : j, (j > i) ? j : i);  // key in preprocessed maps: (i,j) with i<j
    return _computedMixing24Epsilon.at(key);
  }
  /**
   * Returns precomputed mixed squared sigma.
   * @param i typeId index of particle one.
   * @param j typeId index of particle two.
   * @return sigma_ij²
   */
  inline floatType mixingSigmaSquare(intType i, intType j) const {
    auto key = std::make_pair((i < j) ? i : j, (j > i) ? j : i);  // key in preprocessed maps: (i,j) with i<j
    return _computedMixingSigmaSquare.at(key);
  }

 private:
  std::map<intType, floatType> _epsilons;
  std::map<intType, floatType> _sigmas;
  std::map<intType, floatType> _masses;
  std::map<std::pair<intType, intType>, floatType> _computedMixing24Epsilon;
  std::map<std::pair<intType, intType>, floatType> _computedMixingSigmaSquare;
};

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addType(intType typeID, floatType epsilon, floatType sigma,
                                                            floatType mass) {
  _epsilons.emplace(typeID, epsilon);
  for (auto &[indexOfExistingEpsilon, secondEpsilon] : _epsilons) {
    floatType epsilon24 = 24 * sqrt(epsilon * secondEpsilon);
    auto newEntry = std::make_pair(indexOfExistingEpsilon, typeID);
    _computedMixing24Epsilon.emplace(newEntry, epsilon24);
  }

  _sigmas.emplace(typeID, sigma);
  for (auto &[indexOfExistingSigma, existingSigma] : _sigmas) {
    floatType newSigma = (sigma + existingSigma) / 2.0;
    auto newEntry = std::make_pair(indexOfExistingSigma, typeID);
    _computedMixingSigmaSquare.emplace(newEntry, (newSigma * newSigma));
  }

  _masses.emplace(typeID, mass);
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getMass(intType i) const {
  return _masses.at(i);
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::get24Epsilon(intType i) const {
  return _computedMixing24Epsilon.at(std::make_pair(i, i));
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getSigmaSquare(intType i) const {
  return _computedMixingSigmaSquare.at(std::make_pair(i, i));
}
