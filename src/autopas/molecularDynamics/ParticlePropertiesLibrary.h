/**
 * @file ParticlePropertiesLibrary.h
 * @author N. Fottner
 * @date 7/4/19
 */

#pragma once

#include <autopas/utils/ExceptionHandler.h>

#include <cmath>
#include <map>
#include <set>
#include <vector>

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
   * Constructor
   * @param cutoff Cutoff for the Lennard Jones Potential (needed for calculation of shift)
   */
  explicit ParticlePropertiesLibrary(const double cutoff) : _cutoff(cutoff) {}

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

  void calculateMixingCoefficients();

  ~ParticlePropertiesLibrary() = default;

  /**
   * Returns a set of all particle types stored.
   * @return
   */
  std::set<intType> getTypes() const {
    std::set<intType> typeIDs;

    for (size_t index = 0; index < _numRegisteredTypes; ++index) {
      typeIDs.emplace(index);
    }

    return typeIDs;
  }

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
   * @param  i typeId of particle one.
   * @param  j typeId of particle two.
   * @return 24*epsilon_ij
   */
  inline floatType mixing24Epsilon(intType i, intType j) const {
    return _computedMixing24Epsilon.at(i * _numRegisteredTypes + j);
  }

  /**
   * Returns precomputed mixed squared sigma.
   * @param i typeId of particle one.
   * @param j typeId of particle two.
   * @return sigma_ij²
   */
  inline floatType mixingSigmaSquare(intType i, intType j) const {
    return _computedMixingSigmaSquare.at(i * _numRegisteredTypes + j);
  }

  /**
   * Returns precomputed mixed shift 6.
   * @param i typeId of particle one.
   * @param j typeId of particle two.
   * @return shift * 6
   */
  inline floatType mixingShift6(intType i, intType j) const {
    return _computedMixingShift6.at(i * _numRegisteredTypes + j);
  }

  /**
   * Calculate the shift of the lennard jones potential from given cutoff, epsilon, sigma.
   * @param epsilon24
   * @param sigmaSquare
   * @param cutoffSquare squared cutoff of the potential that should be shifted
   * @return shift multiplied by 6
   */
  static double calcShift6(double epsilon24, double sigmaSquare, double cutoffSquare);

 private:
  intType _numRegisteredTypes{0};
  const double _cutoff;

  std::vector<floatType> _epsilons;
  std::vector<floatType> _sigmas;
  std::vector<floatType> _masses;

  std::vector<floatType> _computedMixing24Epsilon;
  std::vector<floatType> _computedMixingSigmaSquare;
  std::vector<floatType> _computedMixingShift6;
};

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::addType(intType typeID, floatType epsilon, floatType sigma,
                                                            floatType mass) {
  if (_numRegisteredTypes != typeID) {
    autopas::utils::ExceptionHandler::exception(
        "ParticlePropertiesLibrary::addType(): trying to register a type with id {}. Please register types "
        "consecutively, starting at id 0. Currently there are {} registered types.",
        typeID, _numRegisteredTypes);
  }
  ++_numRegisteredTypes;
  _epsilons.emplace_back(epsilon);
  _sigmas.emplace_back(sigma);
  _masses.emplace_back(mass);
}

template <typename floatType, typename intType>
void ParticlePropertiesLibrary<floatType, intType>::calculateMixingCoefficients() {
  _computedMixing24Epsilon.resize(_numRegisteredTypes * _numRegisteredTypes);
  _computedMixingSigmaSquare.resize(_numRegisteredTypes * _numRegisteredTypes);
  _computedMixingShift6.resize(_numRegisteredTypes * _numRegisteredTypes);

  auto cutoffSquare = _cutoff * _cutoff;

  for (size_t firstIndex = 0ul; firstIndex < _numRegisteredTypes; ++firstIndex) {
    for (size_t secondIndex = 0ul; secondIndex < _numRegisteredTypes; ++secondIndex) {
      auto globalIndex = _numRegisteredTypes * firstIndex + secondIndex;

      // epsilon
      floatType epsilon24 = 24 * sqrt(_epsilons[firstIndex] * _epsilons[secondIndex]);
      _computedMixing24Epsilon[globalIndex] = epsilon24;

      // sigma
      floatType sigma = (_sigmas[firstIndex] + _sigmas[secondIndex]) / 2.0;
      floatType sigmaSquare = sigma * sigma;
      _computedMixingSigmaSquare[globalIndex] = sigmaSquare;

      // shift6
      floatType newShift6 = calcShift6(epsilon24, sigmaSquare, cutoffSquare);
      _computedMixingShift6[globalIndex] = newShift6;
    }
  }
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getMass(intType i) const {
  return _masses.at(i);
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::get24Epsilon(intType i) const {
  return _computedMixing24Epsilon.at(i * _numRegisteredTypes + i);
}

template <typename floatType, typename intType>
floatType ParticlePropertiesLibrary<floatType, intType>::getSigmaSquare(intType i) const {
  return _computedMixingSigmaSquare.at(i * _numRegisteredTypes + i);
}

template <typename floatType, typename intType>
double ParticlePropertiesLibrary<floatType, intType>::calcShift6(double epsilon24, double sigmaSquare,
                                                                 double cutoffSquare) {
  auto sigmaPow2DivCutoff = sigmaSquare / (cutoffSquare);
  auto sigmaPow6DivCutoff = sigmaPow2DivCutoff * sigmaPow2DivCutoff * sigmaPow2DivCutoff;
  auto shift6 = epsilon24 * (sigmaPow6DivCutoff - sigmaPow6DivCutoff * sigmaPow6DivCutoff);
  return shift6;
}
