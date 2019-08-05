/**
 * @file ParticlePropertiesLibrary.cpp
 * @author N. Fottner
 * @date 7/4/19
 */

#include "ParticlePropertiesLibrary.h"

void ParticlePropertiesLibrary::addType(unsigned long typeID, double epsilon, double sigma, double mass) {
  _epsilons.emplace(typeID, epsilon);
  for (auto &e : _epsilons) {
    unsigned long indexOfExistingEpsilon = std::get<0>(e);
    double secondEpsilon = std::get<1>(e);
    double epsilon24 = 24 * sqrt(epsilon * secondEpsilon);
    auto newEntry = std::make_pair(indexOfExistingEpsilon, typeID);
    _computedMixing24Epsilon.emplace(newEntry, epsilon24);
  }

  _sigmas.emplace(typeID, sigma);
  for (auto &s : _sigmas) {
    unsigned long indexOfExistingSigma = std::get<0>(s);
    double existingSigma = std::get<1>(s);
    double newSigma = (sigma + existingSigma) / 2;
    auto newEntry = std::make_pair(indexOfExistingSigma, typeID);
    _computedMixingSigmaSquare.emplace(newEntry, (newSigma * newSigma));
  }

  _masses.emplace(typeID, mass);
}

ParticlePropertiesLibrary::ParticlePropertiesLibrary(const ParticlePropertiesLibrary &particlePropertiesLibrary) = default;

ParticlePropertiesLibrary &ParticlePropertiesLibrary::operator=(const ParticlePropertiesLibrary &particlePropertiesLibrary) = default;

double ParticlePropertiesLibrary::getMass(unsigned long i) { return _masses.at(i); }

double ParticlePropertiesLibrary::get24Epsilon(unsigned long i) {
  return _computedMixing24Epsilon[std::make_pair(i, i)];
}

double ParticlePropertiesLibrary::getSigmaSquare(unsigned long i) {
  return _computedMixingSigmaSquare[std::make_pair(i, i)];
}
