#include "ParticleClassLibrary.h"

    /** Default constuktor wenn es nur ein Particle Type gibt
    **/
    ParticleClassLibrary::ParticleClassLibrary(double &epsilon, double &sigma, double mass) {
  std::map<unsigned long, double> EMap;
  std::map<unsigned long, double> SMap;
  std::map<unsigned long, double> MMap;
    EMap.emplace(0, epsilon);
    SMap.emplace(0, sigma);
    MMap.emplace(0, mass);
  this->Epsilon = EMap;
  this->Sigma = SMap;
  this->Mass = MMap;
  this->computedMixing24Epsilon.emplace(std::make_pair(0,0),24*epsilon);
  this->computedMixingSigmaSquare.emplace(std::make_pair(0,0),(sigma*sigma));
    }
//@todo need to initilized addType with ascending typeId
//@todo when initializing, and there is only one type -> use appropriate constructor
void ParticleClassLibrary::addType(unsigned long typeID, double epsilon, double sigma,
                                                   double mass) {
    for (auto e: Epsilon) {
        unsigned long indexOfExistingEpsilon = std::get<0>(e);
        double secondEpsilon = std::get<1>(e);
        double epsilon24 = 24 * sqrt(epsilon * secondEpsilon);
        auto newEntry= std::make_pair(indexOfExistingEpsilon, typeID);
        computedMixing24Epsilon.emplace(newEntry,epsilon24);
    }
        Epsilon.emplace(typeID, epsilon);
    for (auto e: Sigma) {
        unsigned long indexOfExistingSigma = std::get<0>(e);
        double existingSigma = std::get<1>(e);
        double newSigma = (sigma + existingSigma)/2;
        auto newEntry= std::make_pair(indexOfExistingSigma, typeID);
        computedMixingSigmaSquare.emplace(newEntry, (newSigma*newSigma));
    }
        Sigma.emplace(typeID, sigma);
        Mass.emplace(typeID, mass);
}


ParticleClassLibrary::ParticleClassLibrary(const ParticleClassLibrary &pcl) = default;


ParticleClassLibrary &ParticleClassLibrary::operator=(const ParticleClassLibrary &pcl) = default;

double ParticleClassLibrary::getMass(unsigned long i) {
return Mass.at(i);
}
double ParticleClassLibrary::get24Epsilon(unsigned long i) { return 24* Epsilon.at(i); }

double ParticleClassLibrary::getSigmaSquare(unsigned long i) {
  return (Sigma.at(i)*Sigma.at(i));
}

double ParticleClassLibrary::mixing24Epsilon(unsigned long i, unsigned long j) {
  auto key = std::make_pair((i<j)?i:j,(j>i)?j:i); //key in preprocessed maps: (i,j) with i<j
  return computedMixing24Epsilon.at(key);
}

double ParticleClassLibrary::mixingSigmaSquare(unsigned long i, unsigned long j) {
    auto key = std::make_pair((i<j)?i:j,(j>i)?j:i); //key in preprocessed maps: (i,j) with i<j
    return computedMixingSigmaSquare.at(key);
}