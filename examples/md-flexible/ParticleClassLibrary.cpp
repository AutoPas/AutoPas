#include "ParticleClassLibrary.h"

    /** Default constuktor wenn es nur ein Particle Type gibt
    **/
    ParticleClassLibrary::ParticleClassLibrary(double &epsilon, double &sigma, double mass) {
  std::map<unsigned long, double> EMap;
  std::map<unsigned long, double> SMap;
  std::map<unsigned long, double> MMap;
    EMap.emplace(0, 24*epsilon);
    SMap.emplace(0, (sigma*sigma));
    MMap.emplace(0, mass);
  this->Epsilon24 = EMap;
  this->SigmaSquare = SMap;
  this->Mass = MMap;
}
//@todo need to initilized addType with ascending typeId

void ParticleClassLibrary::addType(unsigned long typeID, double epsilon, double sigma,
                                                   double mass) {
    for (auto e: Epsilon24) {
        unsigned long indexOfExistingEpsilon = std::get<0>(e);
        double secondEpsilon = std::get<1>(e);
        double epsilon24 = 24 * sqrt(epsilon * secondEpsilon);
        auto newEntry= std::make_pair(indexOfExistingEpsilon, typeID);
        computedMixing24Epsilon.emplace(newEntry,epsilon24);
    }
        Epsilon24.emplace(typeID, 24*epsilon);
    for (auto e: SigmaSquare) {
        unsigned long indexOfExistingSigma = std::get<0>(e);
        double existingSigma = std::get<1>(e);
        double newSigma = (sigma + existingSigma)/2;
        auto newEntry= std::make_pair(indexOfExistingSigma, typeID);
        computedMixing24Epsilon.emplace(newEntry, (newSigma*newSigma));
    }
        SigmaSquare.emplace(typeID, (sigma*sigma));
        Mass.emplace(typeID, mass);
}


ParticleClassLibrary::ParticleClassLibrary(const ParticleClassLibrary &pcl) = default;


ParticleClassLibrary &ParticleClassLibrary::operator=(const ParticleClassLibrary &pcl) = default;

double ParticleClassLibrary::getMass(unsigned long i) {
return Mass.at(i);
}
double ParticleClassLibrary::get24Epsilon(unsigned long i) { return Epsilon24.at(i); }

double ParticleClassLibrary::getSSigma(unsigned long i) {
  return SigmaSquare.at(i);
}

double ParticleClassLibrary::mixing24E(unsigned long i, unsigned long j) {
  auto key = std::make_pair((i<j)?i:j,(j>i)?j:i); //key in preprocessed maps: (i,j) with i<j
  return computedMixing24Epsilon.at(key);
}

double ParticleClassLibrary::mixingSS(unsigned long i, unsigned long j) {
    auto key = std::make_pair((i<j)?i:j,(j>i)?j:i); //key in preprocessed maps: (i,j) with i<j
    return computedMixingSigmaSquare.at(key);
}