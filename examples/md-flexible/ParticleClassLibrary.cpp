#include "ParticleClassLibrary.h"

    /** Default constuktor wenn es nur ein Particle Type gibt
    **/
    template<typename floatType>
    ParticleClassLibrary<floatType>::ParticleClassLibrary(floatType &epsilon, floatType &sigma, floatType mass) {
  std::map<unsigned long, floatType> EMap;
  std::map<unsigned long, floatType> SMap;
  std::map<unsigned long, floatType> MMap;
    EMap.emplace(0, 24*epsilon);
    SMap.emplace(0, (sigma*sigma));
    MMap.emplace(0, mass);
  this->Epsilon24 = EMap;
  this->SigmaSquare = SMap;
  this->Mass = MMap;
}
//@todo need to initilized addType with ascending typeId
template<typename floatType>
void ParticleClassLibrary<floatType>::addType(unsigned long typeID, floatType epsilon, floatType sigma,
                                                   floatType mass) {
    for (auto e: Epsilon24) {
        unsigned long indexOfExistingEpsilon = std::get<0>(e);
        floatType secondEpsilon = std::get<1>(e);
        floatType epsilon24 = 24 * sqrt(epsilon * secondEpsilon);
        auto newEntry= std::make_pair(indexOfExistingEpsilon, typeID);
        computedMixing24Epsilon.emplace(newEntry,epsilon24);
    }
        Epsilon24.emplace(typeID, 24*epsilon);
    for (auto e: SigmaSquare) {
        unsigned long indexOfExistingSigma = std::get<0>(e);
        floatType existingSigma = std::get<1>(e);
        floatType newSigma = (sigma + existingSigma)/2;
        auto newEntry= std::make_pair(indexOfExistingSigma, typeID);
        computedMixing24Epsilon.emplace(newEntry, (newSigma*newSigma));
    }
        SigmaSquare.emplace(typeID, (sigma*sigma));
        Mass.emplace(typeID, mass);
}

template<typename floatType>
ParticleClassLibrary<floatType>::ParticleClassLibrary(const ParticleClassLibrary &pcl) = default;

template<typename floatType>
ParticleClassLibrary<floatType> &ParticleClassLibrary<floatType>::ParticleClassLibrary::operator=(const ParticleClassLibrary &pcl) = default;

template<typename floatType>
floatType ParticleClassLibrary<floatType>::getMass(unsigned long i) {
return Mass.at(i);
}
template<typename floatType>
floatType ParticleClassLibrary<floatType>::get24Epsilon(unsigned long i) { return Epsilon24.at(i); }

template<typename floatType>
floatType ParticleClassLibrary<floatType>::getSSigma(unsigned long i) {
  return SigmaSquare.at(i);
}

template<typename floatType>
floatType ParticleClassLibrary<floatType>::mixing24E(unsigned long i, unsigned long j) {
  auto key = std::make_pair((i<j)?i:j,(j>i)?j:i); //key in preprocessed maps: (i,j) with i<j
  return computedMixing24Epsilon.at(key);
}
template<typename floatType>
floatType ParticleClassLibrary<floatType>::mixingSS(unsigned long i, unsigned long j) {
    auto key = std::make_pair((i<j)?i:j,(j>i)?j:i); //key in preprocessed maps: (i,j) with i<j
    return computedMixingSigmaSquare.at(key);
}