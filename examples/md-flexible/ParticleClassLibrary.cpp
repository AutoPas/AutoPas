#include "ParticleClassLibrary.h"

//@todo entferne unused functions
ParticleClassLibrary::ParticleClassLibrary(std::map<unsigned long, double> &sigma,
                                           std::map<unsigned long, double> &epsilon,
                                           std::map<unsigned long, double> &mass)
    : Epsilon(epsilon), Sigma(sigma), Mass(mass) {}

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
}

ParticleClassLibrary::ParticleClassLibrary() = default;

ParticleClassLibrary::ParticleClassLibrary(const ParticleClassLibrary &pcl) = default;

ParticleClassLibrary &ParticleClassLibrary::operator=(const ParticleClassLibrary &pcl) = default;

double ParticleClassLibrary::getMass(unsigned long i) {
return Mass.at(i);
}


double ParticleClassLibrary::get24Epsilon(unsigned long i) { return 24 * Epsilon.at(i); }

double ParticleClassLibrary::getSSigma(unsigned long i) {
  double sigma = Sigma.at(i);
  return sigma * sigma;
}

double ParticleClassLibrary::mixing24E(unsigned long i, unsigned long j) {
  return 24 * sqrt(Epsilon.at(i) * Epsilon.at(j));
}

double ParticleClassLibrary::mixingSS(unsigned long i, unsigned long j) {
  double mixingS = (Sigma.at(i) + Sigma.at(j)) / 2;
  return mixingS * mixingS;
}