#pragma once

#include <math.h>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"
template <typename floatType>
class ParticleClassLibrary {
 public:
  ParticleClassLibrary()=default;

  ParticleClassLibrary(floatType &epsilon, floatType &sigma, floatType mass);

  ParticleClassLibrary(const ParticleClassLibrary &pcl);

  ParticleClassLibrary &operator=(const ParticleClassLibrary &plc);

  void addType(unsigned long typeID,floatType epsilon, floatType sigma,floatType mass);

  ~ParticleClassLibrary() = default;
  /**Getter for Particle Epsilon*24
   * @param Particle
   * @return Epsilon*24
   */
  floatType get24Epsilon(unsigned long i);
  /**Getter for Particle Square Sigma
   * @param Particle
   * @return Sigma²
   */
  floatType getSSigma(unsigned long i);

  /**Getter for Particle Mass
   * @param Particle
   * @return Sigma
   */
  floatType getMass(unsigned long i);

  /**Returns (Epsilon*24) of the MixingRule of 2 Particles
   * @param Particles; i and j
   * @return 24*(Epsilon of both)
   * */
  floatType mixing24E(unsigned long i, unsigned long j);
  /**Returns Sigma Square of the MixingRule of 2 Particles
   * @param Particles; i and j
   * @return (Sigma of both)²
   * */
  floatType mixingSS(unsigned long i, unsigned long j);

 private:
  std::map<unsigned long, floatType> Epsilon24;
  std::map<unsigned long, floatType> SigmaSquare;
  std::map<unsigned long, floatType> Mass;
  std::map<std::pair<unsigned long,unsigned long>,floatType> computedMixing24Epsilon;
  std::map<std::pair<unsigned long,unsigned long>,floatType> computedMixingSigmaSquare;

};