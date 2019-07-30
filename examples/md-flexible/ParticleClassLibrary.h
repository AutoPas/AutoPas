#pragma once

#include <math.h>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"
//@todo add template parameter -> did resolve into undefined references
class ParticleClassLibrary {
 public:
  ParticleClassLibrary()=default;

  ParticleClassLibrary(double &epsilon, double &sigma, double mass);

  ParticleClassLibrary(const ParticleClassLibrary &pcl);

  ParticleClassLibrary &operator=(const ParticleClassLibrary &plc);

  void addType(unsigned long typeID,double epsilon, double sigma,double mass);

  ~ParticleClassLibrary() = default;
  /**Getter for Particle Epsilon*24
   * @param Particle
   * @return Epsilon*24
   */
  double get24Epsilon(unsigned long i);
  /**Getter for Particle Square Sigma
   * @param Particle
   * @return Sigma²
   */
  double getSSigma(unsigned long i);

  /**Getter for Particle Mass
   * @param Particle
   * @return Sigma
   */
  double getMass(unsigned long i);

  /**Returns (Epsilon*24) of the MixingRule of 2 Particles
   * @param Particles; i and j
   * @return 24*(Epsilon of both)
   * */
  double mixing24E(unsigned long i, unsigned long j);
  /**Returns Sigma Square of the MixingRule of 2 Particles
   * @param Particles; i and j
   * @return (Sigma of both)²
   * */
  double mixingSS(unsigned long i, unsigned long j);

 private:
  std::map<unsigned long, double> Epsilon24;
  std::map<unsigned long, double> SigmaSquare;
  std::map<unsigned long, double> Mass;
  std::map<std::pair<unsigned long,unsigned long>,double> computedMixing24Epsilon;
  std::map<std::pair<unsigned long,unsigned long>,double> computedMixingSigmaSquare;

};