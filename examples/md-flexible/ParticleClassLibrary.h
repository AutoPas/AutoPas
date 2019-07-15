#pragma once

#include <math.h>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"

using namespace std;
using namespace autopas;

class ParticleClassLibrary {
 public:
  ParticleClassLibrary(map<unsigned long, double> &sigma, map<unsigned long, double> &epsilon,
                       map<unsigned long, double> &mass);

  ParticleClassLibrary();

  ParticleClassLibrary(const ParticleClassLibrary &pcl);

  ~ParticleClassLibrary() {}
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
  double getMass(Particle i);

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
  //@TODO static wäre vllt besser ???
  map<unsigned long, double> Epsilon;
  map<unsigned long, double> Sigma;
  map<unsigned long, double> Mass;
};