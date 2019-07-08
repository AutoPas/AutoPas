#pragma once

#include <math.h>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"

using namespace std;
using namespace autopas;

class ParticleClassLibrary {
 public:
  ParticleClassLibrary(map<unsigned long, double> sigma, map<unsigned long, double> epsilon,
                       map<unsigned long, double> mass);

  ParticleClassLibrary();

  ParticleClassLibrary(const ParticleClassLibrary &pcl);

  ~ParticleClassLibrary() {}
  /**Getter for Particle Epsilon
   * @param Particle
   * @return Epsilon
   */
  double getEpsilon(Particle i);
  /**Getter for Particle Sigma
   * @param Particle
   * @return Sigma
   */
  double getSigma(Particle i);
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
  /**Returns the Epsilon of the MixingRule of 2 Particles
   * @param Particles; i and j
   * @return Epsilon of both
   * */
  double mixingE(unsigned long i, unsigned long j);
  /**Returns the Sigma of the MixingRule of 2 Particles
   * @param Particles; i and j
   * @return Sigma of both
   * */
  double mixingS(Particle i, Particle j);

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
