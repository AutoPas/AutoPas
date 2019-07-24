#pragma once

#include <math.h>
#include <map>
#include <vector>
#include "autopas/particles/Particle.h"

class ParticleClassLibrary {
 public:
  ParticleClassLibrary(std::map<unsigned long, double> &sigma, std::map<unsigned long, double> &epsilon,
                       std::map<unsigned long, double> &mass);
  ParticleClassLibrary(double &epsilon, double &sigma, double mass, int numberOfParticles);
  ParticleClassLibrary();

  ParticleClassLibrary(const ParticleClassLibrary &pcl);

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
  template <class Particle>
  double getMass(const Particle &i) {
    return Mass.at(i.getID());
  }

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
  std::map<unsigned long, double> Epsilon;
  std::map<unsigned long, double> Sigma;
  std::map<unsigned long, double> Mass;
};