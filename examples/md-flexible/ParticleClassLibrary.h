
#include "autopas/particles/Particle.h"
#include <map>
#include <vector>
#include <math.h>

using namespace std;
using namespace autopas;

class ParticleClassLibrary {
public:
    ParticleClassLibrary(map<unsigned long,double> sigma,map<unsigned long, double> epsilon, map<unsigned long,double> mass);

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
    double get24Epsilon(Particle i);
    /**Getter for Particle Square Sigma
    * @param Particle
    * @return Sigma²
    */
    double getSSigma(Particle i);

    /**Getter for Particle Mass
    * @param Particle
    * @return Sigma
    */
    double getMass(Particle i);
    /**Returns the Epsilon of the MixingRule of 2 Particles
     * @param Particles; i and j
     * @return Epsilon of both
     * */
    double mixingE(Particle i, Particle j);
    /**Returns the Sigma of the MixingRule of 2 Particles
    * @param Particles; i and j
    * @return Sigma of both
    * */
    double mixingS(Particle i, Particle j);

    /**Returns (Epsilon*24) of the MixingRule of 2 Particles
     * @param Particles; i and j
     * @return 24*(Epsilon of both)
     * */
    double mixing24E(Particle i, Particle j);
    /**Returns Sigma Square of the MixingRule of 2 Particles
    * @param Particles; i and j
    * @return (Sigma of both)²
    * */
    double mixingSS(Particle i, Particle j);



private:
    //@TODO static wäre vllt besser ???
    map<unsigned long, double> Epsilon;
    map<unsigned long, double> Sigma;
    map<unsigned long, double> Mass;
};

ParticleClassLibrary::ParticleClassLibrary(map<unsigned long,double> sigma,map<unsigned long, double> epsilon, map<unsigned long,double> mass):Epsilon(epsilon),Sigma(sigma),Mass(mass) {}


double ParticleClassLibrary::getEpsilon(Particle i) {
    return Epsilon.at(i.getID());
}

 double ParticleClassLibrary::getSigma(Particle i) {
    return Sigma.at(i.getID());
}

double ParticleClassLibrary::get24Epsilon(Particle i) {
    return 24*Epsilon.at(i.getID());
}

double ParticleClassLibrary::getSSigma(Particle i) {
    double sigma=Sigma.at(i.getID());
    return sigma*sigma;
}


double ParticleClassLibrary::getMass(Particle i) {
    return Mass.at(i.getID());
}

double ParticleClassLibrary::mixingE(Particle i, Particle j) {
    return sqrt(Epsilon.at(i.getID())+Epsilon.at(j.getID()));
}

double ParticleClassLibrary::mixingS(Particle i, Particle j) {
    return (Sigma.at(i.getID())+Sigma.at(j.getID()))/2;
}

double ParticleClassLibrary::mixing24E(Particle i, Particle j) {
    return 24 * sqrt(Epsilon.at(i.getID())+Epsilon.at(j.getID()));
}

double ParticleClassLibrary::mixingSS(Particle i, Particle j) {
    double mixingS=(Sigma.at(i.getID())+Sigma.at(j.getID()))/2;
    return mixingS*mixingS;
}