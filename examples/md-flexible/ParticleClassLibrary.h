#include "autopas/AutoPas.h"
#include <map>
#include <vector>
#include <math.h>

using namespace std;
using namespace autopas;

class ParticleClassLibrary {
public:
    ParticleClassLibrary(map<unsigned long,double> sigma,map<unsigned long, double> epsilon, map<unsigned long,double> mass);

    ~ParticleClassLibrary() {}

    double getEpsilon(Particle i);
    double getSigma(Particle i);
    double getMass(Particle i);

    double mixingE(Particle i, Particle j);
    double mixingS(Particle i, Particle j);

    map<unsigned long, double> Epsilon;
    map<unsigned long, double> Sigma;
    map<unsigned long, double> Mass;


};


ParticleClassLibrary::ParticleClassLibrary(map<unsigned long,double> sigma,map<unsigned long, double> epsilon, map<unsigned long,double> mass):Epsilon(epsilon),Sigma(sigma),Mass(mass) {}

double ParticleClassLibrary::getSigma(Particle i) {
    return Epsilon.at(i.getID());
}

double ParticleClassLibrary::getEpsilon(Particle i) {
    return Sigma.at(i.getID());
}

double ParticleClassLibrary::getMass(Particle i) {
    return Mass.at(i.getID());
}

double ParticleClassLibrary::mixingE(Particle i, Particle j) {
    double iEpsi = Epsilon.at(i.getID());
    double jEpsi = Epsilon.at(j.getID());
    return sqrt(iEpsi+jEpsi);
}

double ParticleClassLibrary::mixingS(Particle i, Particle j) {
    double iSig = Sigma.at(i.getID());
    double jSig = Sigma.at(j.getID());
    return (iSig+jSig)/2;
}