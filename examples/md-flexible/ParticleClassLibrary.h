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

private:
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
    return sqrt(Epsilon.at(i.getID())+Epsilon.at(j.getID()));
}

double ParticleClassLibrary::mixingS(Particle i, Particle j) {
    return (Sigma.at(i.getID())+Sigma.at(j.getID()))/2;
}