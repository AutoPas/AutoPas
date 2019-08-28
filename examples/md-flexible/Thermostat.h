/**
 * @file Thermostat.h
 * @author N. Fottner
 * @date 27/8/19
 */

#pragma once
#include "autopas/AutoPas.h"
#include <cstdlib>
#include "autopas/utils/ArrayMath.h"

/**
 * Thermostat to adjust the Temperature of the Simulation
 * WIP
 */
template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
class Thermostat {
public:
    using ThermostatFloatType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryFloatType;
    using ThermostatIntType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryIntType;

    /**
     * Constructor if target Temperature is specified
     * @param t_init initialTemperature of System
     * @param initBM initialization with BM or other Formula
     * @param t_target target Temperature of System
     * @param delta_temp
     * @param ParticlePropertiesLibrary
     * */
     //@todo FRAGE : always use Brownian Motion??
    Thermostat(ThermostatFloatType t_init, bool initBM, ThermostatFloatType t_target, ThermostatFloatType delta_temp, ParticlePropertiesLibraryTemplate &ParticlePropertiesLibrary);

    /**
     * Constructor if only initial Temperature is specified
     * @param t_init
     * @param initBM
     * @param ParticlePropertiesLibrary
     * */
    Thermostat(ThermostatFloatType t_init, bool initBM,ParticlePropertiesLibraryTemplate &ParticlePropertiesLibrary);

    /**
     * Default Destructor
     * */
    ~Thermostat() =default;

    /**
     * Copy Constructor
     * @param ThermostatCopy
     * */
    Thermostat(const Thermostat &ThermostatCopy) = default;

    /**
     * Copy Assignment Constructor
     * @param thermo
     * @return
     * */
    Thermostat &operator=(const Thermostat &thermo) = default;

    /**
     * Initializes the Simulation according to brownian movement.
     * @param ps Particle system to be initialized.
     */
    void initialize(AutoPasTemplate &autopas);

    /**
     * Scales velocity of particles to reach desired temperature.
     * @param ps Particle system
     */
    void applyThermo(AutoPasTemplate &autopas);

    /**
     * Calculates temperature of system.
     * @param c Particle System.
     * @return Temperature of system.
     */
    ThermostatFloatType temperature(AutoPasTemplate &autopas);

    /**
    * Add a random velocity according to the Maxwell-Boltzmann distribution to the
    * particles, with a given mean velocity.
    * code taken from the MolSim Course
    *
    * @param p The particle to initialize.
    * @param factor
     */
    void MaxwellBoltzmannDistribution(autopas::Particle &p,const ThermostatFloatType factor) {
        p.setV(autopas::ArrayMath::addScalar(p.getV(),factor*GaussDeviate()));
    }

    /**
     * helper function for MaxwellBoltzmannDistribution().
     * Generates a gauss deviate, i.e. values according to the normal distribution.
     * code taken from:
     * Griebel et. al.: Numerical Simulation in Molecular Dynamics, p. 427
     */
    static ThermostatFloatType GaussDeviate() {
        ThermostatFloatType a1, a2, s, r, b1;
        static bool iset = false;
        static ThermostatFloatType b2;

        if (!iset) {
            do {
                a1 = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
                a2 = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
                r = a1 * a1 + a2 * a2;
            } while (r >= 1.0);
            s = sqrt(-2.0 * log(r) / r);
            b1 = a1 * s;
            b2 = a2 * s;
            iset = true;
            return b1;
        } else {
            iset = false;
            return b2;
        }
    }

private:

    /**
     * Initial temperature.
     */
    ThermostatFloatType t_init;

    /**
     * Target temperature.
     */
    ThermostatFloatType t_target;

    /**
     * specifies, how the velocity values will be initialized, with Brownian Motion or the new Formula
     *  initBM == true -> we use 0.1 as the factor in the MaxwellBoltzmannDistribution calculation
     * */
    bool initBM;

    /**
     * DELTA T
     * Temperature difference per thermostat application until t_target is reached.
     */
    ThermostatFloatType delta_temp;

    /**
     * ParticlePropertiesLibrary to access Mass of Particles
     * */
     ParticlePropertiesLibraryTemplate _particlePropertiesLibrary;

};
template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void Thermostat<AutoPasTemplate,ParticlePropertiesLibraryTemplate>::initialize(AutoPasTemplate &autopas) {
    ThermostatFloatType kb= 1.38064852 * pow(10, -23); //Boltzmann constant in J/K
    if (initBM) {
        //we consider: "However, the initialization with the Brownian Motion should be optional." means factor=0.1
        #ifdef AUTOPAS_OPENMP
        #pragma omp parallel
        #endif
        for(auto iter=autopas.begin();iter.isValid();++iter) {
            MaxwellBoltzmannDistribution(*iter, 0.1);
        }
    } else {
        ThermostatFloatType currentTempMulKB = temperature(autopas) * kb;
        #ifdef AUTOPAS_OPENMP
        #pragma omp parallel
        #endif
        for(auto iter=autopas.begin();iter.isValid();++iter) {
           ThermostatFloatType factor = pow((currentTempMulKB/_particlePropertiesLibrary.getMass(iter->getTypeId())), 0.5);
            MaxwellBoltzmannDistribution(*iter, factor);
        }
    }
}

template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void Thermostat<AutoPasTemplate,ParticlePropertiesLibraryTemplate>::applyThermo(AutoPasTemplate &autopas) {
    ThermostatFloatType temp = temperature(autopas);
    //wir nehmen an, dass t_target=t_init, da die temperature ja konstant bleiben soll
    static ThermostatFloatType scaling;
    if (t_init == t_target) {
        scaling = pow((this->t_init / temp), 0.5);
    } else {
        t_init = t_init + delta_temp;
        scaling = pow((t_init / temp), 0.5);
    }
    #ifdef AUTOPAS_OPENMP
    #pragma omp parallel
    #endif
    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
        iter->setV(autopas::ArrayMath::mulScalar(iter->getV(), scaling));
    }
}

template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryFloatType Thermostat<AutoPasTemplate,ParticlePropertiesLibraryTemplate>::temperature(AutoPasTemplate &autopas) {
    ThermostatFloatType kinetic = 0;
    //@todo anpassen mit OpenMp
    for(auto iter=autopas.begin();iter.isValid();++iter) {
            auto vel =iter->getV();
            kinetic += _particlePropertiesLibrary.getMass(iter->getTypeId())* /*scalarProduct of Particles velocity:*/ (vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]) / 2;
    }
    ThermostatFloatType temp = (2 * kinetic) / (autopas.getNumberOfParticles() * 3); //AutoPas works always on 3 dimensions
    return temp;
}

template< class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
Thermostat<AutoPasTemplate,ParticlePropertiesLibraryTemplate>::Thermostat(ThermostatFloatType t_init,bool initBM,ThermostatFloatType t_target,ThermostatFloatType delta_temp,ParticlePropertiesLibraryTemplate &ParticlePropertiesLibrary) :
                                                                                                  t_init(t_init),
                                                                                                  t_target(t_target),
                                                                                                  initBM(initBM),
                                                                                                  delta_temp(delta_temp),
                                                                                                  _particlePropertiesLibrary(ParticlePropertiesLibrary){}
template<class AutoPasTemplate,class ParticlePropertiesLibraryTemplate>
Thermostat<AutoPasTemplate,ParticlePropertiesLibraryTemplate>::Thermostat(ThermostatFloatType t_init,bool initBM,ParticlePropertiesLibraryTemplate &ParticlePropertiesLibrary) :  t_init(t_init), initBM(initBM), _particlePropertiesLibrary(ParticlePropertiesLibrary) {
    t_target = t_init;
    delta_temp = 1;
}









