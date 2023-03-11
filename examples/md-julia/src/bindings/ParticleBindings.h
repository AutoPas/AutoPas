#pragma once
#include "../MoleculeJ.h"
#include "autopas/particles/OwnershipState.h"

/**
 * This class helps to wrap the autopas::ParticleBase type
 */
struct WrapParticleBase {
    template<typename T>
    void operator()(T&& particle) {
        using ParticleType = typename T::type;
        
        particle.template constructor<>();

        // wrapped.method("print_particle", &WrapedT::print_particle);
    }
};

/**
 * This class helps to wrap the MoleculeJ type
 */
struct WrapMoleculeJ {
    template<typename T>
    void operator()(T&& particle) {
    
        using ParticleType = typename T::type;

        // default constructor for MoleculeJ
        particle.template constructor<>();

        // typename WrappedT::value_type*
        using ft = typename ParticleType::ft;
        
        // constructor for MoleculeJ with arguments for pos, v, moleculeId, typeId
        particle.template constructor<jlcxx::ArrayRef<ft, 1>, jlcxx::ArrayRef<ft,1>, unsigned long, unsigned long>();

        // setters of MoleculeJ attributes
        particle.method("setPosition", &ParticleType::setPosition);
        particle.method("setVelocity", &ParticleType::setVelocity);
        particle.method("setForce", &ParticleType::setForce);
        particle.method("setOldF", static_cast<void (ParticleType::*) (jlcxx::ArrayRef<double,1>)> (&ParticleType::setOldF));
        particle.method("setID", &ParticleType::setID);
        particle.method("setTypeId", &ParticleType::setTypeId);
        particle.method("setOwnershipState", &ParticleType::setOwnershipState);
        particle.method("setP", static_cast<void (ParticleType::*) (double, int)> (&ParticleType::setPo));
        particle.method("setV", static_cast<void (ParticleType::*) (double, int)> (&ParticleType::setVe));
        particle.method("setF", static_cast<void (ParticleType::*) (double, int)> (&ParticleType::setFo));
        particle.method("setOldF", static_cast<void (ParticleType::*) (double, int)> (&ParticleType::setOldFo));

        // getters of MoleculeJ attributes
        particle.method("getPosition", &ParticleType::getPosition);
        particle.method("getVelocity", &ParticleType::getVelocity);
        particle.method("getForce", &ParticleType::getForce);
        particle.method("getOldF",static_cast<jlcxx::ArrayRef<double,1> (ParticleType::*) ()> (&ParticleType::getOldF));
        particle.method("getID", &ParticleType::getID);
        particle.method("getTypeId", &ParticleType::getTypeId);
        particle.method("getOwnershipState", &ParticleType::getOwnershipState);
        particle.method("getP", static_cast<ft (ParticleType::*)(int)> (&ParticleType::getPo));
        particle.method("getV", static_cast<ft (ParticleType::*)(int)> (&ParticleType::getVe));
        particle.method("getF", static_cast<ft (ParticleType::*)(int)> (&ParticleType::getFo));
        particle.method("getOldF", static_cast<ft (ParticleType::*)(int)> (&ParticleType::getOldFo));

        // add and sub methods of MoleculeJ attributes
        particle.method("addPosition", &ParticleType::addPosition);
        particle.method("addVelocity", &ParticleType::addVelocity);
        particle.method("addForce", &ParticleType::addForce);
        particle.method("subForce", &ParticleType::subForce);

        // get string representation of MoleculeJ object
        particle.method("toString", &ParticleType::toString);
    }
};
//  static_cast<iterator_t (AutoPasType::*)(autopas::options::IteratorBehavior)> (&AutoPasType::begin));

/**
 * This class helps to wrap the autopas::MoleculeLJ class.
 * TODO: if this class is necessary, add bindings which are needed into this wrapper class
 */
struct WrapMoleculeLJ {
    template<typename T>
    void operator()(T&& particle) {
        using ParticleType = typename T::type;

        particle.template constructor<>();
    }
};

/**
 * create Julia Module for Particle types
 */
JLCXX_MODULE define_module_particles(jlcxx::Module& mod) {
    using jlcxx::Parametric;
    using jlcxx::TypeVar;

    /**
     * add OwnershipState values to Julia
     */
    mod.add_bits<autopas::OwnershipState>("OwnershipState", jlcxx::julia_type("CppEnum"));
    mod.set_const("dummyS", autopas::OwnershipState::dummy);    
    mod.set_const("ownedS", autopas::OwnershipState::owned);
    mod.set_const("haloS", autopas::OwnershipState::halo);
    
    /**
     * add ParticleBase type to Julia with template parameters
     * float, int and double, long
     */
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleBase")
            .apply<autopas::ParticleBase<float, int>, autopas::ParticleBase<double, unsigned long>>(WrapParticleBase());

    /**
     * add MoleculeJ type to Julia with template parameters
     * double and float
     */
    mod.add_type<Parametric<TypeVar<1>>>("MoleculeJ")
            .apply<MoleculeJ<double>>(WrapMoleculeJ());

    /**
     * add MoleculeLJ type to Julia with template parameter
     * double
     */
    mod.add_type<Parametric<TypeVar<1>>>("MoleculeLJ")
            .apply<autopas::MoleculeLJ<double>>(WrapMoleculeLJ());
}