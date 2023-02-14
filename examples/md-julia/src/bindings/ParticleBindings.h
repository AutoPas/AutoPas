#pragma once
#include "../MoleculeJ.h"

/**
 * This class helps to wrap the autopas::ParticleBase type
 */
struct WrapParticleBase {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;
        wrapped.template constructor<>();
        // wrapped.method("print_particle", &WrapedT::print_particle);
    }
};

/**
 * This class helps to wrap the MoleculeJ type
 */
struct WrapMoleculeJ {
    template<typename floatType>
    void operator()(floatType&& wrapped) {
    
        using WrappedT = typename floatType::type;

        // default constructor for MoleculeJ
        wrapped.template constructor<>();

        // typename WrappedT::value_type*
        using ft = typename WrappedT::ft;
        
        // constructor for MoleculeJ with arguments for pos, v, moleculeId, typeId
        wrapped.template constructor<jlcxx::ArrayRef<ft, 1>, jlcxx::ArrayRef<ft,1>, unsigned long, unsigned long>();

        // setters of MoleculeJ attributes
        wrapped.method("setPosition", &WrappedT::setPosition);
        wrapped.method("setVelocity", &WrappedT::setVelocity);
        wrapped.method("setForce", &WrappedT::setForce);
        wrapped.method("setOldF", &WrappedT::setOldF);
        wrapped.method("setID", &WrappedT::setID);
        wrapped.method("setTypeId", &WrappedT::setTypeId);

        // getters of MoleculeJ attributes
        wrapped.method("getPosition", &WrappedT::getPosition);
        wrapped.method("getVelocity", &WrappedT::getVelocity);
        wrapped.method("getForce", &WrappedT::getForce);
        wrapped.method("getOldF", &WrappedT::getOldF);
        wrapped.method("getID", &WrappedT::getID);
        wrapped.method("getTypeId", &WrappedT::getTypeId);

        // add and sub methods of MoleculeJ attributes
        wrapped.method("addPosition", &WrappedT::addPosition);
        wrapped.method("addVelocity", &WrappedT::addVelocity);
        wrapped.method("addForce", &WrappedT::addForce);
        wrapped.method("subForce", &WrappedT::subForce);

        // get string representation of MoleculeJ object
        wrapped.method("toString", &WrappedT::toString);
    }
};

/**
 * This class helps to wrap the autopas::MoleculeLJ class.
 * TODO: if this class is necessary, add bindings which are needed into this wrapper class
 */
struct WrapMoleculeLJ {
    template<typename floatType>
    void operator()(floatType&& wrapped) {
        using WrappedT = typename floatType::type;

        wrapped.template constructor<>();
    }
};

/**
 * create Julia Module for Particle types
 */
JLCXX_MODULE define_module_particles(jlcxx::Module& mod) {
    using jlcxx::Parametric;
    using jlcxx::TypeVar;
    
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