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
        
        // constructor for MoleculeJ with arguments for pos, v, moleculeId, typeId
        // wrapped.template constructor<jlcxx::ArrayRef<double, 1>, jlcxx::ArrayRef<double,1>, unsigned long, unsigned long>();

        // setters of MoleculeJ attributes
        wrapped.method("setPos", &WrappedT::setPos);
        wrapped.method("setV", &WrappedT::setV);
        wrapped.method("setF", &WrappedT::setF);
        wrapped.method("setOldF", &WrappedT::setOldF);
        wrapped.method("setID", &WrappedT::setID);

        // getters of MoleculeJ attributes
        wrapped.method("getPos", &WrappedT::getPos);
        wrapped.method("getV", &WrappedT::getV);
        wrapped.method("getF", &WrappedT::getF);
        wrapped.method("getOldF", &WrappedT::getOldF);
        wrapped.method("getID", &WrappedT::getID);

        // add and sub methods of MoleculeJ attributes
        wrapped.method("addPos", &WrappedT::addPos);
        wrapped.method("addV", &WrappedT::addV);
        wrapped.method("addF", &WrappedT::addF);
        wrapped.method("subF", &WrappedT::subF);

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
            .apply<autopas::ParticleBase<float, int>, autopas::ParticleBase<double, long>>(WrapParticleBase());

    /**
     * add MoleculeJ type to Julia with template parameters
     * double and float
     */
    mod.add_type<Parametric<TypeVar<1>>>("MoleculeJ")
            .apply<MoleculeJ<double>, MoleculeJ<float>>(WrapMoleculeJ());

    /**
     * add MoleculeLJ type to Julia with template parameter
     * double
     */
    mod.add_type<Parametric<TypeVar<1>>>("MoleculeLJ")
            .apply<autopas::MoleculeLJ<double>>(WrapMoleculeLJ());
}