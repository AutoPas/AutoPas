#pragma once
#include <jlcxx/jlcxx.hpp>
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"

/**
 * This class is a helper class to wrap the ParticlePropertiesLibrary with CxxWrap
 */
struct WrapParticlePropertiesLibrary{
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;

        // constructor of ParticlePropertyLibrary to Julia
        wrapped.template constructor<double>();
        
        // add another type of particle to particleproperties
        wrapped.method("addType", &WrappedT::addType);

        // calculate 
        wrapped.method("calculateMixingCoefficients", &WrappedT::calculateMixingCoefficients);
        wrapped.method("calcShift6", &WrappedT::calcShift6);        
        
        // get 24epsilon of one type 
        wrapped.method("get24Epsilon", &WrappedT::get24Epsilon);
        // get sigma * sigma of one type
        wrapped.method("getSigmaSquare", &WrappedT::getSigmaSquare);
        // get mass of one type
        wrapped.method("getMass", &WrappedT::getMass);

    	// get mixed 24epsilon of two types i and j
        wrapped.method("mixing24Epsilon", &WrappedT::mixing24Epsilon);
        // get mixed sigma*sigma of two types i and j
        wrapped.method("mixingSigmaSquare", &WrappedT::mixingSigmaSquare);
        // get mixed shift6 of two types i and j
        wrapped.method("mixingShift6", &WrappedT::mixingShift6);

        // TDOO: check if necessary
        // wrapped.method("getMixingData", &WrappedT::getMixingData); // return type?
        // wrapped.method("getMixingDataPtr", &WrappedT::getMixingDataPtr); // return type?
    }
};

/**
 * create Julia Module with ParticlePropertyLibrary type
 */
JLCXX_MODULE define_module_properties(jlcxx::Module& mod) {
    using jlcxx::Parametric;
    using jlcxx::TypeVar;

    /**
     * add ParticlePropertiesLibrary type to Julia (with template parameters floatType = double and intType = int)
     */
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticlePropertiesLibrary")
            .apply<ParticlePropertiesLibrary<double, int>>(WrapParticlePropertiesLibrary());
}