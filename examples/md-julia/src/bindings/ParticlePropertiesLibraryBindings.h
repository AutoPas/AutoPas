#pragma once
#include <jlcxx/jlcxx.hpp>
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"

/**
 * This class is a helper class to wrap the ParticlePropertiesLibrary with CxxWrap
 */
struct WrapParticlePropertiesLibrary{
    template<typename T>
    void operator()(T&& pplibrary) {
        using PplType = typename T::type;

        // constructor of ParticlePropertyLibrary to Julia
        pplibrary.template constructor<double>();
        
        // add another type of particle to particleproperties
        pplibrary.method("addType", &PplType::addType);

        // calculate 
        pplibrary.method("calculateMixingCoefficients", &PplType::calculateMixingCoefficients);
        pplibrary.method("calcShift6", &PplType::calcShift6);        
        
        // get 24epsilon of one type 
        pplibrary.method("get24Epsilon", &PplType::get24Epsilon);
        // get sigma * sigma of one type
        pplibrary.method("getSigmaSquare", &PplType::getSigmaSquare);
        // get mass of one type
        pplibrary.method("getMass", &PplType::getMass);

    	// get mixed 24epsilon of two types i and j
        pplibrary.method("mixing24Epsilon", &PplType::mixing24Epsilon);
        // get mixed sigma*sigma of two types i and j
        pplibrary.method("mixingSigmaSquare", &PplType::mixingSigmaSquare);
        // get mixed shift6 of two types i and j
        pplibrary.method("mixingShift6", &PplType::mixingShift6);

        // TDOO: check if necessary
        // pplibrary.method("getMixingData", &PplType::getMixingData); // return type?
        // pplibrary.method("getMixingDataPtr", &PplType::getMixingDataPtr); // return type?
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
            .apply<ParticlePropertiesLibrary<double, size_t>>(WrapParticlePropertiesLibrary());
}