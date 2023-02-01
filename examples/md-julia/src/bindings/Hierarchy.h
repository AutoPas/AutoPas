/**
 * TODO: think of a great comment
 */

namespace jlcxx
{
    using jlcxx::Parametric;
    
    /**
     *specifying BuildParameterList for TODO: good comment
     */
    
    /**
     * BuildParameterList for ParticleIteratorInterface
     */
    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::ParticleIteratorInterface<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };

    /**
     * BuildParameterList for ParticleIteratorWrapper
     */
    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::ParticleIteratorWrapper<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };

    /**
     * BuildParameterList for ParticlePropertiesLibrary
     */
    template<typename floatType, typename intType>
    struct BuildParameterList<ParticlePropertiesLibrary<floatType, intType>> {
        typedef ParameterList<floatType, intType> type;
    };

    /**
     * define inheritance relations between types
     */

    /**
     * define ParticleIteratorWrapper as derived type of ParticleIteratorInterface
     */
    template<typename Particle, bool modifiable> struct SuperType<autopas::ParticleIteratorWrapper<Particle, modifiable>> {typedef autopas::ParticleIteratorInterface<Particle, modifiable> type;};


    // adding templates for types
    /**
     * TODO: good comment for template parameters
     */
    template<typename floatType> struct IsMirroredType<autopas::MoleculeLJ<floatType>> : std::false_type { };
    template<typename floatType> struct IsMirroredType<MoleculeJ<floatType>> : std::false_type { };
    template<typename Particle> struct IsMirroredType<autopas::AutoPas<Particle>> : std::false_type { };
    template<typename T, typename V> struct IsMirroredType<autopas::ParticleBase<T, V>> : std::false_type { };
    template<typename Particle, bool modifiable> struct IsMirroredType<autopas::ParticleIteratorInterface<Particle, modifiable>> : std::false_type { };
    template<typename Particle, bool modifiable> struct IsMirroredType<autopas::ParticleIteratorWrapper<Particle, modifiable>> : std::false_type { };
    template<typename floatType, typename intType> struct IsMirroredType<ParticlePropertiesLibrary<floatType, intType>> : std::false_type { };
}