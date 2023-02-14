/**
 * This class helps to wrap the ParticleIteratorInterface type
 */
struct WrapIteratorInterface {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;        
    }
};

/**
 * This class helps to wrap the ParticleInteratorWrapper type
 */
struct WrapIteratorWrapper {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;
        // check if iterator is still valid/next element exists
        wrapped.method("isValid", &WrappedT::isValid);
        // increase iterator / ++ operator
        wrapped.method("++", &WrappedT::operator++);
        // dereference iterator / * operator
        wrapped.method("*", &WrappedT::operator*);
    }
};

/**
 * define Julia module wtih Iterators for AutoPas container objects
 */

JLCXX_MODULE define_module_iterators(jlcxx::Module& mod) {
    using jlcxx::Parametric;
    using jlcxx::TypeVar;

    /**
     * add ParticleIteratorInterface type to Julia with template parameters:
     * 1) autopas::MoleculeLJ<double>, true
     * 2) MoleculeJ<double>, true
     */
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorInterface")
            .apply<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>, autopas::ParticleIteratorInterface<MoleculeJ<double>, true>>(WrapIteratorInterface());

    /**
     * add ParticleIteratorWrapper type to Julia with template parameters:
     * 1) autopas::MoleculeLJ<double>, true
     * 2) MoleculeJ<double>, true 
     */
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorWrapper", jlcxx::julia_base_type<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>())
            .apply<autopas::ParticleIteratorWrapper<autopas::MoleculeLJ<double>, true>, autopas::ParticleIteratorWrapper<MoleculeJ<double>, true>>(WrapIteratorWrapper());
    
}