#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <type_traits>
#include "jlcxx/functions.hpp"
#include "autopas/molecularDynamics/MoleculeLJ.h"
// #include "autopas/AutoPas.h"
#include "TypeDef.h"
// #include "Simulation.h"
#include "autopas/AutoPas.h"
#include "Hydrogen.h"
#include <iostream>

// extern template class autopas::AutoPas<ParticleType>;
// template class autopas::AutoPas<ParticleType>;


struct WrapMoleculeLJ {
    template<typename floatType>
    void operator()(floatType&& wrapped) {
        // typedef typename T::type WrappedT;
        using WrappedT = typename floatType::type;
        wrapped.template constructor<>();
        wrapped.template constructor<jlcxx::ArrayRef<double,1>, jlcxx::ArrayRef<double,1>, int, int>();
        wrapped.method("setPos", &WrappedT::setPos);
        wrapped.method("setV", &WrappedT::setV);
        wrapped.method("toString", &WrappedT::toString);
        // wrapped.method("print_particle", &WrappedT::print_particle);
    }
};

struct WrapParticleBase {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;
        wrapped.template constructor<>();
        // wrapped.method("print_particle", &WrapedT::print_particle);
    }
};

struct WrapIteratorInterface {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrapedT = typename T::type;
        // wrapped.template constructor<>();
        
    }
};

/*
struct WrapIteratorInterfaceImpl {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;
        // wrapped.template constructor<>();
    }

};
*/

struct WrapIteratorWrapper {
    template<typename T>
    void operator()(T&& wrapped) {
        using WrappedT = typename T::type;
        // wrapped.template constructor<>();
        // wrapped.method("setIterator", &WrappedT::setInterator);
        wrapped.method("isValid", &WrappedT::isValid);
        wrapped.method("inc", &WrappedT::operator++);
        wrapped.method("deref", &WrappedT::operator*);
        // TODO: wrap operator
    }
};

namespace jlcxx
{
    using jlcxx::Parametric;
    // specifying BuildParameterList
    
    
    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::ParticleIteratorInterface<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };

    /*
    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::internal::ParticleIteratorInterfaceImpl<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };
    */

    template<typename Particle, bool modifiable>
    struct BuildParameterList<autopas::ParticleIteratorWrapper<Particle, modifiable>> {
        typedef ParameterList<Particle, bool> type;
    };
    
    // adding inheritance for types
    // template<typename Particle, bool modifiable> struct SuperType<autopas::internal::ParticleIteratorInterfaceImpl<Particle, modifiable>> {typedef autopas::ParticleIteratorInterface<Particle, modifiable> type;};
    template<typename Particle, bool modifiable> struct SuperType<autopas::ParticleIteratorWrapper<Particle, modifiable>> {typedef autopas::ParticleIteratorInterface<Particle, modifiable> type;};
    // template<> struct SuperType<autopas::options::IteratorBehavior> {typedef autopas::options::Option<autopas::options::IteratorBehavior> type;};

    // adding templates for types
    template<typename floatType> struct IsMirroredType<autopas::MoleculeLJ<floatType>> : std::false_type { };
    template<typename Particle> struct IsMirroredType<autopas::AutoPas<Particle>> : std::false_type { };
    template<typename T, typename V> struct IsMirroredType<autopas::ParticleBase<T, V>> : std::false_type { };
    template<typename Particle, bool modifiable> struct IsMirroredType<autopas::ParticleIteratorInterface<Particle, modifiable>> : std::false_type { };
    // template<typename Particle, bool modifiable> struct IsMirroredType<autopas::internal::ParticleIteratorInterfaceImpl<Particle, modifiable>> : std::false_type{ };
    template<typename Particle, bool modifiable> struct IsMirroredType<autopas::ParticleIteratorWrapper<Particle, modifiable>> : std::false_type { };
    // template<typename actualOption> struct IsMirroredType<autopas::options::Option<actualOption>> : std::false_type { };

}

struct WrapAutoPas {
    // using iterator_t = typename autopas::IteratorTraits<autopas::MoleculeLJ<double>>::iterator_t;
    template<typename Particle>
    void operator()(Particle&& wrapped) {
        // typedef typename T::type WrappedT;
        using WrappedT = typename Particle::type;
        wrapped.template constructor<>();
        wrapped.method("init", &WrappedT::init);
        // resizeBox
        // force retune
        wrapped.method("finalize", &WrappedT::finalize);
        wrapped.method("updateContainer", &WrappedT::updateContainer);
        wrapped.method("addParticle", &WrappedT::addParticle);
        // addhaloParticle
        wrapped.method("deleteAllParticles", &WrappedT::deleteAllParticles);
        // .method("open", static_cast<void (IO::LCReader::*)(const std::string&)>(&IO::LCReader::open));
        //  t0.method("get_x_pos", static_cast<int (WrapClass::*)() >(&WrapClass::get_x_pos));
        //   t0.method("set_x_pos", static_cast<void (WrapClass::*)(int) >(&WrapClass::set_x_pos));
        // wrapped.method("begin", static_cast<iterator_t (WrappedT::*)(autopas::IteratorBehavior) > (&WrappedT::begin));
        // wrapped.method("begin", static_cast<autopas::ParticleIteratorWrapper<Particle, bool> (WrappedT::*)() > (&WrappedT::begin));
        wrapped.method("begin", &WrappedT::begin);
        // wrapped.method("begin", static_cast<autopas::ParticleIteratorWrapper<autopas::MoleculeLJ<double>, true> (WrappedT::*)() >(&WrappedT::begin));
        // wrapped.method("deleteParticle", &Wrapped::deleteParticle); may not be wrapped
        wrapped.method("printBoxSize", &WrappedT::printBoxSize);
    }
};

/*
struct WrapOptions {
    template<typename actualOption>
    void operator()(actualOption&& wrapped) {
        using WrappedT = typename actualOption::type;
    }
};
*/


std::vector<double> get_vec(double x1, double x2, double x3) {
    return {x1, x2, x3};
}


void setBox(autopas::AutoPas<autopas::MoleculeLJ<double>>& ap) {
    const std::array<double, 3> min_{0,0,0};
    const std::array<double, 3> max_{7.5,7.5,7.5};
    
    ap.setBoxMin(min_);
    ap.setBoxMax(max_);
    /*
    auto m1 = ap.getBoxMin();
    auto m2 = ap.getBoxMax();
    std::cout << "min: " << m1.at(0) << ", " << m1.at(1) << ", " << m1.at(2) << "\n";
    std::cout << "min: " << m2.at(0) << ", " << m2.at(1) << ", " << m2.at(2) << "\n";
    */
}

void getBox(autopas::AutoPas<autopas::MoleculeLJ<double>>& ap) {
    auto m1 = ap.getBoxMin();
    auto m2 = ap.getBoxMax();
    std::cout << "get_min: " << m1.at(0) << ", " << m1.at(1) << ", " << m1.at(2) << "\n";
    std::cout << "get_min: " << m2.at(0) << ", " << m2.at(1) << ", " << m2.at(2) << "\n";
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    using jlcxx::Parametric;
    using jlcxx::TypeVar;
    
    // add type MoleculeLJ
    mod.add_type<Parametric<TypeVar<1>>>("MoleculeLJ")
            .apply<autopas::MoleculeLJ<float>, autopas::MoleculeLJ<double>>(WrapMoleculeLJ());

    mod.add_type<Parametric<TypeVar<1>>>("ParticleBase")
            .apply<autopas::ParticleBase<float, int>, autopas::ParticleBase<double, long>>(WrapParticleBase());

    mod.add_type<Hydrogen>("Hydrogen")
            .constructor<>()
            .constructor<int, std::vector<double>, std::vector<double>>()
            .method("update_x", &Hydrogen::update_x)
            .method("update_v", &Hydrogen::update_v);

    /*
    mod.add_type<autopas::options::IteratorBehavior>("IteratorBehavior", jlcxx::julia_base_type<autopas::options::Option<autopas::options::IteratorBehavior>>())
            .constructor<>();

    mod.add_type<Parametric<TypeVar<1>>>("Option")
            .apply<autopas::options::Option<autopas::options::IteratorBehavior>>(WrapOptions());

    */
    // add iterator types necessary to iterate over the AutoPasContainer
    // add type ParticleIteratorInterface
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorInterface")
            .apply<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>(WrapIteratorInterface());

    
    /*
    // add type ParticleIteratorInterfaceImpl
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorInterfaceImpl", jlcxx::julia_base_type<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>())
            .apply<autopas::internal::ParticleIteratorInterfaceImpl<autopas::MoleculeLJ<double>, true>>(WrapIteratorInterfaceImpl());
    */

    // add type ParticleIteratorWrapper which is used to iterate over the AutoPasContainer
    mod.add_type<Parametric<TypeVar<1>, TypeVar<2>>>("ParticleIteratorWrapper", jlcxx::julia_base_type<autopas::ParticleIteratorInterface<autopas::MoleculeLJ<double>, true>>())
            .apply<autopas::ParticleIteratorWrapper<autopas::MoleculeLJ<double>, true>>(WrapIteratorWrapper());
    
    // add type AutoPas
    // mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
    //         .apply<autopas::AutoPas<autopas::ParticleBase<float, int>>, autopas::AutoPas<autopas::MoleculeLJ<double>>>(WrapAutoPas());
    
    mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
            .apply<autopas::AutoPas<autopas::MoleculeLJ<double>>>(WrapAutoPas());

    // mod.add_type<Parametric<TypeVar<1>>>("AutoPas")
    //         .apply<autopas::AutoPas<Hydrogen>, autopas::AutoPas<ParticleType>>(WrapAutoPas());

    mod.method("get_vec", &get_vec);
    mod.method("setBox", &setBox);
    mod.method("getBox", &getBox);
}

// cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DAUTOPAS_BUILD_TESTS ..