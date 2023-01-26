#pragma once
#include "../MoleculeJ.h"

struct WrapMoleculeJ {
    template<typename floatType>
    void operator()(floatType&& wrapped) {
        // typedef typename T::type WrappedT;
        using WrappedT = typename floatType::type;

        // constructors for MoleculeLJ
        wrapped.template constructor<>();
        // wrapped.template constructor<jlcxx::ArrayRef<double, 1>, jlcxx::ArrayRef<double,1>, unsigned long, unsigned long>();

        // setters of MoleculeJ attributes
        wrapped.method("setPos", &WrappedT::setPos);
        wrapped.method("setV", &WrappedT::setV);
        wrapped.method("setF", &WrappedT::setF);
        wrapped.method("setOldF", &WrappedT::setOldF);

        // getters of MoleculeLJ attributes
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

        // further methods
        wrapped.method("toString", &WrappedT::toString);
    }
};