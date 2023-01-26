#pragma once
#include "autopas/particles/Particle.h"
#include <jlcxx/jlcxx.hpp>

/**
 * Molecule class for md simulation in Julia
 */

template<typename floatType>
class MoleculeJ : public autopas::Particle {
    public:
    /**
     * default constructor
     */
    MoleculeJ() {}

    /**
     * Constructor 
     * @param position position of the molecule
     * @param velocity velocity of the molecule
     * @param moleculeId Id ot the molecule
     * @param typeId indicating the type of the molecule
     */
    MoleculeJ(jlcxx::ArrayRef<floatType,1> pos, jlcxx::ArrayRef<floatType,1> v,
        unsigned long moleculeId, unsigned long typeId = 0) : _typeId(typeId) {
        // : Particle(pos, v, moleculeId),
        std::array<floatType, 3> pos_{pos[0], pos[1], pos[2]};
        std::array<floatType, 3> v_{v[0], v[1], v[2]};
        autopas::Particle(pos_, v_, moleculeId);
    }

    /**
     * set position of the molecule
     */
    void setPos(jlcxx::ArrayRef<floatType,1> pos) {
        ParticleBase::setR({pos[0], pos[1], pos[2]});
    }

    /**
     * set velocity of the molecule
     */
    void setV(jlcxx::ArrayRef<floatType,1> v) {
        ParticleBase::setV({v[0], v[1], v[2]});
    }

    /**
     * set force of the molecule
     */
    void setF(jlcxx::ArrayRef<floatType,1> f) {
        ParticleBase::setF({f[0], f[1], f[2]});
    }

    /**
     * set oldForce of the molecule
     */
    void setOldF(jlcxx::ArrayRef<floatType,1> oldF) {
        _oldF = {oldF[0], oldF[1], oldF[2]};
    }

    /**
     * get position of the molecule
     */
    jlcxx::ArrayRef<double,1> getPos() {
        return {_r.data(), _r.size()};
    }

    /**
     * get velocity of the molecule
     */
    jlcxx::ArrayRef<double,1> getV() {
        return {_v.data(), _v.size()};
    }

    /**
     * get force of the molecule
     */
    jlcxx::ArrayRef<double,1> getF() {
        return {_f.data(), _f.size()};
    }

    /**
     * get oldForce of the molecule
     */
    jlcxx::ArrayRef<floatType,1> getOldF() {
        return {_oldF.data(), _oldF.size()};
    }

    /**
     * add position to current position of the molecule
     */
    void addPos(jlcxx::ArrayRef<floatType,1> pos) {
        ParticleBase::addR({pos[0], pos[1], pos[2]});
    }

    /**
     * add velocity to current velocity of the molecule
     */
    void addV(jlcxx::ArrayRef<floatType,1> v) {
        ParticleBase::addV({v[0], v[1], v[2]});
    }

    /**
     * add force to current force of the molecule
     */
    void addF(jlcxx::ArrayRef<floatType,1> f) {
        ParticleBase::addF({f[0], f[1], f[2]});
    }

    /**
     * substract force from current force of the molecule
     */
    void subF(jlcxx::ArrayRef<floatType,1> f) {
        ParticleBase::subF({f[0], f[1], f[2]});
    }

    private:
    /**
     * Particle type id.
     */
    size_t _typeId = 0;
    
    /**
     * Old Force of the particle experiences as 3D vector.
     */
    std::array<floatType,3> _oldF = {0.0, 0.0, 0.0};
};