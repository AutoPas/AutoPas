#pragma once
#include "autopas/particles/ParticleBase.h"
#include "autopas/molecularDynamics/MoleculeLJ.h"
#include <jlcxx/jlcxx.hpp>

/**
 * Molecule class for md simulation in Julia
 */

template<typename floatType>
struct MoleculeJ : public autopas::MoleculeLJ<floatType> {
    public:

    typedef floatType ft;

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
        unsigned long moleculeId, unsigned long typeId = 0) {
        // : Particle(pos, v, moleculeId),
        std::array<floatType, 3> pos_{pos[0], pos[1], pos[2]};
        std::array<floatType, 3> v_{v[0], v[1], v[2]};
        setPosition(pos);
        setVelocity(v);
        setID(moleculeId);
        setTypeId(typeId);

        // autopas::Particle(pos_, v_, moleculeId);
    }
    
    MoleculeJ(std::array<floatType, 3> pos, std::array<floatType, 3> v, unsigned long moleculeId,
                      unsigned long typeId = 0)
      : autopas::Particle(pos, v, moleculeId) {}

    /**
     * set position of the molecule
     */
    void setPosition(jlcxx::ArrayRef<floatType,1> pos) {
        autopas::MoleculeLJ::setR({pos[0], pos[1], pos[2]});
    }

    /**
     * set velocity of the molecule
     */
    void setVelocity(jlcxx::ArrayRef<floatType,1> v) {
        autopas::MoleculeLJ::setV({v[0], v[1], v[2]});
    }

    /**
     * set force of the molecule
     */
    void setForce(jlcxx::ArrayRef<floatType,1> f) {
        autopas::MoleculeLJ::setF({f[0], f[1], f[2]});
    }

    /**
     * set oldForce of the molecule
     */
    void setOldForce(jlcxx::ArrayRef<floatType,1> oldF) {
        // _oldF = {oldF[0], oldF[1], oldF[2]};
    }

    /**
     * get position of the molecule
     */
    jlcxx::ArrayRef<double,1> getPosition() {
        return {autopas::MoleculeLJ::_r.data(), autopas::MoleculeLJ::_r.size()};
    }

    /**
     * get velocity of the molecule
     */
    jlcxx::ArrayRef<double,1> getVelocity() {
        return {_v.data(), _v.size()};
    }

    /**
     * get force of the molecule
     */
    jlcxx::ArrayRef<double,1> getForce() {
        return {_f.data(), _f.size()};
    }

    /**
     * get oldForce of the molecule
     */
    jlcxx::ArrayRef<floatType,1> getOldForce() {
        return {_oldF.data(), _oldF.size()};
    }

    // size_t getTypeId() const { return _typeId;}

    /**
     * add position to current position of the molecule
     */
    void addPosition(jlcxx::ArrayRef<floatType,1> pos) {
        addR({pos[0], pos[1], pos[2]});
    }

    /**
     * add velocity to current velocity of the molecule
     */
    void addVeloctiy(jlcxx::ArrayRef<floatType,1> v) {
        addV({v[0], v[1], v[2]});
    }

    /**
     * add force to current force of the molecule
     */
    void addForce(jlcxx::ArrayRef<floatType,1> f) {
        addF({f[0], f[1], f[2]});
    }

    /**
     * substract force from current force of the molecule
     */
    void subForce(jlcxx::ArrayRef<floatType,1> f) {
        subF({f[0], f[1], f[2]});
    }

    private:
    /**
     * Particle type id.
     */
    // size_t _typeId = 0;
    
    /**
     * Old Force of the particle experiences as 3D vector.
     */
    // std::array<floatType,3> _oldF = {0.0, 0.0, 0.0};
};