/**
 * @file MoleculeInterface.h
 * @date 10/03/2022
 * @author S. Newcome
 */

#pragma once

#include "autopas/particles/Particle.h"

namespace MDLibrary {

/**
 * This interface serves as a common parent class for all particle types provided in the Molecular Dynamics Library for AutoPas
 */

class MoleculeInterface : public autopas::Particle {
 public:
  MoleculeInterface() = default;

  /**
   * Constructor of molecule, with initialization of typeID
   * @param pos Position of molecule.
   * @param vel Velocity of molecule.
   * @param moleculeID Unique id of the molecule.
   * @param typeID Id of type of molecule (dictating any features e.g. lennard-jones parameters, site positions)
   */
  MoleculeInterface(std::array<double,3> pos, std::array<double,3> vel, unsigned long moleculeId, unsigned long typeId = 0)
      : autopas::Particle(pos, vel, moleculeId), _typeId(typeId) {}

  /**
   * Destructor of the Molecule Interface
   */
  virtual ~MoleculeInterface() = default;



  /**
   * Get the old force. (Used for some integrator)
   * @return old force.
   */
  [[nodiscard]] virtual std::array<double, 3> getOldF() const = 0;

  /**
   * Set the old force. (Used for some integrator)
   * @param oldForce
   */
   virtual void setOldF(const std::array<double, 3> &oldForce) {}

   /**
   * Get the quaternion defining rotation. (Used for rotation-dependant simulations)
   * @return quaternion defining rotation
    */
   [[nodiscard]] virtual const std::array<double, 4> &getQ() = 0;

   /**
   * Set the quaternion defining rotation. (Used for rotation-dependant simulations)
   * @param q quaternion defining rotation
    */
   virtual void setQ(const std::array<double, 4> &q) {}

   /**
   * Get the angular velocity. (Used for rotation-dependant simulations)
   * @return angular velocity
    */
   [[nodiscard]] virtual const std::array<double, 3> &getAngularVel() = 0;

   /**
   * Set the angular velocity. (Used for rotation-dependant simulations)
   * @param angularVelocity
    */
   virtual void setAngularVel(const std::array<double, 3> &angularVel) {}

   /**
    * Get TypeId. Used for differentiate between different types of molecules.
    * @return typeId
    */
   [[nodiscard]] size_t getTypeId() const { return _typeId; }

   /**
   * Set the type id of the Molecule. Used for differentiate between different types of molecules.
   * @param typeId
    */
   void setTypeId(size_t typeId) { _typeId = typeId; }

  private:
   /**
    * Particle type id.
    */
   size_t _typeId = 0;

};
}
