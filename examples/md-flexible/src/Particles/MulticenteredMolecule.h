/**
 * @file MulticenteredParticleBase.h
 * @date 14/02/2022
 * @author S. Newcome
 */

#pragma once

#include "autopas/particles/ParticleBase.h"
#include "autopas/particles/OwnershipState.h"

/**
 * Standard multi-centre particle type. Consists of 4 site types: Lennard-Jones, Coulomb, Dipole, Quadrupole
 *
 * The molecule is treated as a single particle for the purposes of cutoffs and containers, with a quaternion for
 * angular direction, a 3D vector-array for angular velocity, and a vectors of site positions relative to the CoM and
 * angular direction.
 *
 * @tparam floatType
 * @tparam idType
 */
template <typename floatType, typename idType>
class MulticenteredMolecule : autopas::ParticleBase<floatType, idType> {
 public:
  MulticenteredMolecule() = default;

  /**
   * Constructor of the MulticenteredParticle Class
   * @param r Position of the particle.
   * @param v Velocity of the particle.
   * @param q Quaternion defining rotation of particle.
   * @param D Rotational velocity of the particle.
   * @param sites Vector of sites of the particle.
   * @param id Id of the particle.
   */
   MulticenteredMolecule(std::array<floatType, 3> r, std::array<floatType, 3> v, std::array<floatType, 4> q,
                            std::array<floatType, 3> D, std::vector<std::array<floatType,3>> sitePosLJ,
                            std::vector<std::array<floatType,3>> sitePosC, std::vector<std::array<floatType,3>> sitePosD,
                            std::vector<std::array<floatType,3>> sitePosQ, idType id)
                            : _r(r), _v(v), _q(q), _D(D), _sitePosLJ(sitePosLJ), _sitePosC(sitePosC),
                              _sitePosD(sitePosD), _sitePosQ(sitePosQ, _id) {}

   /**
    * Destructor of the MulticenteredParticle class.
    */
    virtual ~MulticenteredMolecule() = default;

   protected:
   /**
    * (Centre of) Particle position as 3D coords
    */
   std::array<floatType,3> _r;

   /**
    * Velocity of particle.
    */
   std::array<floatType,3> _v;

   /**
    * Force experienced by particle.
    */
   std::array<floatType,3> _f;

   /**
    * Rotational direction of particle as quaternion.
    */
   std::array<floatType,4> _q;

   /**
    * Rotational velocity of the particle
    */
   std::array<floatType,3> _D;

   /**
    * Angular momentum of the particle.
    */


   /**
    * Net torque applied to particle. (+ve = counterclockwise)
    */
   floatType _t;

   /**
    * Position of Lennard-Jones sites relative to CoM.
    */
   std::vector<std::array<floatType,3>> _sitePosLJ;

   /**
    * Position of Coulombic sites relative to CoM.
    */
   std::vector<std::array<floatType,3>> _sitePosC;

   /**
    * Position of Dipole sites relative to CoM.
    */
   std::vector<std::array<floatType,3>> _sitePosD;

   /**
    * Position of Quadrupole sites relative to CoM.
    */
   std::vector<std::array<floatType,3>> _sitePosQ;

   /**
    * Particle id.
    */
   idType _id;

   /**
   * Defines the state of the ownership of the particle.
    */
   autopas::OwnershipState _ownershipState{autopas::OwnershipState::owned};

  public:
   /**
   * get the force acting on the particle
   * @return force
    */
   [[nodiscard]] const std::array<floatType, 3> &getF() const { return _f; }

   /**
   * Set the force acting on the particle
   * @param f force
    */
   void setF(const std::array<floatType, 3> &f) { _f = f; }

   /**
   * Add a partial force to the force acting on the particle
   * @param f partial force to be added
    */
   void addF(const std::array<floatType, 3> &f) { _f = autopas::utils::ArrayMath::add(_f, f); }

   /**
   * Substract a partial force from the force acting on the particle
   * @param f partial force to be substracted
    */
   void subF(const std::array<floatType, 3> &f) { _f = autopas::utils::ArrayMath::sub(_f, f); }

   /**
   * Get the id of the particle
   * @return id
    */
   idType getID() const { return _id; }

   /**
   * Set the id of the particle
   * @param id id
    */
   void setID(idType id) { _id = id; }

   /**
   * Get the position of the particle
   * @return current position
    */
   [[nodiscard]] const std::array<floatType, 3> &getR() const { return _r; }

   /**
   * Set the position of the particle
   * @param r new position
    */
   void setR(const std::array<floatType, 3> &r) { _r = r; }

   /**
   * Add a distance vector to the position of the particle
   * @param r vector to be added
    */
   void addR(const std::array<floatType, 3> &r) { _r = autopas::utils::ArrayMath::add(_r, r); }

   /**
   * Get the velocity of the particle
   * @return current velocity
    */
   [[nodiscard]] const std::array<floatType, 3> &getV() const { return _v; }

   /**
   * Set the velocity of the particle
   * @param v new velocity
    */
   void setV(const std::array<floatType, 3> &v) { _v = v; }

   /**
   * Add a vector to the current velocity of the particle
   * @param v vector to be added
    */
   void addV(const std::array<floatType, 3> &v) { _v = autopas::utils::ArrayMath::add(_v, v); }

   /**
   * Get the quaternion defining rotation
   * @return quaternion defining rotation
    */
   [[nodiscard]] const std::array<floatType, 4> &getQ() const { return _q; }

   /**
   * Set the quaternion defining rotation
   * @param q quaternion defining rotation
    */
   void setQ(const std::array<floatType, 4> &q) { _q = q; }

   /**
   * Get the rotational velocity
   * @return rotational velocity
    */
   [[nodiscard]] const std::array<floatType, 3> &getD() const { return _D; }

   /**
   * Set the rotational velocity
   * @param D rotational velocity
    */
   void setD(const std::array<floatType, 3> &D) { _D = D; }

   /**
    * Add Lennard-Jones site
    * @param sitePos relative position of LJ site
    */
    void addSiteLJ(const std::array<floatType,3> &sitePos) {_sitePosLJ.push_back(sitePos); }

    /**
    * Set Lennard-Jones sites
    * @param sitePosLJ relative position of LJ sites
     */
    void setSitesLJ(const std::vector<std::array<floatType,3>> &sitePosLJ) {_sitePosLJ = sitePosLJ; }

    /**
     * Get Lennard-Jones site relative positions.
     * @return LJ site relative positions
     */
     [[nodiscard]] const std::vector<std::array<floatType,3>> &getSitesLJ() const { return _sitePosLJ; }

    /**
    * Add Coulombic site
    * @param relative position of Coulombic site
     */
    void addSiteC(const std::array<floatType,3> &sitePos) {_sitePosC.push_back(sitePos); }

    /**
    * Set Coulombic sites
    * @param sitePosC relative position of Coloumbic sites
     */
    void setSitesC(const std::vector<std::array<floatType,3>> &sitePosC) {_sitePosC = sitePosC; }

    /**
     * Get Coulombic site relative positions.
     * @return Coulombic site relative positions
     */
    [[nodiscard]] const std::vector<std::array<floatType,3>> &getSitesC() const { return _sitePosC; }

    /**
    * Add Dipole site
    * @param relative position of Dipole site
     */
    void addSiteD(const std::array<floatType,3> &sitePos) {_sitePosD.push_back(sitePos); }

    /**
    * Set Dipole sites
    * @param sitePosD relative position of Dipole sites
     */
    void setSitesD(const std::vector<std::array<floatType,3>> &sitePosD) {_sitePosD = sitePosD; }

    /**
     * Get Dipole site relative positions.
     * @return Dipole site relative positions
     */
    [[nodiscard]] const std::vector<std::array<floatType,3>> &getSitesD() const { return _sitePosD; }

    /**
    * Add Quadrupole site
    * @param relative position of Quadrupole site
     */
    void addSiteQ(const std::array<floatType,3> &sitePos) {_sitePosQ.push_back(sitePos); }

    /**
    * Set Quadrupole sites
    * @param sitePosQ relative position of Quadrupole sites
     */
    void setSitesQ(const std::vector<std::array<floatType,3>> &sitePosQ) {_sitePosQ = sitePosQ; }

    /**
     * Get Quadrupole site relative positions.
     * @return Quadrupole site relative positions
     */
    [[nodiscard]] const std::vector<std::array<floatType,3>> &getSitesQ() const { return _sitePosQ; }

    /**
   * Creates a string containing all data of the particle.
   * @return String representation.
     */
    [[nodiscard]] virtual std::string toString() const {
      std::ostringstream text;
      std::ostringstream lj_str, c_str, d_str, q_str;
      if (_sitePosLJ.empty) {lj_str << ""; } else {
        lj_str << "\n  Lennard-Jones: " << autopas::utils::ArrayUtils::to_string(_sitePosLJ);
      }
      if (_sitePosC.empty) {c_str << ""; } else {
        c_str << "\n  Coulomb       : " << autopas::utils::ArrayUtils::to_string(_sitePosC);
      }
      if (_sitePosD.empty) {d_str << ""; } else {
        d_str << "\n  Dipole        : " << autopas::utils::ArrayUtils::to_string(_sitePosD);
      }
      if (_sitePosLJ.empty) {lj_str << ""; } else {
        q_str << "\n  Quadrupole    : " << autopas::utils::ArrayUtils::to_string(_sitePosQ);
      }
      // clang-format off
      text << "Particle"
         << "\nID                 : " << _id
         << "\nPosition           : "
         << autopas::utils::ArrayUtils::to_string(_r)
         << "\nVelocity           : "
         << autopas::utils::ArrayUtils::to_string(_v)
         << "\nForce              : "
         << autopas::utils::ArrayUtils::to_string(_f)
         << "\nQuaternion         : "
         << autopas::utils::ArrayUtils::to_string(_q)
         << "\nRotational Velocity: "
         << autopas::utils::ArrayUtils::to_string(_D)
         << "\nRelative Site Positions:"
         << lj_str << c_str << d_str << q_str
         << "\nOwnershipState     : "
         << _ownershipState;
      // clang-format on
      return text.str();
    }

    using SoAArraysType =
        typename autopas::utils::SoAType<MulticenteredMolecule<floatType,idType> *, idType /*id*/ , floatType /*x*/,
          floatType /*y*/, floatType /*z*/, floatType /*fx*/, floatType /*fy*/, floatType /*fz*/,
          autopas::OwnershipState /*ownershipState*/>::Type;


};

