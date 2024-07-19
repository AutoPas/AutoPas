/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>

#include "autopas/particles/Particle.h"
#include "autopas/utils/ExceptionHandler.h"

namespace mdLib {

/**
 * Molecule class for the LJFunctor.
 */
class MoleculeLJ_NoPPL : public autopas::Particle {
 public:
  MoleculeLJ_NoPPL() = default;

  /**
   * Constructor of lennard jones molecule with initialization of typeID.
   * @param pos Position of the molecule.
   * @param v Velocity of the molecule.
   * @param moleculeId Unique Id of the molecule.
   * @param typeId TypeId of the molecule.
   * @param squareRootEpsilon sqrt(epsilon of molecule)
   * @param sigmaDiv2 sigma of molecule/2
   */
  MoleculeLJ_NoPPL(const std::array<double, 3> &pos, const std::array<double, 3> &v, unsigned long moleculeId,
                   double squareRootEpsilon = 1., double sigmaDiv2 = 0.5);

  ~MoleculeLJ_NoPPL() = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int {
    ptr,
    id,
    posX,
    posY,
    posZ,
    velocityX,
    velocityY,
    velocityZ,
    forceX,
    forceY,
    forceZ,
    oldForceX,
    oldForceY,
    oldForceZ,
    squareRootEpsilon,
    sigmaDiv2,
    mass,
    ownershipState
  };

  /**
   * The type for the SoA storage.
   *
   * @note The attribute owned is of type float but treated as a bool.
   * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
   * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<MoleculeLJ_NoPPL *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/,
                                       double /*vx*/, double /*vy*/, double /*vz*/, double /*fx*/, double /*fy*/,
                                       double /*fz*/, double /*oldFx*/, double /*oldFy*/, double /*oldFz*/,
                                       double /*squareRootEpsilon*/, double /*sigmaDiv2*/, double /*mass*/,
                                       autopas::OwnershipState /*ownershipState*/>::Type;

  /**
   * Non-const getter for the pointer of this object.
   * @tparam attribute Attribute name.
   * @return this.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }
  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   * @note Moving this function to the .cpp leads to undefined references
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::velocityX) {
      return getV()[0];
    } else if constexpr (attribute == AttributeNames::velocityY) {
      return getV()[1];
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      return getV()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      return getOldF()[0];
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      return getOldF()[1];
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      return getOldF()[2];
    } else if constexpr (attribute == AttributeNames::squareRootEpsilon) {
      return getSquareRootEpsilon();
    } else if constexpr (attribute == AttributeNames::sigmaDiv2) {
      return getSigmaDiv2();
    } else if constexpr (attribute == AttributeNames::mass) {
      return getMass();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
   * @note Moving this function to the .cpp leads to undefined references
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::id) {
      setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::velocityX) {
      _v[0] = value;
    } else if constexpr (attribute == AttributeNames::velocityY) {
      _v[1] = value;
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      _v[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::oldForceX) {
      _oldF[0] = value;
    } else if constexpr (attribute == AttributeNames::oldForceY) {
      _oldF[1] = value;
    } else if constexpr (attribute == AttributeNames::oldForceZ) {
      _oldF[2] = value;
    } else if constexpr (attribute == AttributeNames::squareRootEpsilon) {
      _squareRootEpsilon = value;
    } else if constexpr (attribute == AttributeNames::sigmaDiv2) {
      _sigmaDiv2 = value;
    } else if constexpr (attribute == AttributeNames::mass) {
      _mass = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the old force.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getOldF() const;

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<double, 3> &oldForce);

  /** ---- Mixing Parameters (Non-Efficient) ---- **/

  /**
   * Get epsilon.
   * @warning this molecule stores sqrt(epsilon) and so this function performs a square on the stored result. Use
   * getSquareRootEpsilon in performance critical regions which require the square root of epsilon.
   * @return epsilon
   */
  [[nodiscard]] double getEpsilon() const;

  /**
   * Set epsilon.
   * @warning this molecule stores sqrt(epsilon) and so this function performs a square root and should not be used in
   * performance critical regions.
   * @param epsilon
   */
  void setEpsilon(const double &epsilon);

  /**
   * Get sigma.
   * @warning this molecule stores sigma/2 and so this function performs a multiplication on the stored result. Use
   * getSigmaDiv2 in performance critical regions which require sigma/2.
   * @return sigma
   */
  [[nodiscard]] double getSigma() const;

  /**
   * Set sigma.
   * @warning this molecule stores sigma/2 and so this function performs a division and should not be used in
   * performance critical regions.
   * @param sigma
   */
  void setSigma(const double &sigma);

  /** ---- Mixing Parameters (Efficient) ---- **/

  /**
   * Get sqrt(epsilon).
   * @return sqrt(epsilon)
   */
  [[nodiscard]] const double &getSquareRootEpsilon() const;

  /**
   * Set sqrt(epsilon)
   * @param squareRootEpsilon
   */
  void setSquareRootEpsilon(const double &squareRootEpsilon);

  /**
   * Get sigma/2
   * @return sigma/2
   */
  [[nodiscard]] const double &getSigmaDiv2() const;

  /**
   * Set sigma/2
   * @param sigmaDiv2
   */
  void setSigmaDiv2(const double &sigmaDiv2);

  /**
   * Get mass
   * @return mass
   */
  [[nodiscard]] const double &getMass() const;

  /**
   * Set mass
   * @param mass
   */
  void setMass(const double &mass);

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const override;

 protected:
  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> _oldF = {0., 0., 0.};

 private:
  /**
   * sqrt(epsilon). Used directly in force calculation to avoid square roots in kernel.
   */
  double _squareRootEpsilon{1.};

  /**
   * sigma/2. Used directly in force calculation to avoid divisions in kernel.
   */
  double _sigmaDiv2{0.5};

  /**
   * mass of molecule
   */
  double _mass{1.};
};

}  // namespace mdLib
