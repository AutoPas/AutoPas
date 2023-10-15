/**
* @file AbsoluteMultiSiteMoleculeLJ.h
* @date 10/10/2023
* @author Johannes Riemenschneider
 */

#pragma once

#include "MoleculeLJ.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/particles/ParticleBase.h"

namespace mdLib {
/**
* Standard multi-site LJ molecules.
*
* The molecule is treated as a single particle for the purposes of cutoffs and containers, with a quaternion for
* angular direction, a 3D vector-array for angular velocity, and a vectors of site positions relative to the center of
* mass and angular direction.
*
*/
class AbsoluteMultiSiteMoleculeLJ : public mdLib::MoleculeLJ {
 using idType = size_t;

public:
 AbsoluteMultiSiteMoleculeLJ() = default;

 /**
  * Constructor of the AbsoluteMultiSiteMoleculeLJ Class
  * This Constructor does NOT initialize the absolute Site positions. Therefore these need to be initialized later in the program
  * @param r Position of the particle.
  * @param v Velocity of the particle.d
  * @param angularVel Rotational velocity of the particle.
  * @param moleculeId Id of the particle.
  * @param typeId Id of the type of the particle. Used in conjunction with ParticlePropertiesLibrary to access
  * molecular information such as site types and relative site positions.
  */
 AbsoluteMultiSiteMoleculeLJ(std::array<double, 3> r, std::array<double, 3> v,
                     std::array<double, 3> angularVel, unsigned long moleculeId, unsigned long typeId = 0);

 /**
  * Destructor of the AbsoluteMultiSiteMoleculeLJ class.
  */
 ~AbsoluteMultiSiteMoleculeLJ() override = default;

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
   //quaternion0,
   //quaternion1,
   //quaternion2,
   //quaternion3,
   absoluteSitePositionsX,
   absoluteSitePositionsY,
   absoluteSitePositionsZ,
   angularVelX,
   angularVelY,
   angularVelZ,
   torqueX,
   torqueY,
   torqueZ,
   typeId,
   ownershipState
 };

 /**
  * The type for the SoA storage.
  *
  * @note The attribute owned is of type float but treated as a bool.
  * This means it shall always only take values 0.0 (=false) or 1.0 (=true).
  * The reason for this is the easier use of the value in calculations (See LJFunctor "energyFactor")
  */
 // clang-format off
 using SoAArraysType = typename autopas::utils::SoAType<
     AbsoluteMultiSiteMoleculeLJ *,
     size_t, // id
     double, // x
     double, // y
     double, // z
     double, // vx
     double, // vy
     double, // vz
     double, // fx
     double, // fy
     double, // fz
     double, // oldFx
     double, // oldFy
     double, // oldFz
     //double, // q0
     //double, // q1
     //double, // q2
     //double, // q3
     std::vector<double>, //absSitePosX
     std::vector<double>, //absSitePosY
     std::vector<double>, //absSitePosZ
     double, // angVx
     double, // angVy
     double, // angVz
     double, // tx
     double, // ty
     double, // tz
     size_t, // typeid
     autopas::OwnershipState //ownerState
 >::Type;
 // clang-format on

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
     //}else if constexpr (attribute == AttributeNames::quaternion0) {
     //  return getQuaternion()[0];
     //} else if constexpr (attribute == AttributeNames::quaternion1) {
     //  return getQuaternion()[1];
     //} else if constexpr (attribute == AttributeNames::quaternion2) {
     //  return getQuaternion()[2];
     //} else if constexpr (attribute == AttributeNames::quaternion3) {
     //  return getQuaternion()[3];
   }else if constexpr (attribute == AttributeNames::absoluteSitePositionsX){
     return getAbsoluteSitePositionsX();
   }else if constexpr (attribute == AttributeNames::absoluteSitePositionsY){
     return getAbsoluteSitePositionsY();
   }else if constexpr (attribute == AttributeNames::absoluteSitePositionsZ){
     return getAbsoluteSitePositionsZ();
   } else if constexpr (attribute == AttributeNames::angularVelX) {
     return getAngularVel()[0];
   } else if constexpr (attribute == AttributeNames::angularVelY) {
     return getAngularVel()[1];
   } else if constexpr (attribute == AttributeNames::angularVelZ) {
     return getAngularVel()[2];
   } else if constexpr (attribute == AttributeNames::torqueX) {
     return getTorque()[0];
   } else if constexpr (attribute == AttributeNames::torqueY) {
     return getTorque()[1];
   } else if constexpr (attribute == AttributeNames::torqueZ) {
     return getTorque()[2];
   } else if constexpr (attribute == AttributeNames::typeId) {
     return getTypeId();
   } else if constexpr (attribute == AttributeNames::ownershipState) {
     return this->_ownershipState;
   } else {
     autopas::utils::ExceptionHandler::exception("AbsoluteMultiSiteMoleculeLJ::get() unknown attribute {}", attribute);
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
   //} else if constexpr (attribute == AttributeNames::quaternion0) {
   //  _q[0] = value;
   //} else if constexpr (attribute == AttributeNames::quaternion1) {
   //  _q[1] = value;
   //} else if constexpr (attribute == AttributeNames::quaternion2) {
   //  _q[2] = value;
   //} else if constexpr (attribute == AttributeNames::quaternion3) {
   //  _q[3] = value;
   } else if constexpr (attribute == AttributeNames::absoluteSitePositionsX) {
     _absSitePositionsX = value;
   } else if constexpr (attribute == AttributeNames::absoluteSitePositionsY) {
     _absSitePositionsY = value;
   } else if constexpr (attribute == AttributeNames::absoluteSitePositionsZ) {
     _absSitePositionsZ = value;
   } else if constexpr (attribute == AttributeNames::angularVelX) {
     _angularVel[0] = value;
   } else if constexpr (attribute == AttributeNames::angularVelY) {
     _angularVel[1] = value;
   } else if constexpr (attribute == AttributeNames::angularVelZ) {
     _angularVel[2] = value;
   } else if constexpr (attribute == AttributeNames::torqueX) {
     _torque[0] = value;
   } else if constexpr (attribute == AttributeNames::torqueY) {
     _torque[1] = value;
   } else if constexpr (attribute == AttributeNames::torqueZ) {
     _torque[2] = value;
   } else if constexpr (attribute == AttributeNames::typeId) {
     _typeId = value;
   } else if constexpr (attribute == AttributeNames::ownershipState) {
     this->_ownershipState = value;
   } else {
     autopas::utils::ExceptionHandler::exception("AbsoluteMultiSiteMoleculeLJ::set() unknown attribute {}", attribute);
   }
 }

 ///**
 // * Get the quaternion defining rotation
 // * @return quaternion defining rotation
 // */
 //[[nodiscard]] const std::array<double, 4> &getQuaternion() const;

 ///**
 // * Set the quaternion defining rotation
 // * @param q quaternion defining rotation
 // */
 //void setQuaternion(const std::array<double, 4> &q);

 /**
  * Set the x-component of the absoluteSitePositions
  */
 void setAbsoluteSitePositionsX(const std::vector<double> sitesX);

 /**
  * Set the y-component of the absoluteSitePositions
  */
 void setAbsoluteSitePositionsY(const std::vector<double> sitesY);

 /**
  * Set the z-component of the absoluteSitePositions
  */
 void setAbsoluteSitePositionsZ(const std::vector<double> sitesZ);


 /**
  * Returns the number of sites of that molecule.
  */
  int getNumberOfSites();

 /**
  * Get the x-components of all Sites in that molecule.
  * @return Vector containing the absolute x-positions of all Sites
  */
 [[nodiscard]] const std::vector<double> &getAbsoluteSitePositionsX() const;

 /**
  * Get the y-components of all Sites in that molecule.
  * @return Vector containing the absolute y-positions of all Sites
  */
 [[nodiscard]] const std::vector<double> &getAbsoluteSitePositionsY() const;

 /**
  * Get the z-components of all Sites in that molecule.
  * @return Vector containing the absolute z-positions of all Sites
  */
 [[nodiscard]] const std::vector<double> &getAbsoluteSitePositionsZ() const;

 /**
  * Get the angular velocity
  * @return angular velocity
  */
 [[nodiscard]] const std::array<double, 3> &getAngularVel() const;

 /**
  * Set the angular velocity
  * @param angularVel
  */
 void setAngularVel(const std::array<double, 3> &angularVel);

 /**
  * Adds given angular velocity to the particle's angular velocity.
  * @param angularVel angular velocity to be added
  */
 void addAngularVel(const std::array<double, 3> &angularVel);

 /**
  * Get the torque.
  * @return torque
  */
 [[nodiscard]] const std::array<double, 3> &getTorque() const;

 /**
  * Set the torque.
  * @param torque
  */
 void setTorque(const std::array<double, 3> &torque);

 /**
  * Adds given torque to the particle's torque.
  * @param torque torque to be added
  */
 void addTorque(const std::array<double, 3> &torque);

 /**
  * Subracts given torque to the particle's torque.
  * @param torque torque to be subtracted
  */
 void subTorque(const std::array<double, 3> &torque);

 /**
  * Creates a string containing all data of the particle.
  * @return String representation.
  */
 [[nodiscard]] std::string toString() const override;

protected:

 /**
  * Helper function for toString method
  * @return
  */
 template<typename T> std::string  vectorToString(std::vector<T> v) const;

 // /**
 // * Rotational direction of particle as quaternion.
 // */
 //std::array<double, 4> _q{};

 /**
  * Vector storing the x-component the absolute Site positions.
  * Since the number of sites in this molecule stay constant the size of this vector will also remain constant after initialization.
  */
 std::vector<double> _absSitePositionsX{};

 /**
  * Vector storing the y-component the absolute Site positions.
  * Since the number of sites in this molecule stay constant the size of this vector will also remain constant after initialization.
  */

 std::vector<double> _absSitePositionsY{};
 /**
  * Vector storing the z-component the absolute Site positions.
  * Since the number of sites in this molecule stay constant the size of this vector will also remain constant after initialization.
  */
 std::vector<double> _absSitePositionsZ{};

 /**
  * Angular velocity of the particle
  */
 std::array<double, 3> _angularVel{};

 /**
  * Torque applied to particle.
  */
 std::array<double, 3> _torque{};
};

}  // namespace mdLib