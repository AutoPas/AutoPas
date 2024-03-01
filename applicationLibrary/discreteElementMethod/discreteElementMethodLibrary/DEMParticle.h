/**
 * @file DEMParticle.h
 * @author Hoppe (hoppef@hsu-hh.de)
 * @brief Main file for all DEM functionality
 * @version 0.1
 * @date 2022-10-20
 * 
 */

#pragma once

#include <vector>

#include "autopas/particles/Particle.h"
#include "autopas/utils/ExceptionHandler.h"

namespace demLib {

/**
 * @brief Particle class for the DEM functor
 *
 */
class DEMParticle final : public autopas::Particle {
 public:

  /**
   * @brief Default constructor for a new Particle DEM object.
   * 
   */
  DEMParticle() = default;

  /**
   * @brief Construct a new Particle DEM object.
   * 
   * @param pos Particle position
   * @param v Particle velocity
   * @param particleId Particle ID number
   * @param rad Particle radius
   * @param mass Particle mass
   * @param young Particle Young's modulus
   * @param poisson Particle Poisson's ratio
   */
  explicit DEMParticle(std::array<double, 3> pos, std::array<double, 3> v, unsigned long particleId, 
                        double rad = 0.0, double mass = 1.0, double young = 0.0, double poisson=0.0)
      : autopas::Particle(pos, v, particleId), _radius(rad), _mass(mass), _young(young), _poisson(poisson) {}


  /**
   * @brief Destroy the Particle DEM object
   * 
   */
  ~DEMParticle() final = default;


  /**
   * @brief Enums used as ids for accessing and creating a dynamically sized SoA.
   * 
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
    rad,
    mass,
    young,
    poisson,
    typeId,
    ownershipState
  };


  /**
   * @brief The type for the SoA storage.
   * 
   */
  using SoAArraysType = typename autopas::utils::SoAType<
      DEMParticle *, size_t /*id*/, 
      double /*posX*/, double /*posY*/, double /*posZ*/,
      double /*velocityX*/, double /*velocityY*/, double /*velocityZ*/, 
      double /*forceX*/, double /*forceY*/, double /*forceZ*/, 
      double /*oldForceX*/, double /*oldForceY*/, double /*oldForceZ*/, 
      double /*rad*/, double /*mass*/, double /*young*/, double /*poisson*/,
      size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;


  /**
   * @brief Non-const getter for the pointer of this object.
   * 
   * @tparam attribute Attribute name.
   * @return this. 
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }

  /**
   * @brief Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * 
   * @tparam attribute Attribute name of the requested attribute.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
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
    } else if constexpr (attribute == AttributeNames::rad) {
      return getRad();
    } else if constexpr (attribute == AttributeNames::mass) {
      return getMass();
    } else if constexpr (attribute == AttributeNames::young) {
      return getYoung();
    } else if constexpr (attribute == AttributeNames::poisson) {
      return getPoisson();
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("DEMParticle::get() unknown attribute {}", attribute);
    }
  }

  /**
   * @brief Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * 
   * @tparam attribute attribute Attribute name.
   * @param value New value of the requested attribute.
   * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
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
    } else if constexpr (attribute == AttributeNames::rad) {
      _radius = value;
    } else if constexpr (attribute == AttributeNames::mass) {
      _mass = value;
    } else if constexpr (attribute == AttributeNames::young) {
      _young = value;
    } else if constexpr (attribute == AttributeNames::poisson) {
      _poisson = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      setTypeId(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("DEMParticle::set() unknown attribute {}", attribute);
    }
  }

  /**
   * @brief Get the old force
   * 
   * @return Force of previous step
   */
  [[nodiscard]] std::array<double, 3> getOldF() const { return _oldF; }

  /**
   * @brief Set the old force
   * 
   * @param New force of previous time step
   */
  void setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

  /**
   * @brief Get the particle radius
   * 
   * @return Current particle radius
   */
  [[nodiscard]] double getRad() const { return _radius; }

  /**
   * @brief Set the particle radius
   * 
   * @param New particle radius
   */
  void setRad(const double &radius) { _radius = radius; }

  /**
   * @brief Get the particle mass
   * 
   * @return Current particle mass
   */
  [[nodiscard]] double getMass() const { return _mass; }

  /**
   * @brief Set the particle mass
   * 
   * @param New particle mass
   */
  void setMass(const double &mass) { _mass = mass; }

  /**
   * @brief Get the Young's modulus of the particle
   * 
   * @return Current Young's modulus
   */
  [[nodiscard]] double getYoung() const { return _young; }

  /**
   * @brief Set the Young's modulus of the particle
   * 
   * @param New Young's modulus
   */
  void setYoung(const double &young) { _young = young; }

  /**
   * @brief Get the particle's Poisson ratio
   * 
   * @return Current Poisson ratio
   */
  [[nodiscard]] double getPoisson() const { return _poisson; }

  /**
   * @brief Set the particle's Poisson ratio
   * 
   * @param New Poisson ratio
   */
  void setPoisson(const double &poisson) { _poisson = poisson; }

  [[nodiscard]] size_t getTypeId() const { return _typeId; }

  void setTypeId(const size_t & typeId) { _typeId = typeId; }


 private:
  
  /**
   * @brief Particle type ID
   * 
   */
  size_t _typeId = 0;
  
  /**
   * @brief Particle radius
   * 
   */
  double _radius = 0.0;

  /**
   * @brief Particle mass
   * 
   */
  double _mass = 1.0;

  /**
   * @brief Young's modulus E of the particle
   * 
   * Mechanical property measuring the tensile or compressive stiffness of a solid material 
   * when force is applied lengthwise.
   * Quantifies the relationship between tensile/compressive stress sigma (force per unit area) 
   * and axial strain epsilon (proportional deformation).
   * Higher E means stiffer material.
   * 
   */
  double _young = 0.0;

  /**
   * @brief Poisson's ratio nu
   * 
   * Measure of the deformation (expansion or contraction) of a material in directions 
   * perpendicular to the specific direction of loading.
   * Values range between 0 (= no deformation) to 0.5 (= strong/easy deformation).
   * Solid materials are in the range of nu = 0.2 to 0.3.
   * 
   */
  double _poisson = 0.0;

  /**
   * @brief Old Force of the particle experiences as 3D vector.
   * 
   */
  std::array<double, 3> _oldF = {0., 0., 0.};

};

}  // namespace demLib
