/**
 * @file MoleculeLJ.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <vector>
#include "autopas/particles/Particle.h"

namespace autopas {

/**
 * lennard jones molecule class
 */
class MoleculeLJ : public Particle {
 public:
  MoleculeLJ() = default;

  /**
   * constructor of a lennard jones molecule
   * @param r position of the molecule
   * @param v velocity of the molecule
   * @param id id of the molecule
   */
  explicit MoleculeLJ(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : Particle(r, v, id) {}

  ~MoleculeLJ() override = default;

  /**
   * get epsilon (characteristic energy of the lj potential)
   * @return epsilon
   */
  static double getEpsilon() { return EPSILON; }

  /**
   * set epsilon (characteristic energy of the lj potential)
   * @param epsilon
   */
  static void setEpsilon(double epsilon) { EPSILON = epsilon; }

  /**
   * get sigma (characteristic length of the lj potential)
   * @return sigma
   */
  static double getSigma() { return SIGMA; }

  /**
   * set sigma (characteristic length of the lj potential)
   * @param sigma
   */
  static void setSigma(double sigma) { SIGMA = sigma; }

  /**get mass
   * @return MASS
   */
  static double getMass() { return MASS; }

  /**set mass
   * @param mass
   */
  static void setMass(double mass) { MASS = mass; }


    /**
     * Enums used as ids for accessing and creating a dynamically sized SoA.
     */
    enum AttributeNames : int { id, typeId, posX, posY, posZ, forceX, forceY, forceZ, owned };

  /**
   * the type for the soa storage
   */
    typedef autopas::utils::SoAType<size_t,size_t, double, double, double, double, double, double, double>::Type SoAArraysType;


    /**
     * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
     * @tparam attribute Attribute name.
     * @return Value of the requested attribute.
     * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
     */
    template <AttributeNames attribute>
    constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
        switch (attribute) {
            case AttributeNames::id:
                return getID();
            case AttributeNames::typeId:
                return getTypeId();
            case AttributeNames::posX:
                return getR()[0];
            case AttributeNames::posY:
                return getR()[1];
            case AttributeNames::posZ:
                return getR()[2];
            case AttributeNames::forceX:
                return getF()[0];
            case AttributeNames::forceY:
                return getF()[1];
            case AttributeNames::forceZ:
                return getF()[2];
            case AttributeNames::owned:
                return isOwned() ? 1. : 0.;
            default:
                utils::ExceptionHandler::exception("ParticleBase::get: unknown attribute");
                return 0;
        }
    }

    /**
     * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
     * @tparam attribute Attribute name.
     * @param value New value of the requested attribute.
     * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
     */
    template <AttributeNames attribute>
    constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
        switch (attribute) {
            case AttributeNames::id:
                setID(value);
                break;
            case AttributeNames::typeId:
                setTypeId(value);
                break;
            case AttributeNames::posX:
                _r[0] = value;
                break;
            case AttributeNames::posY:
                _r[1] = value;
                break;
            case AttributeNames::posZ:
                _r[2] = value;
                break;
            case AttributeNames::forceX:
                _f[0] = value;
                break;
            case AttributeNames::forceY:
                _f[1] = value;
                break;
            case AttributeNames::forceZ:
                _f[2] = value;
                break;
            case AttributeNames::owned:
                setOwned(value == 1.);
                break;
        }
    }


  /**get OldForce
   * @return OLDF
   */
  std::array<double, 3> getOldf() const { return OLDF; }

  /**set OldForce
   * @param oldf
   */
  void setOldf(const std::array<double, 3> &oldf) { OLDF = oldf; }

  /**get TypeId
   * @return _typeId
   * */
  size_t getTypeId() const { return _typeId; }
  /**set _TypeId of Particle
   * @param typeId
   * */
  void setTypeId(size_t typeId) { _typeId = typeId; }

 private:
  static double EPSILON, SIGMA, MASS;

  /**
   * Particle type id.
   */
  size_t _typeId = 0;

 private:
  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> OLDF = {0., 0., 0.};
};

}  // namespace autopas