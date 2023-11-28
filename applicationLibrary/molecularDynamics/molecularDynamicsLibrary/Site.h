//
// Created by johnny on 14.11.23.
//

#pragma once
#include "MoleculeLJ.h"

namespace mdLib {
class Site : public autopas::Particle{
 public:

  Site() = default;

 /**
 * Constructor initializing Site completely
 * @param pos Absolute position of the site
 * @param v Velocity of the site
 * @param siteId Id uniquely identifying the site
 * @param typeId TypeId of the site.
 * @param moleculeId Identifier of the molecule this site is a part of. Sites are in a many-to-one relationship with Molecules.
   */
  Site(const std::array<double, 3> &pos,const std::array<double, 3> &v, unsigned long siteId,
             unsigned long moleculeId, unsigned long indexInsideMolecule);

  void setMoleculeId(unsigned long moleculeId);

  [[nodiscard]] unsigned long getMoleculeId() const;

  [[nodiscard]] size_t getIndexInsideMolecule() const;

  void setIndexInsideMolecule(size_t indexInsideMolecule);

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int {
    ptr,
    id,
    owningMoleculeId,
    indexInsideMolecule,
    posX,
    posY,
    posZ,
    typeId,
    forceX,
    forceY,
    forceZ,
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
      typename autopas::utils::SoAType<Site *, size_t /*id*/, unsigned long /*owningMoleculeId*/, unsigned long /*indexInsideMolecule*/,
                                       double /*x*/, double /*y*/, double /*z*/, unsigned long /*typeid*/,
                                       double /*fx*/, double /*fy*/, double /*fz*/,
                                        autopas::OwnershipState /*ownershipState*/>::Type;

  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::owningMoleculeId){
      return getMoleculeId();
    } else if constexpr (attribute == AttributeNames::indexInsideMolecule){
      return getIndexInsideMolecule();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("Site::get() unknown attribute {}", attribute);
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
    }else if constexpr (attribute == AttributeNames::indexInsideMolecule){
      setIndexInsideMolecule(value);
    } else if constexpr (attribute == AttributeNames::owningMoleculeId){
      setMoleculeId(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::typeId) {
      setTypeId(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("Site::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get TypeId.
   * @return
   */
  [[nodiscard]] unsigned long getTypeId() const;

  /**
   * Set typeId of this individual Site.
   * @param typeId
   */
  void setTypeId(unsigned long typeId);

 private:
  /**
   * The moleculeId identifies the MultiSite Molecule this Site belongs to. It corresponds to the moleculeId of MultiSiteMoleculeLJ
   */
  unsigned long _moleculeID;

  /**
   * The indexInsideMolecule defines which Site this site is inside of its respective molecule.
   * For example: particlePropertiesLibrary.getSitePositions(_moleculeID)[_indexInsideMolecule] gets SitePosition of this site
   */
   unsigned long _indexInsideMolecule;

   /**
   * The typeId is the typeId of this individual Site (in contrast to the typeId of the entire molecule)
   */
   unsigned long _typeId;
};
}