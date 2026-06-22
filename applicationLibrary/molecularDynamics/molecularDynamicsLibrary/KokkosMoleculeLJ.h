/**
 * @file KokkosMoleculeLJ.h
 * @date 23 Jan 2026 copied from MoleculeLJ.h
 * @author Luis Gall
 */

#pragma once

#include <vector>

#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utilsKokkos/KokkosSoA.h"

namespace mdLib {

/**
 * Molecule class for the LJFunctor.
 */
class KokkosMoleculeLJ {
 public:

  using ParticleSoAFloatPrecision = float;

  KOKKOS_INLINE_FUNCTION KokkosMoleculeLJ() : KokkosMoleculeLJ(Kokkos::Array<ParticleSoAFloatPrecision, 3>{0., 0., 0.}, Kokkos::Array<ParticleSoAFloatPrecision, 3>{0., 0., 0.}, 0, 0) {}

  KOKKOS_INLINE_FUNCTION KokkosMoleculeLJ(const Kokkos::Array<ParticleSoAFloatPrecision, 3> &pos, const Kokkos::Array<ParticleSoAFloatPrecision, 3> &v, unsigned long moleculeId,
                         unsigned long typeId) :
                            _r{pos}, _v{v}, _id{moleculeId}, _typeId {typeId}, _mass{0}, _ownershipState{autopas::OwnershipState::owned}
                            , _f{0., 0., 0.}
                          #ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
                            , _rAtRebuild{0., 0., 0.}
                          #endif
                            {}

  KokkosMoleculeLJ(const std::array<ParticleSoAFloatPrecision, 3> &pos, const std::array<ParticleSoAFloatPrecision, 3> &v, unsigned long moleculeId,
                         unsigned long typeId)
    : _r{pos.at(0), pos.at(1), pos.at(2)}
  , _v{v.at(0), v.at(1), v.at(2)}
  , _f{0., 0., 0.}
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  , _rAtRebuild{0., 0., 0.}
#endif
  , _id(moleculeId)
  , _mass(0)
  , _ownershipState{autopas::OwnershipState::owned}
  , _typeId(typeId) {}

  ~KokkosMoleculeLJ() = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : size_t {
    ptr,
    id,
    posX,
    posY,
    posZ,
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    rebuildX,
    rebuildY,
    rebuildZ,
#endif
    velocityX,
    velocityY,
    velocityZ,
    forceX,
    forceY,
    forceZ,
    oldForceX,
    oldForceY,
    oldForceZ,
    mass,
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
  using SoAArraysType = autopas::utils::SoAType<
      KokkosMoleculeLJ *, size_t /*id*/, ParticleSoAFloatPrecision /*x*/, ParticleSoAFloatPrecision /*y*/,
      ParticleSoAFloatPrecision /*z*/,
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
      ParticleSoAFloatPrecision /*rebuildX*/, ParticleSoAFloatPrecision /*rebuildY*/,
      ParticleSoAFloatPrecision /*rebuildZ*/,
#endif
      ParticleSoAFloatPrecision /*vx*/, ParticleSoAFloatPrecision /*vy*/, ParticleSoAFloatPrecision /*vz*/,
      ParticleSoAFloatPrecision /*fx*/, ParticleSoAFloatPrecision /*fy*/, ParticleSoAFloatPrecision /*fz*/,
      ParticleSoAFloatPrecision /*oldFx*/, ParticleSoAFloatPrecision /*oldFy*/, ParticleSoAFloatPrecision /*oldFz*/,
      ParticleSoAFloatPrecision /*mass*/, size_t /*typeid*/, autopas::OwnershipState /*ownershipState*/>::Type;

  using KokkosSoAArraysType =
      autopas::utilsKokkos::KokkosSoA<size_t * /*id*/, ParticleSoAFloatPrecision * /*x*/, ParticleSoAFloatPrecision * /*y*/,
                                ParticleSoAFloatPrecision * /*z*/,
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
                                ParticleSoAFloatPrecision * /*rebuildX*/, ParticleSoAFloatPrecision * /*rebuildY*/,
                                ParticleSoAFloatPrecision * /*rebuildZ*/,
#endif
                                ParticleSoAFloatPrecision * /*vx*/, ParticleSoAFloatPrecision * /*vy*/,
                                ParticleSoAFloatPrecision * /*vz*/, ParticleSoAFloatPrecision * /*fx*/,
                                ParticleSoAFloatPrecision * /*fy*/, ParticleSoAFloatPrecision * /*fz*/,
                                ParticleSoAFloatPrecision * /*oldFx*/, ParticleSoAFloatPrecision * /*oldFy*/,
                                ParticleSoAFloatPrecision * /*oldFz*/, ParticleSoAFloatPrecision * /*mass*/,
                                size_t * /*typeid*/, autopas::OwnershipState * /*ownershipState*/>;

  template <AttributeNames attribute>
  constexpr decltype(auto) operator()() {
    if constexpr (attribute == ptr) {
      return this;
    } else {
      return getRef<attribute>();
    }
  }

  template <AttributeNames attribute>
  constexpr auto operator()() const {
    return get<attribute>();
  }

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
  template <AttributeNames attribute, std::enable_if_t<attribute != ptr, bool> = true>
  constexpr std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    if constexpr (attribute == id) {
      return _id;
    } else if constexpr (attribute == posX) {
      return _r[0];
    } else if constexpr (attribute == posY) {
      return _r[1];
    } else if constexpr (attribute == posZ) {
      return _r[2];
    }
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    else if constexpr (attribute == rebuildX) {
      return _rAtRebuild[0];
    } else if constexpr (attribute == rebuildY) {
      return _rAtRebuild[1];
    } else if constexpr (attribute == rebuildZ) {
      return _rAtRebuild[2];
    }
#endif
    else if constexpr (attribute == velocityX) {
      return _v[0];
    } else if constexpr (attribute == velocityY) {
      return _v[1];
    } else if constexpr (attribute == velocityZ) {
      return _v[2];
    } else if constexpr (attribute == forceX) {
      return _f[0];
    } else if constexpr (attribute == forceY) {
      return _f[1];
    } else if constexpr (attribute == forceZ) {
      return _f[2];
    } else if constexpr (attribute == oldForceX) {
      return _oldF[0];
    } else if constexpr (attribute == oldForceY) {
      return _oldF[1];
    } else if constexpr (attribute == oldForceZ) {
      return _oldF[2];
    } else if constexpr (attribute == mass) {
      return _mass;
    } else if constexpr (attribute == typeId) {
      return _typeId;
    } else if constexpr (attribute == ownershipState) {
      return _ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  template <AttributeNames attribute, std::enable_if_t<attribute != ptr, bool> = true>
  constexpr decltype(auto) getRef() {
    if constexpr (attribute == id) {
      return (_id);
    } else if constexpr (attribute == posX) {
      return (_r[0]);
    } else if constexpr (attribute == posY) {
      return (_r[1]);
    } else if constexpr (attribute == posZ) {
      return (_r[2]);
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    } else if constexpr (attribute == rebuildX) {
      return (_rAtRebuild[0]);
    } else if constexpr (attribute == rebuildY) {
      return (_rAtRebuild[1]);
    } else if constexpr (attribute == rebuildZ) {
      return (_rAtRebuild[2]);
#endif
    } else if constexpr (attribute == velocityX) {
      return (_v[0]);
    } else if constexpr (attribute == velocityY) {
      return (_v[1]);
    } else if constexpr (attribute == velocityZ) {
      return (_v[2]);
    } else if constexpr (attribute == forceX) {
      return (_f[0]);
    } else if constexpr (attribute == forceY) {
      return (_f[1]);
    } else if constexpr (attribute == forceZ) {
      return (_f[2]);
    } else if constexpr (attribute == oldForceX) {
      return (_oldF[0]);
    } else if constexpr (attribute == oldForceY) {
      return (_oldF[1]);
    } else if constexpr (attribute == oldForceZ) {
      return (_oldF[2]);
    } else if constexpr (attribute == mass) {
      return (_mass);
    } else if constexpr (attribute == typeId) {
      return (_typeId);
    } else if constexpr (attribute == ownershipState) {
      return (_ownershipState);
    } else {
      autopas::utils::ExceptionHandler::exception("KokkosMoleculeLJ::getRef() unknown attribute {}", attribute);
    }
  }

  template <AttributeNames attribute, std::enable_if_t<attribute != ptr, bool> = true>
  constexpr const std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == id) {
      return _id;
    } else if constexpr (attribute == posX) {
      return _r[0];
    } else if constexpr (attribute == posY) {
      return _r[1];
    } else if constexpr (attribute == posZ) {
      return _r[2];
    }
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    else if constexpr (attribute == rebuildX) {
      return _rAtRebuild[0];
    } else if constexpr (attribute == rebuildY) {
      return _rAtRebuild[1];
    } else if constexpr (attribute == rebuildZ) {
      return _rAtRebuild[2];
    }
#endif
    else if constexpr (attribute == velocityX) {
      return _v[0];
    } else if constexpr (attribute == velocityY) {
      return _v[1];
    } else if constexpr (attribute == velocityZ) {
      return _v[2];
    } else if constexpr (attribute == forceX) {
      return _f[0];
    } else if constexpr (attribute == forceY) {
      return _f[1];
    } else if constexpr (attribute == forceZ) {
      return _f[2];
    } else if constexpr (attribute == oldForceX) {
      return _oldF[0];
    } else if constexpr (attribute == oldForceY) {
      return _oldF[1];
    } else if constexpr (attribute == oldForceZ) {
      return _oldF[2];
    } else if constexpr (attribute == mass) {
      return _mass;
    } else if constexpr (attribute == typeId) {
      return _typeId;
    } else if constexpr (attribute == ownershipState) {
      return _ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::get() unknown attribute {}", attribute);
    }
  }

  template <AttributeNames attribute>
  constexpr void set(std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == id) {
      _id = value;
    } else if constexpr (attribute == posX) {
      _r[0] = value;
    } else if constexpr (attribute == posY) {
      _r[1] = value;
    } else if constexpr (attribute == posZ) {
      _r[2] = value;
    }
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    else if constexpr (attribute == rebuildX) {
      _rAtRebuild[0] = value;
    } else if constexpr (attribute == rebuildY) {
      _rAtRebuild[1] = value;
    } else if constexpr (attribute == rebuildZ) {
      _rAtRebuild[2] = value;
    }
#endif
    else if constexpr (attribute == velocityX) {
      _v[0] = value;
    } else if constexpr (attribute == velocityY) {
      _v[1] = value;
    } else if constexpr (attribute == velocityZ) {
      _v[2] = value;
    } else if constexpr (attribute == forceX) {
      _f[0] = value;
    } else if constexpr (attribute == forceY) {
      _f[1] = value;
    } else if constexpr (attribute == forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == oldForceX) {
      _oldF[0] = value;
    } else if constexpr (attribute == oldForceY) {
      _oldF[1] = value;
    } else if constexpr (attribute == oldForceZ) {
      _oldF[2] = value;
    } else if constexpr (attribute == mass) {
      _mass = value;
    } else if constexpr (attribute == typeId) {
      _typeId = value;
    } else if constexpr (attribute == ownershipState) {
      _ownershipState = static_cast<autopas::OwnershipState>(value);
    } else {
      autopas::utils::ExceptionHandler::exception("MoleculeLJ::set() unknown attribute {}", attribute);
    }
  }

  void setID(unsigned long id) {
    _id = id;
  };

  void setOwnershipState(autopas::OwnershipState value) {
    _ownershipState = value;
  }

  void setV(const std::array<ParticleSoAFloatPrecision, 3>& v) {
    const Kokkos::Array temp {v.at(0), v.at(1), v.at(2)};
    _v = temp;
  }

  void setMass(ParticleSoAFloatPrecision mass) {
    _mass = mass;
  }

  void setF(const std::array<ParticleSoAFloatPrecision, 3>& f) {
    const Kokkos::Array temp {f.at(0), f.at(1), f.at(2)};
    _f = temp;
  }

  void setR(const std::array<ParticleSoAFloatPrecision, 3>& r) {
    const Kokkos::Array temp {r.at(0), r.at(1), r.at(2)};
    _r = temp;
  }

  void setOldF(const std::array<ParticleSoAFloatPrecision, 3>& oldForce) {
    const Kokkos::Array temp {oldForce.at(0), oldForce.at(1), oldForce.at(2)};
    _oldF = temp;
  }

  void markAsDeleted () {
    setOwnershipState(autopas::OwnershipState::dummy);
  }

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  void resetRAtRebuild() {
    setRAtRebuild(getR());
  }

  std::array<ParticleSoAFloatPrecision, 3> getRAtRebuild() const {
    return {_rAtRebuild[0], _rAtRebuild[1], _rAtRebuild[2]};
  }

  void setRAtRebuild(const std::array<ParticleSoAFloatPrecision, 3>& r) {
    _rAtRebuild[0] = r[0];
    _rAtRebuild[1] = r[1];
    _rAtRebuild[2] = r[2];
  }

  std::array<ParticleSoAFloatPrecision, 3> calculateDisplacementSinceRebuild () const {
      return {
        _rAtRebuild[0] - _r[0], _rAtRebuild[1] - _r[1], _rAtRebuild[2] - _r[2]
      };
   }
#endif

  /**
 * Defines whether the particle is owned by the current AutoPas object (aka (MPI-)process)
 * @return true if the particle is owned by the current AutoPas object, false otherwise
 */
  [[nodiscard]] bool isOwned() const { return _ownershipState == autopas::OwnershipState::owned; }

  /**
   * Defines whether the particle is a halo particle, i.e., not owned by the current AutoPas object (aka (MPI-)process)
   * @return true if the particle is not owned by the current AutoPas object, false otherwise.
   * @note when a
   */
  [[nodiscard]] bool isHalo() const { return _ownershipState == autopas::OwnershipState::halo; }

  /**
   * Returns whether the particle is a dummy particle.
   * @return true if the particle is a dummy.
   */
  [[nodiscard]] bool isDummy() const { return _ownershipState == autopas::OwnershipState::dummy; }

  unsigned long getID() const { return _id; }

  std::array<ParticleSoAFloatPrecision, 3> getV() const { return {_v[0], _v[1], _v[2] }; }

  std::array<ParticleSoAFloatPrecision, 3> getR() const { return {_r[0], _r[1], _r[2]} ; }

  std::array<ParticleSoAFloatPrecision, 3> getF() const { return {_f[0], _f[1], _f[2] }; }

  std::array<ParticleSoAFloatPrecision, 3> getOldF() const { return {_oldF[0], _oldF[1], _oldF[2] }; }

  autopas::OwnershipState getOwnershipState() const { return _ownershipState; }

  void addV(const std::array<ParticleSoAFloatPrecision, 3>& increment) {
    _v[0] += increment[0];
    _v[1] += increment[1];
    _v[2] += increment[2];
  }

  void addF(const std::array<ParticleSoAFloatPrecision, 3>& increment) {
    _f[0] += increment[0];
    _f[1] += increment[1];
    _f[2] += increment[2];
  }

  void addR(const std::array<ParticleSoAFloatPrecision, 3>& increment) {
    _r[0] += increment[0];
    _r[1] += increment[1];
    _r[2] += increment[2];
  }

  /**
   * Get TypeId.
   * @return
   */
  [[nodiscard]] size_t getTypeId() const;

  /**
   * Set the type id of the Molecule.
   * @param typeId
   */
  void setTypeId(size_t typeId);

  /**
   * Creates a string containing all data of the particle.
   * @return String representation.
   */
  [[nodiscard]] std::string toString() const;

  friend std::ostream &operator<<(std::ostream &os, const KokkosMoleculeLJ &particle) {
    using autopas::utils::ArrayUtils::operator<<;
    os << "Particle"
       << "\nID      : " << particle._id << "\nOwnershipState : " << particle._ownershipState;
    // clang-format on
    return os;
  }

 protected:
  /**
   * Molecule type id. In single-site simulations, this is used as a siteId to look up site attributes in the particle
   * properties library.
   *
   * In multi-site simulations, where a multi-site molecule class inheriting from this class is used, typeId is used as
   * a molId to look up molecular attributes (including siteIds of the sites).
   */
  size_t _typeId = 0;

  Kokkos::Array<ParticleSoAFloatPrecision, 3> _r;

  Kokkos::Array<ParticleSoAFloatPrecision, 3> _v;

  Kokkos::Array<ParticleSoAFloatPrecision, 3> _f;

  ParticleSoAFloatPrecision _mass;

  autopas::OwnershipState _ownershipState;

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  Kokkos::Array<ParticleSoAFloatPrecision, 3> _rAtRebuild;
#endif

  unsigned long _id = 0;

  /**
   * Old Force of the particle experiences as 3D vector.
   */
  Kokkos::Array<ParticleSoAFloatPrecision, 3> _oldF = {0., 0., 0.};
};

}  // namespace mdLib
