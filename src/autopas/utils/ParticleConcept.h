/**
 * @file ParticleConcept.h
 * @author S. Newcome
 * @date 28.07.26
 */

#pragma once

#include <array>
#include <concepts>
#include <ostream>
#include <string>
#include <type_traits>

#include "autopas/particles/OwnershipState.h"

namespace autopas::utils {

/**
 * Concept describing everything AutoPas requires of a particle.
 *
 * This mirrors the public interface of autopas::ParticleBase, which is the recommended base class for any particle
 * used with AutoPas. Satisfying this concept without inheriting from ParticleBase is possible but not encouraged.
 *
 * @note Particles are deliberately not required to be default constructible, only copyable, to avoid imposing any
 * requirement upon the constructor.
 * @note Only the SoA attributes that every particle defines (ptr, posX, posY, posZ and ownershipState) are required
 * here. Attributes such as id or the force components are not universal, e.g. sphLib::SPHParticle defines neither.
 *
 * @warning In addition to this concept, we require a markAsDeleted() function on the particle. As this is a private
 * function, accessed only through a friend function, we cannot check it in this concept. If developing a custom
 * particle derived from ParticleBase, this function will already be included, but a purely custom particle must
 * implement it and befriend internal::markParticleAsDeleted():
 * @code
 * template <autopas::utils::ParticleType T>
 * friend void autopas::internal::markParticleAsDeleted(T &);
 * @endcode
 *
 * @tparam Particle_T Type to check.
 */
template <class Particle_T>
concept ParticleType =
    std::copyable<Particle_T> and

    requires(Particle_T p, const Particle_T cp, const std::array<double, 3> &vec, double maxDistSquared,
             std::ostream &os) {
      // Types needed to build and address the SoAs.
      typename Particle_T::SoAArraysType;
      typename Particle_T::AttributeNames;
      requires std::is_enum_v<typename Particle_T::AttributeNames>;

      // Position.
      { cp.getR() } -> std::same_as<const std::array<double, 3> &>;
      { p.setR(vec) } -> std::same_as<void>;
      { p.addR(vec) } -> std::same_as<void>;
      { p.setRDistanceCheck(vec, maxDistSquared) } -> std::same_as<bool>;
      { p.addRDistanceCheck(vec, maxDistSquared) } -> std::same_as<bool>;

      // Velocity.
      { cp.getV() } -> std::same_as<const std::array<double, 3> &>;
      { p.setV(vec) } -> std::same_as<void>;
      { p.addV(vec) } -> std::same_as<void>;

      // Force.
      { cp.getF() } -> std::same_as<const std::array<double, 3> &>;
      { p.setF(vec) } -> std::same_as<void>;
      { p.addF(vec) } -> std::same_as<void>;
      { p.subF(vec) } -> std::same_as<void>;

      // Id. Phrased so that the id type does not have to be named.
      { cp.getID() } -> std::integral;
      { p.setID(cp.getID()) } -> std::same_as<void>;

      // Ownership.
      { cp.isOwned() } -> std::same_as<bool>;
      { cp.isHalo() } -> std::same_as<bool>;
      { cp.isDummy() } -> std::same_as<bool>;
      { cp.getOwnershipState() } -> std::same_as<OwnershipState>;
      { p.setOwnershipState(OwnershipState::owned) } -> std::same_as<void>;

      // String representation.
      { cp.toString() } -> std::same_as<std::string>;
      { os << cp } -> std::same_as<std::ostream &>;

      // SoA access.
      { p.template get<Particle_T::AttributeNames::ptr>() } -> std::same_as<Particle_T *>;
      { cp.template get<Particle_T::AttributeNames::posX>() } -> std::floating_point;
      { cp.template get<Particle_T::AttributeNames::posY>() } -> std::floating_point;
      { cp.template get<Particle_T::AttributeNames::posZ>() } -> std::floating_point;
      { cp.template get<Particle_T::AttributeNames::ownershipState>() } -> std::same_as<OwnershipState>;
      { p.template set<Particle_T::AttributeNames::posX>(double{}) } -> std::same_as<void>;
      { p.template set<Particle_T::AttributeNames::posY>(double{}) } -> std::same_as<void>;
      { p.template set<Particle_T::AttributeNames::posZ>(double{}) } -> std::same_as<void>;
      { p.template set<Particle_T::AttributeNames::ownershipState>(OwnershipState::owned) } -> std::same_as<void>;
    }

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    and requires(Particle_T p, const Particle_T cp, const std::array<double, 3> &vec) {
      { cp.getRAtRebuild() } -> std::same_as<const std::array<double, 3> &>;
      { p.setRAtRebuild(vec) } -> std::same_as<void>;
      { p.resetRAtRebuild() } -> std::same_as<void>;
      // Returned by value as a const array, so this can not be matched with std::same_as.
      { cp.calculateDisplacementSinceRebuild() } -> std::convertible_to<std::array<double, 3>>;
    }
#endif
    ;

}  // namespace autopas::utils
