/**
 * @file Object.h
 * @author N. Fottner
 * @date 1/8/19
 */
#pragma once

#include <array>
#include <iomanip>
#include <iosfwd>
#include <vector>

#include "AbsoluteMultiSiteMoleculeInitializer.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/Quaternion.h"
#include "src/MoleculeContainer.h"
#include "src/TypeDefinitions.h"

/**
 * Base class for describing objects made of particles.
 */
class Object {
 public:
  /**
   * Constructor that should be used by inheriting types.
   * @param velocity
   * @param typeId
   */
  Object(const std::array<double, 3> &velocity, unsigned long typeId) : _velocity(velocity), _typeId(typeId) {}

  virtual ~Object() = default;

  /**
   * Generate the object in the given AutoPas container.
   * @param particles The container to which the new particles will be appended to.
   */
#if (MD_FLEXIBLE_MODE != MULTISITE) or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
  virtual void generate(std::vector<ParticleType> &particles, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl) const = 0;
#else
  virtual void generate(std::vector<ParticleType> &particles, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl,
                        MoleculeContainer& moleculeContainer) const = 0;
#endif

#if (MD_FLEXIBLE_MODE != MULTISITE) or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
  void insertMolecule(const std::array<double, 3>& position, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl,
                                    std::vector<ParticleType> &particles) {
    //the typeID set in getDummyParticle is the ID of the MOLECULE. Since in this configuration a Particle is just a single
    //site we need to overwrite this typeID with the proper siteID
    ParticleType particle = getDummyParticle(particles.size());
    particle.setR(position);
#if defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
    const ParticlePropertiesLibraryType * ppl_ptr = ppl.get();  //black magic to get const reference out of const std::shared_ptr
    particle.setQuaternion(particle.getQuaternion(), *ppl_ptr);
#endif
    particles.push_back(particle);
  }
#else
  void insertMolecule(const std::array<double, 3>& position, const std::shared_ptr<const ParticlePropertiesLibraryType> ppl,
                                    std::vector<ParticleType> &particles, MoleculeContainer& moleculeContainer) const {
    //insert molecule by inserting individual sites at their respective place and inserting the bundled molecule in moleculeContainer
    const auto moleculeID = moleculeContainer.size();
    MoleculeType molecule = getDummyMolecule(moleculeID);
    molecule.setR(position);
    moleculeContainer.add(std::move(molecule));

    //the typeID set in getDummyParticle is the ID of the MOLECULE. Since in this configuration a Particle is just a single
    //site we need to overwrite this typeID with the proper siteID
    ParticleType particle = getDummyParticle(particles.size());
    const auto siteTypes = ppl->getSiteTypes(_typeId);
    const auto unrotatedRelativeSitePositions = ppl->getSitePositions(_typeId);
    const auto rotatedRelativeSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(molecule.getQuaternion(), unrotatedRelativeSitePositions);

    //insert individual sites
    for(size_t siteIndex=0; siteIndex<siteTypes.size(); siteIndex++) {
      const auto& relSitePos = rotatedRelativeSitePositions[siteIndex];
      const auto siteType = siteTypes[siteIndex];

      particle.setR(autopas::utils::ArrayMath::add(position, relSitePos));
      particle.setTypeId(siteType); //overwrite previously stored moleculeID with the actual SiteID //when i am actually storing the siteType here then i don't have  a place to store the reference to the molecule. That's why i am using the typeID as moleculeID instead
      particle.setID(particle.getID() + 1);
      particle.setMoleculeId(moleculeID);
      particles.push_back(particle);
    }
  }
#endif

  /**
   * Create a particle that acts as blueprint for all particles to be created for the object.
   * @param particleId: Defines the id of the generated dummy particle.
   * @return a particle initialized with default values.
   */
  [[nodiscard]] ParticleType getDummyParticle(const size_t &particleId) const {
    ParticleType particle{};
    particle.setID(particleId);
    particle.setTypeId(_typeId);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    particle.setV(_velocity);
    particle.setF({0.0, 0.0, 0.0});
    particle.setOldF({0.0, 0.0, 0.0});

#if MD_FLEXIBLE_MODE == MULTISITE and not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
#if not defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
    particle.setQuaternion({1.0, 0.0, 0.0, 0.0});  // todo: add option for this to be set randomly
#else
    particle.setOnlyQuaternion({1.0, 0.0, 0.0, 0.0});
#endif
    particle.setAngularVel({0.0, 0.0, 0.0});
    particle.setTorque({0.0, 0.0, 0.0});
#endif

    return particle;
  }

#if defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH) and (MD_FLEXIBLE_MODE == MULTISITE)
  [[nodiscard]] MoleculeType getDummyMolecule(const size_t &particleId) const {
    MoleculeType molecule{};
    molecule.setID(particleId);
    molecule.setTypeId(_typeId);
    molecule.setOwnershipState(autopas::OwnershipState::owned);
    molecule.setV(_velocity);
    molecule.setF({0.0, 0.0, 0.0});
    molecule.setOldF({0.0, 0.0, 0.0});
    molecule.setQuaternion({1.0, 0.0, 0.0, 0.0});

    return molecule;
  }
#endif

  /**
   * Getter for Velocity
   * @return velocity
   */
  [[nodiscard]] const std::array<double, 3> &getVelocity() const { return _velocity; }

  /**
   * Getter for typeId of Particles in Objet
   * @return typeId
   */
  [[nodiscard]] unsigned long getTypeId() const { return _typeId; }

  /**
   * Getter for the smallest x,y,z coordinates for Object
   * @return BoxMin of Cube
   */
  [[nodiscard]] virtual std::array<double, 3> getBoxMin() const = 0;

  /**
   * Getter for the highest x,y,z coordinates for Object
   * @return BoxMax of Cube
   */
  [[nodiscard]] virtual std::array<double, 3> getBoxMax() const = 0;

  /**
   * Returns the total amount of Particles in the Object
   * @return ParticlesTotal
   */
  [[nodiscard]] virtual size_t getParticlesTotal() const = 0;

  /**
   * Getter for ParticleSpacing.
   * Objects that are not based on a grid return 0 since this is the minimal guaranteed spacing.
   * @return particleSpacing
   */
  [[nodiscard]] virtual double getParticleSpacing() const { return 0; }

  /**
   * String description string of the object.
   * @return multiline std::string
   */
  [[nodiscard]] virtual std::string to_string() const {
    std::ostringstream output;
    output << std::setw(_valueOffset) << std::left << "velocity"
           << ":  " << autopas::utils::ArrayUtils::to_string(_velocity) << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-type-id"
           << ":  " << _typeId << std::endl;
    return output.str();
  };

  /**
   * Stream operator
   * @param os
   * @param object
   * @return
   */
  friend std::ostream &operator<<(std::ostream &os, const Object &object) {
    os << object.to_string();
    return os;
  }

 protected:
  /**
   * Velocity of every particle in the object.
   */
  std::array<double, 3> _velocity;
  /**
   * Type of every particle in the object. For single-site simulations, this refers directly to the siteId. For
   * multi-site simulations, this refers to the molId.
   */
  unsigned long _typeId;
  /**
   * valueOffset of MDFlexConfig - expected indent
   */
  static constexpr size_t _valueOffset = 33 - 6;
};
