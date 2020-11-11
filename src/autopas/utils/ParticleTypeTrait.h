/**
 * @file ParticleTypeTrait.h
 * @author seckler
 * @date 27.11.19
 */

#pragma once
namespace autopas {
template <class Particle>
class AutoPas;
namespace utils {

/**
 * ParticleTypeTrait class.
 * This provides the particle type for any container given.
 * The default version is supposed to take a container (LinkedCells, DirectSum, ...)
 * @tparam Container
 */
template <class Container>
struct ParticleTypeTrait {
  /**
   * The Particle Type.
   */
  using value = typename Container::ParticleType;
};

/**
 * Specialization for vectors of ParticleCells.
 * @tparam ParticleCell
 */
template <class ParticleCell>
struct ParticleTypeTrait<std::vector<ParticleCell>> {
  /**
   * The Particle Type.
   */
  using value = typename ParticleCell::ParticleType;
};

/**
 * Specialization for the AutoPas class.
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle>
struct ParticleTypeTrait<autopas::AutoPas<Particle>> {
  /**
   * The Particle Type.
   */
  using value = Particle;
};
}  // namespace utils
}  // namespace autopas
