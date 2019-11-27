/**
 * @file ParticleTypeTrait.h
 * @author seckler
 * @date 27.11.19
 */

#pragma once
namespace autopas {
template <class Particle, class ParticleCell>
class AutoPas;
namespace utils {

template <class Container>
struct ParticleTypeTrait {
  using value = typename Container::ParticleType;
};

template <class ParticleCell>
struct ParticleTypeTrait<std::vector<ParticleCell>> {
  using value = typename ParticleCell::ParticleType;
};

template <class Particle, class ParticleCell>
struct ParticleTypeTrait<autopas::AutoPas<Particle, ParticleCell>> {
  using value = Particle;
};
}  // namespace utils
}  // namespace autopas
