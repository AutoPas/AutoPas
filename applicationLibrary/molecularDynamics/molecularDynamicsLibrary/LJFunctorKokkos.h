#pragma once
#include "autopas/baseFunctors/KokkosFunctor.h"
#include "autopas/utils/kokkos/ArrayUtils.h"

namespace mdLib {

template <typename Particle_T>
class LJFunctorKokkos final : public autopas::KokkosFunctor<Particle_T, LJFunctorKokkos<Particle_T>> {
 public:
  using SoAArraysType = Particle_T::SoAArraysType;

  explicit LJFunctorKokkos(double cutoff, const double epsilon, const double sigma)
      : autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>(cutoff),
        _cutoff_sq(cutoff * cutoff),
        _epsilon(epsilon),
        _sigma(sigma) {}

  bool allowsNewton3() { return true; }
  bool allowsNonNewton3() { return false; }
  bool isRelevantForTuning() { return true; }
  std::string getName() { return "Lennard-Jones Kokkos"; }

  KOKKOS_FUNCTION
  std::array<double, 3> KokkosLJ(const double self_posx, const double self_posy, const double self_posz,
                                 const double other_posx, const double other_posy, const double other_posz) const {
    const double epsilon24 = 24.0 * _epsilon;
    const double sigmaSquared = _sigma * _sigma;

    std::array<double, 3> force_vector{};
    const auto dx = self_posx - other_posx;
    const auto dy = self_posy - other_posy;
    const auto dz = self_posz - other_posz;

    const auto dr2 = dx * dx + dy * dy + dz * dz;

    double invdr2 = 1. / dr2;
    double lj2 = sigmaSquared * invdr2;
    double lj6 = lj2 * lj2 * lj2;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;

    double force_magnitude_full = epsilon24 * (lj12 + lj12m6) * invdr2;

    // Apply cutoff mask and calculate force vector
    const auto mask = (dr2 <= _cutoff_sq);
    const auto force_magnitude = force_magnitude_full * mask;

    force_vector[0] = force_magnitude * dx;
    force_vector[1] = force_magnitude * dy;
    force_vector[2] = force_magnitude * dz;

    return force_vector;
  }

  KOKKOS_FUNCTION void KokkosSoAFunctor(auto &team, auto &_soa, auto &block, size_t b1_start, size_t b1_end) const {
    for (uint64_t i = b1_start; i < b1_end; ++i) {
      autopas::utils::kokkos::ArrayUtils::Vector3 accumulated_force;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, b1_start, b1_end),
          [&](const int64_t j, autopas::utils::kokkos::ArrayUtils::Vector3 &local_sum) {
            if (i != j) {
              auto force = KokkosLJ(block(i - b1_start, 0), block(i - b1_start, 1), block(i - b1_start, 2),
                                    block(j - b1_start, 0), block(j - b1_start, 1), block(j - b1_start, 2));

              local_sum[0] += force[0];
              local_sum[1] += force[1];
              local_sum[2] += force[2];
            }
          },
          Kokkos::Sum<autopas::utils::kokkos::ArrayUtils::Vector3>(accumulated_force));

      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceX>().d_view(i), accumulated_force[0]);
        Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceY>().d_view(i), accumulated_force[1]);
        Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceZ>().d_view(i), accumulated_force[2]);
      });
    }
  };

  KOKKOS_FUNCTION void KokkosSoAFunctorPairwise(auto &team, auto &_soa, auto &block1, size_t i, size_t b1_start,
                                                size_t b2_start, size_t b2_end) const {
    autopas::utils::kokkos::ArrayUtils::Vector3 accumulated_force;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, b2_start, b2_end),
        [&](const int64_t j, autopas::utils::kokkos::ArrayUtils::Vector3 &local_sum) {
          auto pos_x = block1(i - b1_start, 0);
          auto pos_y = block1(i - b1_start, 1);
          auto pos_z = block1(i - b1_start, 2);

          const auto force =
              KokkosLJ(pos_x, pos_y, pos_z, _soa.template get<Particle_T::AttributeNames::posX>().d_view(j),
                       _soa.template get<Particle_T::AttributeNames::posY>().d_view(j),
                       _soa.template get<Particle_T::AttributeNames::posZ>().d_view(j));
          const double Fx = force[0];
          const double Fy = force[1];
          const double Fz = force[2];

          Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceX>().d_view(j), -Fx);
          Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceY>().d_view(j), -Fy);
          Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceZ>().d_view(j), -Fz);

          local_sum[0] += Fx;
          local_sum[1] += Fy;
          local_sum[2] += Fz;
        },
        Kokkos::Sum<autopas::utils::kokkos::ArrayUtils::Vector3>(accumulated_force));

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceX>().d_view(i), accumulated_force[0]);
      Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceY>().d_view(i), accumulated_force[1]);
      Kokkos::atomic_add(&_soa.template get<Particle_T::AttributeNames::forceZ>().d_view(i), accumulated_force[2]);
    });
  };

  KOKKOS_FUNCTION
  void load_particle(auto &soa, Particle_T &particle, size_t index) const {
    Particle_T *p = &particle;
    soa.template get<Particle_T::AttributeNames::ptr>().h_view(index) = reinterpret_cast<uintptr_t>(p);

    soa.template get<Particle_T::AttributeNames::forceZ>().h_view(index) =
        particle.template get<Particle_T::AttributeNames::forceX>();

    soa.template get<Particle_T::AttributeNames::forceY>().h_view(index) =
        particle.template get<Particle_T::AttributeNames::forceY>();

    soa.template get<Particle_T::AttributeNames::forceZ>().h_view(index) =
        particle.template get<Particle_T::AttributeNames::forceZ>();

    soa.template get<Particle_T::AttributeNames::posX>().h_view(index) =
        particle.template get<Particle_T::AttributeNames::posX>();
    soa.template get<Particle_T::AttributeNames::posY>().h_view(index) =
        particle.template get<Particle_T::AttributeNames::posY>();
    soa.template get<Particle_T::AttributeNames::posZ>().h_view(index) =
        particle.template get<Particle_T::AttributeNames::posZ>();
  }

  KOKKOS_FUNCTION
  void store_particle(auto &soa, size_t index) const {
    auto particle = reinterpret_cast<Particle_T *>(soa.template get<Particle_T::AttributeNames::ptr>().h_view(index));

    particle->template set<Particle_T::AttributeNames::forceX>(
        soa.template get<Particle_T::AttributeNames::forceX>().h_view(index));
    particle->template set<Particle_T::AttributeNames::forceY>(
        soa.template get<Particle_T::AttributeNames::forceY>().h_view(index));
    particle->template set<Particle_T::AttributeNames::forceZ>(
        soa.template get<Particle_T::AttributeNames::forceZ>().h_view(index));

    particle->template set<Particle_T::AttributeNames::blockId>(
      soa.template get<Particle_T::AttributeNames::blockId>().h_view(index));
  }

 private:
  const double _cutoff_sq, _epsilon, _sigma;
};

}  // namespace mdLib