#pragma once
#include "autopas/baseFunctors/KokkosFunctor.h"
#include "autopas/utils/kokkos/ArrayUtils.h"

namespace mdLib {

template <typename Particle_T>
class LJFunctorKokkos final : public autopas::KokkosFunctor<Particle_T, LJFunctorKokkos<Particle_T>> {
 public:
  using SoAArraysType = autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>::SoAArraysType;

  explicit LJFunctorKokkos(double cutoff) : autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>(cutoff) {}

  bool allowsNewton3() { return true; }
  bool allowsNonNewton3() { return false; }
  bool isRelevantForTuning() { return true; }
  std::string getName() { return "Lennard-Jones Kokkos"; }

  KOKKOS_FUNCTION
  std::array<double, 3> KokkosLJ(const double self_posx, const double self_posy, const double self_posz,
                                 const double other_posx, const double other_posy, const double other_posz,
                                 const double cutoff_radius_sq) const {
    const double sigma_6 = sigma * sigma * sigma * sigma * sigma * sigma;
    const double sigma_12 = sigma_6 * sigma_6;

    std::array<double, 3> force_vector{};
    const auto dx = self_posx - other_posx;
    const auto dy = self_posy - other_posy;
    const auto dz = self_posz - other_posz;

    const auto r_sq = dx * dx + dy * dy + dz * dz;

    const auto r2_inv = 1.0 / r_sq;
    const auto r6_inv = r2_inv * r2_inv * r2_inv;  // 1/r^6
    const auto r12_inv = r6_inv * r6_inv;          // 1/r^12

    const auto force_magnitude_full =
        r2_inv * ((48.0 * epsilon * sigma_12) * r12_inv - (24.0 * epsilon * sigma_6) * r6_inv);

    // this avoids branching and ensures all threads
    // follow the same instruction path on GPU
    const auto mask = (r_sq <= cutoff_radius_sq);
    const auto force_magnitude = force_magnitude_full * mask;

    force_vector[0] = force_magnitude * dx;
    force_vector[1] = force_magnitude * dy;
    force_vector[2] = force_magnitude * dz;

    return force_vector;
  }

  KOKKOS_FUNCTION static void KokkosSoAFunctor(auto &team, SoAArraysType &_soa, auto &block, size_t b1_start,
                                               size_t b1_end) {
    for (uint64_t i = b1_start; i < b1_end; ++i) {
      autopas::utils::kokkos::ArrayUtils::Vector3 accumulated_force;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, b1_start, b1_end),
          [&](const int64_t j, autopas::utils::kokkos::ArrayUtils::Vector3 &local_sum) {
            if (i != j) {
              const auto force = std::array<double, 3>{0.0, 0.0, 0.0};

              KokkosLJ()(block(i - b1_start, 0), block(i - b1_start, 1), block(i - b1_start, 2), block(j - b1_start, 0),
                         block(j - b1_start, 1), block(j - b1_start, 2), getCutoff() * getCutoff());

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

  KOKKOS_FUNCTION static void KokkosSoAFunctorPairwise(auto &team, SoAArraysType &_soa, auto &block1, size_t i,
                                                       size_t b1_start, size_t b2_start, size_t b2_end) {
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
                       _soa.template get<Particle_T::AttributeNames::posZ>().d_view(j), getCutoff() * getCutoff());
          const double Fx = force[0];
          const double Fy = force[1];
          const double Fz = force[2];

          // TODO accumulate locally then update global
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
};

}  // namespace mdLib