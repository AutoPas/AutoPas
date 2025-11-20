#pragma once
#include "autopas/baseFunctors/KokkosFunctor.h"

namespace array_utils {
template <class ScalarType, int N>
struct array_type {
  std::array<ScalarType, N> the_array;

  KOKKOS_INLINE_FUNCTION
  array_type() {
    for (int i = 0; i < N; i++) {
      the_array[i] = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  array_type(const array_type &rhs) {
    for (int i = 0; i < N; i++) {
      the_array[i] = rhs.the_array[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  ScalarType &operator[](int i) { return the_array[i]; }

  KOKKOS_INLINE_FUNCTION
  const ScalarType &operator[](int i) const { return the_array[i]; }

  KOKKOS_INLINE_FUNCTION
  array_type &operator+=(const array_type &src) {
    for (int i = 0; i < N; i++) {
      the_array[i] += src.the_array[i];
    }
    return *this;
  }
};

typedef array_type<double, 3> Vector3;
}  // namespace array_utils

namespace Kokkos {
template <>
struct reduction_identity<array_utils::Vector3> {
  KOKKOS_FORCEINLINE_FUNCTION static array_utils::Vector3 sum() { return {}; }
};


}  // namespace Kokkos



namespace mdLib {

template <typename Particle_T>
class LJFunctorKokkos final : public autopas::KokkosFunctor<Particle_T, LJFunctorKokkos<Particle_T>> {


  public:
  template <typename Space>
  using SoAArraysType = autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>::template SoAArraysType<Space>;

  explicit LJFunctorKokkos(double cutoff) : autopas::KokkosFunctor<Particle_T, LJFunctorKokkos>(cutoff) {}

  bool allowsNewton3() { return true; }
  bool allowsNonNewton3() { return false; }
  bool isRelevantForTuning() { return true; }
  std::string getName() { return "Lennard-Jones Kokkos"; }

  template <typename Space>
  KOKKOS_FUNCTION static void KokkosSoAFunctor(SoAArraysType<Space> _soa, size_t index) {
    std::swap(_soa.force, _soa.old_force);

    const auto cutoff_sq = _cutoff * _cutoff;

    Kokkos::parallel_for(
        "clear force", Kokkos::RangePolicy<autopas::utils::kokkos::DeviceSpace>(0, _soa.size()),
        KOKKOS_LAMBDA(const int64_t i) {
          _soa.force(i, X_AXIS) = 0;
          _soa.force(i, Y_AXIS) = 0;
          _soa.force(i, Z_AXIS) = 0;
        });

    typedef Kokkos::View<double *[3], typename autopas::utils::kokkos::DeviceSpace::scratch_memory_space,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        block;

    auto tp = Kokkos::TeamPolicy<autopas::utils::kokkos::DeviceSpace>(_num_pairs, Kokkos::AUTO)
                  .set_scratch_size(0, Kokkos::PerTeam(block::shmem_size(BLOCK_SIZE)));

    Kokkos::parallel_for(
        "update force", tp,
        KOKKOS_CLASS_LAMBDA(const typename Kokkos::TeamPolicy<utils::kokkos::DeviceSpace>::member_type &team) {
          const auto pair_id = team.league_rank();
          const auto b1 = _pairs.d_view(pair_id, 0);
          const auto b2 = _pairs.d_view(pair_id, 1);

          auto b1_start = b1 * BLOCK_SIZE;
          auto b1_end = Kokkos::min((b1 + 1) * BLOCK_SIZE, _soa.size());

          auto b2_start = b2 * BLOCK_SIZE;
          auto b2_end = Kokkos::min((b2 + 1) * BLOCK_SIZE, _soa.size());
          ;

          block block1(team.team_scratch(0), BLOCK_SIZE);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, b1_start, b1_end), [&](const int64_t j) {
            block1(j - b1_start, X_AXIS) = _soa.position(j, X_AXIS);
            block1(j - b1_start, Y_AXIS) = _soa.position(j, Y_AXIS);
            block1(j - b1_start, Z_AXIS) = _soa.position(j, Z_AXIS);
          });
          team.team_barrier();
          if (b1 == b2) {
            // TODO optimized loop
            for (int64_t i = b1_start; i < b1_end; ++i) {
              auto type1 = _soa.typeId(i);
              array_utils::Vector3 accumulated_force;
              Kokkos::parallel_reduce(
                  Kokkos::TeamThreadRange(team, b1_start, b1_end),
                  [&](const int64_t j, array_utils::Vector3 &local_sum) {
                    if (i != j) {
                      const auto type2 = _soa.typeId(j);
                      const auto force = f(block1(i - b1_start, X_AXIS), block1(i - b1_start, Y_AXIS),
                                           block1(i - b1_start, Z_AXIS), block1(j - b1_start, X_AXIS),
                                           block1(j - b1_start, Y_AXIS), block1(j - b1_start, Z_AXIS), cutoff_sq);

                      local_sum[0] += force[0];
                      local_sum[1] += force[1];
                      local_sum[2] += force[2];
                    }
                  },
                  Kokkos::Sum<array_utils::Vector3>(accumulated_force)

              );

              Kokkos::single(Kokkos::PerTeam(team), [&]() {
                Kokkos::atomic_add(&_soa.force(i, X_AXIS), accumulated_force[0]);
                Kokkos::atomic_add(&_soa.force(i, Y_AXIS), accumulated_force[1]);
                Kokkos::atomic_add(&_soa.force(i, Z_AXIS), accumulated_force[2]);
              });
            }
          } else {
            for (int64_t i = b1_start; i < b1_end; ++i) {
              auto mask = _pairs.d_view(pair_id, 2);
              if (mask & (1L << (i - b1_start))) {
                auto type1 = _soa.typeId(i);
                array_utils::Vector3 accumulated_force;
                Kokkos::parallel_reduce(
                    Kokkos::TeamThreadRange(team, b2_start, b2_end),
                    [&](const int64_t j, array_utils::Vector3 &local_sum) {
                      auto type2 = _soa.typeId(j);
                      auto pos_x = block1(i - b1_start, X_AXIS);
                      auto pos_y = block1(i - b1_start, Y_AXIS);
                      auto pos_z = block1(i - b1_start, Z_AXIS);

                      const auto force = f(pos_x, pos_y, pos_z, _soa.position(j, X_AXIS), _soa.position(j, Y_AXIS),
                                           _soa.position(j, Z_AXIS), cutoff_sq);
                      const double Fx = force[0];
                      const double Fy = force[1];
                      const double Fz = force[2];
                      // TODO accumulate locally then update global
                      Kokkos::atomic_add(&_soa.force(j, X_AXIS), -Fx);
                      Kokkos::atomic_add(&_soa.force(j, Y_AXIS), -Fy);
                      Kokkos::atomic_add(&_soa.force(j, Z_AXIS), -Fz);

                      local_sum[0] += Fx;
                      local_sum[1] += Fy;
                      local_sum[2] += Fz;
                    },
                    Kokkos::Sum<array_utils::Vector3>(accumulated_force));

                Kokkos::single(Kokkos::PerTeam(team), [&]() {
                  Kokkos::atomic_add(&_soa.force(i, X_AXIS), accumulated_force[0]);
                  Kokkos::atomic_add(&_soa.force(i, Y_AXIS), accumulated_force[1]);
                  Kokkos::atomic_add(&_soa.force(i, Z_AXIS), accumulated_force[2]);
                });
              }
            }
          }
        });
  };
};

}  // namespace mdLib