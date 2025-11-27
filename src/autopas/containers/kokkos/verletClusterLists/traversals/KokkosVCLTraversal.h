#pragma once

#include <filesystem>

#include "KokkosTraversalInterface.h"
#include "Kokkos_Sort.hpp"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/utils/kokkos/ArrayUtils.h"
#include "autopas/utils/kokkos/KokkosSoA.h"

constexpr auto normalize(const double coord, const double cell_size) -> int64_t {
  const double scaled_distance = std::max(0.0, coord) / cell_size;
  return static_cast<int64_t>(scaled_distance);
}

KOKKOS_INLINE_FUNCTION
constexpr void axisTranspose(int X[3], int b) {
  const int M = 1u << (b - 1);
  int P, Q, t;

  for (Q = M; Q > 1u; Q >>= 1) {
    P = Q - 1u;
    for (int i = 0; i < 3; ++i) {
      if (X[i] & Q)
        X[0] ^= P;
      else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }
  X[1] ^= X[0];
  X[2] ^= X[1];

  t = 0;
  for (Q = M; Q > 1u; Q >>= 1) {
    if (X[2] & Q) t ^= (Q - 1u);
  }
  X[0] ^= t;
  X[1] ^= t;
  X[2] ^= t;
}

KOKKOS_INLINE_FUNCTION
constexpr auto hilbert_index_3d(int x, int y, int z) -> unsigned long {
  // return morton_code(x,y,z);

  constexpr int BITS = 21;

  int X[3] = {x, y, z};
  axisTranspose(X, BITS);

  unsigned long index = 0;
  for (int i = 0; i < BITS; ++i) {
    unsigned long bx = (X[0] >> i) & 1u;
    unsigned long by = (X[1] >> i) & 1u;
    unsigned long bz = (X[2] >> i) & 1u;
    index |= (bx << (3 * i + 2)) | (by << (3 * i + 1)) | (bz << (3 * i + 0));
  }
  return index;
}

namespace autopas::containers::kokkos::traversal {

template <class Particle_T, typename Functor_T>
class KokkosVCLTraversal : public autopas::TraversalInterface, public KokkosTraversalInterface<Particle_T> {
 public:
  KokkosVCLTraversal(Functor_T *f, DataLayoutOption dataLayout, bool newton3, double cutoff,
                     std::array<double, 3> cellSize, double skinSize = 0.0)
      : TraversalInterface(dataLayout, newton3),
        _functor( *f ),
        _cutoff(cutoff),
        _cellSize(cellSize),
        _skinSize(skinSize) {}

  ~KokkosVCLTraversal() override = default;

  [[nodiscard]] autopas::TraversalOption getTraversalType() const override { return autopas::TraversalOption::kk_vcl; }

  [[nodiscard]] bool isApplicable() const override { return true; }
  void initTraversal() override {
    // TODO make dynamic
    clusterParticles();
    updateBoundingBox();
    updatePairs();
  }
  void traverseParticles() override { calculateForce(); }
  void endTraversal() override {
    // TODO copy back
  }

  void rebuild(std::vector<Particle_T> &particles) override {
    // TODO add particles
    clusterParticles();
    updateBoundingBox();
    updatePairs();
  }

 private:
  constexpr static int X_AXIS = 0;
  constexpr static int Y_AXIS = 1;
  constexpr static int Z_AXIS = 2;
  constexpr static int AXIS_MIN = 0;
  constexpr static int AXIS_MAX = 1;

  constexpr static int BLOCK_SIZE = 32;


  void calculateForce() {
    Kokkos::parallel_for(
        "clear force", Kokkos::RangePolicy<autopas::utils::kokkos::DeviceSpace>(0, _soa.size()),
        KOKKOS_CLASS_LAMBDA(const int64_t i) {
          _soa.template get<Particle_T::AttributeNames::forceX>().d_view(i) = 0;
          _soa.template get<Particle_T::AttributeNames::forceY>().d_view(i) = 0;
          _soa.template get<Particle_T::AttributeNames::forceZ>().d_view(i) = 0;
        });

    typedef Kokkos::View<double *[3], typename autopas::utils::kokkos::DeviceSpace::scratch_memory_space,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        block;

    auto tp = Kokkos::TeamPolicy<autopas::utils::kokkos::DeviceSpace>(_num_pairs, Kokkos::AUTO)
                  .set_scratch_size(0, Kokkos::PerTeam(block::shmem_size(BLOCK_SIZE)));

    Kokkos::parallel_for(
        "update force", tp,
        KOKKOS_CLASS_LAMBDA(const typename Kokkos::TeamPolicy<utils::kokkos::DeviceSpace>::member_type &team) {
          const size_t pair_id = team.league_rank();
          const size_t b1 = _pairs.d_view(pair_id, 0);
          const size_t b2 = _pairs.d_view(pair_id, 1);

          auto b1_start = b1 * BLOCK_SIZE;
          auto b1_end = Kokkos::min((b1 + 1) * BLOCK_SIZE, _soa.size());

          auto b2_start = b2 * BLOCK_SIZE;
          auto b2_end = Kokkos::min((b2 + 1) * BLOCK_SIZE, _soa.size());

          block block1(team.team_scratch(0), BLOCK_SIZE);

          // TODO Maybe theres a copy for that?
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, b1_start, b1_end), [&](const int64_t j) {
            block1(j - b1_start, X_AXIS) = _soa.template get<Particle_T::AttributeNames::posX>().d_view(j);
            block1(j - b1_start, Y_AXIS) = _soa.template get<Particle_T::AttributeNames::posY>().d_view(j);
            block1(j - b1_start, Z_AXIS) = _soa.template get<Particle_T::AttributeNames::posZ>().d_view(j);
          });
          team.team_barrier();
          if (b1 == b2) {
            _functor.KokkosSoAFunctor(team, _soa, block1, b1_start, b1_end);
          } else {
            for (uint64_t i = b1_start; i < b1_end; ++i) {
              auto mask = _pairs.d_view(pair_id, 2);
              if (mask & (1L << (i - b1_start))) {
                _functor.KokkosSoAFunctorPairwise(team, _soa, block1, i, b1_start, b2_start, b2_end);
              }
            }
          }
        });
  }

  KOKKOS_INLINE_FUNCTION
  double bbox_distance_sq(const double minx1, const double maxx1, const double miny1, const double maxy1,
                          const double minz1, const double maxz1, const double minx2, const double maxx2,
                          const double miny2, const double maxy2, const double minz2, const double maxz2) const {
    const auto axis_gap = [](const double min1, const double max1, const double min2, const double max2) {
      const double d1 = min2 - max1;
      const double d2 = min1 - max2;
      return Kokkos::max(0.0, Kokkos::max(d1, d2));
    };
    const double dx = axis_gap(minx1, maxx1, minx2, maxx2);
    const double dy = axis_gap(miny1, maxy1, miny2, maxy2);
    const double dz = axis_gap(minz1, maxz1, minz2, maxz2);
    return dx * dx + dy * dy + dz * dz;
  }

  void updatePairs() {
    auto tp1 = Kokkos::TeamPolicy(_num_blocks, Kokkos::AUTO);
    Kokkos::parallel_reduce(
        "count interacting _pairs", tp1,
        KOKKOS_CLASS_LAMBDA(const typename decltype(tp1)::member_type &team, uint64_t &local_sum) {
          const size_t b1 = team.league_rank();
          const size_t start = b1;
          const size_t end = _num_blocks;
          uint64_t sum = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, start, end),
              [&](const int64_t b2, uint64_t &inner_sum) {
                // TODO Don't calculate if b1 ==b2
                auto distance = bbox_distance_sq(
                    _bounding_box.d_view(b1, X_AXIS, AXIS_MIN), _bounding_box.d_view(b1, X_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b1, Y_AXIS, AXIS_MIN), _bounding_box.d_view(b1, Y_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b1, Z_AXIS, AXIS_MIN), _bounding_box.d_view(b1, Z_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b2, X_AXIS, AXIS_MIN), _bounding_box.d_view(b2, X_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b2, Y_AXIS, AXIS_MIN), _bounding_box.d_view(b2, Y_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b2, Z_AXIS, AXIS_MIN), _bounding_box.d_view(b2, Z_AXIS, AXIS_MAX));
                if (distance <= _cutoff * _cutoff) {
                  auto b1_start = b1 * BLOCK_SIZE;
                  auto b1_end = Kokkos::min((b1 + 1) * BLOCK_SIZE, _soa.size());
                  uint64_t mask = 0l;
                  for (int64_t i = b1_start; i < b1_end; ++i) {
                    auto pos_x = _soa.template get<Particle_T::AttributeNames::posX>().d_view(i);
                    auto pos_y = _soa.template get<Particle_T::AttributeNames::posY>().d_view(i);
                    auto pos_z = _soa.template get<Particle_T::AttributeNames::posZ>().d_view(i);
                    auto dist = bbox_distance_sq(
                        pos_x, pos_x, pos_y, pos_y, pos_z, pos_z, _bounding_box.d_view(b2, X_AXIS, AXIS_MIN),
                        _bounding_box.d_view(b2, X_AXIS, AXIS_MAX), _bounding_box.d_view(b2, Y_AXIS, AXIS_MIN),
                        _bounding_box.d_view(b2, Y_AXIS, AXIS_MAX), _bounding_box.d_view(b2, Z_AXIS, AXIS_MIN),
                        _bounding_box.d_view(b2, Z_AXIS, AXIS_MAX));
                    if (dist <= _cutoff * _cutoff) {
                      mask |= 1L << (i - b1_start);
                    }
                  }
                  int64_t signed_mask = 0L;
                  std::memcpy(&signed_mask, &mask, sizeof(uint64_t));
                  if (signed_mask > 0L) {
                    inner_sum += 1;
                  }
                }
              },
              Kokkos::Sum<uint64_t>(sum));

          Kokkos::single(Kokkos::PerTeam(team), [&]() {
            local_sum += sum;
            _block_pairs_count.d_view(b1) = sum;
          });
        },
        Kokkos::Sum<uint64_t>(_num_pairs));

    if (_pairs.extent(0) < _num_pairs) {
      printf("resizing to  %lu\n", _num_pairs);
      Kokkos::resize(_pairs, _num_pairs, 3);
    }

    Kokkos::parallel_scan(
        "calculate_start_indices", Kokkos::RangePolicy<utils::kokkos::DeviceSpace>(0, _num_blocks),
        KOKKOS_CLASS_LAMBDA(const int64_t b1, uint64_t &update, const bool final) {
          if (final) {
            _block_pairs_start_index.d_view(b1) = update;
          }
          const uint64_t block_count = _block_pairs_count.d_view(b1);
          update += block_count;
        });
    Kokkos::parallel_for(
        "write_interacting_pairs", Kokkos::TeamPolicy<typename utils::kokkos::DeviceSpace>(_num_blocks, Kokkos::AUTO),
        KOKKOS_CLASS_LAMBDA(const typename Kokkos::TeamPolicy<typename utils::kokkos::DeviceSpace>::member_type &team) {
          size_t b1 = team.league_rank();
          uint64_t current_index = _block_pairs_start_index.d_view(b1);

          Kokkos::parallel_scan(
              Kokkos::TeamThreadRange(team, b1, _num_blocks), [&](const int64_t b2, int64_t &num, bool final) {
                auto distance = bbox_distance_sq(
                    _bounding_box.d_view(b1, X_AXIS, AXIS_MIN), _bounding_box.d_view(b1, X_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b1, Y_AXIS, AXIS_MIN), _bounding_box.d_view(b1, Y_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b1, Z_AXIS, AXIS_MIN), _bounding_box.d_view(b1, Z_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b2, X_AXIS, AXIS_MIN), _bounding_box.d_view(b2, X_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b2, Y_AXIS, AXIS_MIN), _bounding_box.d_view(b2, Y_AXIS, AXIS_MAX),
                    _bounding_box.d_view(b2, Z_AXIS, AXIS_MIN), _bounding_box.d_view(b2, Z_AXIS, AXIS_MAX));
                if (distance <= _cutoff * _cutoff) {
                  size_t b1_start = b1 * BLOCK_SIZE;
                  size_t b1_end = Kokkos::min((b1 + 1) * BLOCK_SIZE, _soa.size());
                  uint64_t mask = 0l;

                  for (uint64_t i = b1_start; i < b1_end; ++i) {
                    auto pos_x = _soa.template get<Particle_T::AttributeNames::posX>().d_view(i);
                    auto pos_y = _soa.template get<Particle_T::AttributeNames::posY>().d_view(i);
                    auto pos_z = _soa.template get<Particle_T::AttributeNames::posZ>().d_view(i);
                    auto dist = bbox_distance_sq(
                        pos_x, pos_x, pos_y, pos_y, pos_z, pos_z, _bounding_box.d_view(b2, X_AXIS, AXIS_MIN),
                        _bounding_box.d_view(b2, X_AXIS, AXIS_MAX), _bounding_box.d_view(b2, Y_AXIS, AXIS_MIN),
                        _bounding_box.d_view(b2, Y_AXIS, AXIS_MAX), _bounding_box.d_view(b2, Z_AXIS, AXIS_MIN),
                        _bounding_box.d_view(b2, Z_AXIS, AXIS_MAX));
                    if (dist <= _cutoff * _cutoff) {
                      mask |= 1L << (i - b1_start);
                    }
                  }
                  int64_t signed_mask = 0L;
                  std::memcpy(&signed_mask, &mask, sizeof(uint64_t));

                  if (signed_mask > 0L) {
                    if (final) {
                      _pairs.d_view(current_index + num, 0) = b1;
                      _pairs.d_view(current_index + num, 1) = b2;
                      _pairs.d_view(current_index + num, 2) = signed_mask;
                    }
                    num += 1;
                  }
                }
              });
          team.team_barrier();
        });
  }

  template <size_t axis>
  KOKKOS_INLINE_FUNCTION constexpr void min_max_reduction(auto &team, auto start, auto end,
                                                          Kokkos::MinMax<double>::value_type &minmax) const {
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, start, end),
        [&](const uint64_t i, Kokkos::MinMax<double>::value_type &local_minmax) {
          const double coordinate = _soa.template get<axis>().d_view(i);
          local_minmax.min_val = Kokkos::min(local_minmax.min_val, coordinate);
          local_minmax.max_val = Kokkos::max(local_minmax.max_val, coordinate);
        },
        Kokkos::MinMax<double>(minmax));
  };

  void updateBoundingBox() {
    auto tp = Kokkos::TeamPolicy(_num_blocks, Kokkos::AUTO);
    Kokkos::parallel_for(
        "update bounding boxes", tp, KOKKOS_CLASS_LAMBDA(const typename decltype(tp)::member_type &team) {
          const auto block_id = team.league_rank();
          const size_t start = block_id * BLOCK_SIZE;
          const size_t end = Kokkos::min(start + BLOCK_SIZE, _soa.size());

          Kokkos::MinMax<double>::value_type minmax_x;
          Kokkos::MinMax<double>::value_type minmax_y;
          Kokkos::MinMax<double>::value_type minmax_z;

          min_max_reduction<Particle_T::AttributeNames::posX>(team, start, end, minmax_x);
          min_max_reduction<Particle_T::AttributeNames::posY>(team, start, end, minmax_y);
          min_max_reduction<Particle_T::AttributeNames::posZ>(team, start, end, minmax_z);

          Kokkos::single(Kokkos::PerTeam(team), [&]() {
            _bounding_box.d_view(block_id, X_AXIS, AXIS_MIN) = minmax_x.min_val - _skinSize;
            _bounding_box.d_view(block_id, X_AXIS, AXIS_MAX) = minmax_x.max_val + _skinSize;
            _bounding_box.d_view(block_id, Y_AXIS, AXIS_MIN) = minmax_y.min_val - _skinSize;
            _bounding_box.d_view(block_id, Z_AXIS, AXIS_MIN) = minmax_z.min_val - _skinSize;
            _bounding_box.d_view(block_id, Z_AXIS, AXIS_MAX) = minmax_z.max_val + _skinSize;
          });
        });

    _bounding_box.modify_device();
  }

  void clusterParticles() {
    Kokkos::parallel_for(
        "generate indices", Kokkos::RangePolicy(0, _soa.size()), KOKKOS_CLASS_LAMBDA(const int64_t &i) {
          std::array<int64_t, 3> index = {0, 0, 0};

          // auto cell_size = optimize_cell_size(_cutoff);

          index[0] = normalize(_soa.template get<Particle_T::AttributeNames::posX>().d_view(i), _cellSize[0]);
          index[1] = normalize(_soa.template get<Particle_T::AttributeNames::posY>().d_view(i), _cellSize[1]);
          index[2] = normalize(_soa.template get<Particle_T::AttributeNames::posZ>().d_view(i), _cellSize[2]);

          auto index_code = hilbert_index_3d(index[0], index[1], index[2]);
          _spatial_index.d_view(i) = index_code;
        });

    auto policy = Kokkos::RangePolicy(0, _spatial_index.extent(0));
    Kokkos::parallel_for(
        "initialize permutation index", policy, KOKKOS_CLASS_LAMBDA(const int64_t i) { _permutations.d_view(i) = i; });

    Kokkos::Experimental::sort_by_key(space, _spatial_index.d_view, _permutations.d_view);

    _soa.sort(space, _permutations.d_view);

    _rebuild_position.modify_device();
    _spatial_index.modify_device();
    _permutations.modify_device();
  }

  Functor_T _functor;
  double _cutoff;
  std::array<double, 3> _cellSize;
  double _skinSize;

  utils::kokkos::DeviceSpace space;

  template <typename scalar>
  using view_type = Kokkos::DualView<scalar, Kokkos::LayoutLeft, utils::kokkos::DeviceSpace>;

  utils::kokkos::KokkosSoA<typename Particle_T::SoAArraysType> _soa;

  view_type<double *[3][2]> _bounding_box;
  // BLOCK sized
  view_type<uint64_t *> _block_pairs_count;
  // BLOCK sized
  view_type<uint64_t *> _block_pairs_start_index;
  // DYNAMICALLY sized: maximal BLOCK^2 but usually linearly in number of interactions
  // Third stores mask
  view_type<int64_t *[3]> _pairs;
  // PARTICLE sized
  view_type<uint64_t *> _spatial_index;
  // PARTICLE sized
  view_type<double *[3]> _rebuild_position;
  // PARTICLE sized
  view_type<double *> _displacement;
  // PARTICLE sized
  view_type<int64_t *> _permutations;
  // Spans from 0.0, to value in each dimension
  view_type<double[3]> _boundaries;

  uint64_t _num_blocks;
  uint64_t _num_pairs;
};

}  // namespace autopas::containers::kokkos::traversal
