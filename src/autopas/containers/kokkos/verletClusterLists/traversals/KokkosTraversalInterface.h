#pragma once

namespace autopas::containers::kokkos::traversal {
template <typename Particle_T>
class KokkosTraversalInterface {
 public:
  virtual ~KokkosTraversalInterface() = default;
  virtual void rebuild(std::vector<Particle_T> &particles) = 0;
};
}  // namespace autopas::containers::kokkos::traversal