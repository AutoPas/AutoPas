#pragma once

namespace autopas::containers::kokkos::traversal {
template <typename Particle_T>
class KokkosTraversalInterface {
 public:
  virtual ~KokkosTraversalInterface() = default;

 private:
  virtual void rebuild() = 0;
};
}  // namespace autopas::containers::kokkos::traversal