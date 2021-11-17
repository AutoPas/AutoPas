/**
 * @file ParticleView.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include <Kokkos_Core.hpp>

/**
 * ParticleView class.
 * This class uses a Kokkos::view to store particles and allows access to its references.
 * It also keeps the order of particles based on the currently used container.
 * @tparam ParticleType Type of the Particle that is being stored
 */
template <class ParticleType>
class ParticleView {
 public:
  ParticleView<ParticleType>() : _particleListImp("ParticleView", _capacity){};

  void addParticle(const ParticleType &p) {
    _particleListLock.lock();

    if (_size == _capacity) {
      _capacity *= 2;
      Kokkos::resize(_particleListImp, _capacity);
    }
    deep_copy(subview(_particleListImp, _size), p);
    _size++;

    _particleListLock.unlock();
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, std::string label = "") {
    Kokkos::parallel_for(
        label, Kokkos::RangePolicy<>(0, _size), KOKKOS_LAMBDA(const size_t &i) { forEachLambda(_particleListImp[i]); });
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, size_t begin, size_t cellSize, std::string label = "") {
    Kokkos::parallel_for(
        label, Kokkos::RangePolicy<>(begin, begin + cellSize),
        KOKKOS_LAMBDA(const size_t &i) { forEachLambda(_particleListImp[i]); });
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, std::string label = "") {
    Kokkos::parallel_reduce(
        label, Kokkos::RangePolicy<>(0, _size),
        KOKKOS_LAMBDA(const size_t &i, A &a) { reduceLambda(_particleListImp[i], a); }, Kokkos::Sum<A>(result));
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, size_t begin, size_t cellSize, std::string label = "") {
    Kokkos::parallel_reduce(
        label, Kokkos::RangePolicy<>(begin, begin + cellSize),
        KOKKOS_LAMBDA(const size_t &i, A &a) { reduceLambda(_particleListImp[i], a); }, Kokkos::Sum<A>(result));
  }

  void deleteAll() {}

  size_t getSize() const { return _size; }

 private:
  /**
   * Flag indicating whether there are out-of-date references in the vector.
   */
  bool _dirty = false;
  //  size_t _dirtyIndex{0};

  size_t _capacity{8};
  size_t _size{0};

  autopas::AutoPasLock _particleListLock;
  Kokkos::View<ParticleType *> _particleListImp;
};
