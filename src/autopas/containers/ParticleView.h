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
    std::array<double, 3> dummy{};
    _forEach<false, false>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy, 0, _size, dummy, dummy, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label = "") {
    _forEach<false, true>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy, 0, _size, lowerCorner, higherCorner, label);
  }

  // TODO (lgaertner): temporary solution until container binning is implemented (maybe keep as backup?)
  template <typename Lambda>
  void forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior,
               std::string label = "ParticleView::forEach(behavior)") {
    std::array<double, 3> dummy{};
    _forEach<true, false>(forEachLambda, behavior, 0, _size, dummy, dummy, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, autopas::IteratorBehavior behavior, std::string label = "") {
    _forEach<true, true>(forEachLambda, behavior, 0, _size, lowerCorner, higherCorner, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, size_t begin, size_t cellSize, std::string label = "") {
    std::array<double, 3> dummy{};
    _forEach<false, false>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy, begin, cellSize, dummy, dummy, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior, size_t begin, size_t cellSize, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label = "") {
    _forEach<true, true>(forEachLambda, behavior, begin, cellSize, lowerCorner, higherCorner, label);
  }





  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, std::string label = "") {
    std::array<double, 3> dummy{};
    _reduce<false, false>(reduceLambda, result, autopas::IteratorBehavior::ownedOrHalo, 0ul, _size, dummy, dummy,
                          label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
              std::string label = "") {
    _reduce<false, true>(reduceLambda, result, autopas::IteratorBehavior::ownedOrHalo, 0ul, _size, lowerCorner,
                         higherCorner, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, std::string label = "") {
    std::array<double, 3> dummy{};
    _reduce<true, false>(reduceLambda, result, behavior, 0ul, _size, dummy, dummy, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, std::array<double, 3> lowerCorner,
              std::array<double, 3> higherCorner, std::string label = "") {
    _reduce<true, true>(reduceLambda, result, behavior, 0ul, _size, lowerCorner, higherCorner, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, size_t begin, size_t cellSize, std::string label = "") {
    std::array<double, 3> dummy{};
    _reduce<false, false>(reduceLambda, result, autopas::IteratorBehavior::ownedOrHalo, begin, cellSize, dummy, dummy,
                          label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, size_t begin, size_t cellSize,
              std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label = "") {
    _reduce<true, true>(reduceLambda, result, behavior, begin, cellSize, lowerCorner, higherCorner, label);
  }

  // TODO (lgaertner) reduce in region

  void deleteAll() {}

  size_t getSize() const { return _size; }

 private:
  template <bool ownershipCheck, bool regionCheck, typename Lambda, typename A>
  void _reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, size_t begin, size_t cellSize,
               std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label) {
    Kokkos::parallel_reduce(
        label, Kokkos::RangePolicy<>(begin, begin + cellSize),
        KOKKOS_LAMBDA(const size_t &i, A &a) {
          if ((not ownershipCheck) or behavior.contains(_particleListImp[i])) {
            if ((not regionCheck) or autopas::utils::inBox(_particleListImp[i].getR(), lowerCorner, higherCorner)) {
              reduceLambda(_particleListImp[i], a);
            }
          }
        },
        Kokkos::Sum<A>(result));
  }

  template <bool ownershipCheck, bool regionCheck, typename Lambda>
  void _forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior, size_t begin, size_t cellSize,
               std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label) {
    Kokkos::parallel_for(
        label, Kokkos::RangePolicy<>(begin, begin + cellSize),
        KOKKOS_LAMBDA(const size_t &i) {
          if ((not ownershipCheck) or behavior.contains(_particleListImp[i])) {
            if ((not regionCheck) or autopas::utils::inBox(_particleListImp[i].getR(), lowerCorner, higherCorner)) {
              forEachLambda(_particleListImp[i]);
            }
          }
        });
  }

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
