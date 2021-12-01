/**
 * @file ParticleView.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

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

    _dirty = true;
    _particleListLock.unlock();
  }

  template <typename Lambda>
  void binParticles(Lambda &particleBinningLambda, Kokkos::View<size_t *> &begin, Kokkos::View<size_t *> &cellsize,
                    std::string label = "") {
//    Kokkos::View<ParticleType> particleListImp = this->_particleListImp;

    if (_dirty) {
      // scan to create buckets
      auto nBuckets = begin.size();
      Kokkos::RangePolicy<> bucketRange(0, nBuckets);
      Kokkos::RangePolicy<> particleRange(0, _size);

      Kokkos::View<size_t *> counts("Temporary view to count nParticles per bucket.", nBuckets);

      auto counts_sv = Kokkos::Experimental::create_scatter_view(counts);

      // TODO
      Kokkos::parallel_for(
          particleRange,
          KOKKOS_LAMBDA(const size_t &i) {
            auto counts_data = counts_sv.access();
            counts_data(particleBinningLambda(_particleListImp[i])) += 1;
          },
          label + "cell_count");

      Kokkos::fence();
      Kokkos::Experimental::contribute(counts, counts_sv);

      //
      auto offset_scan = KOKKOS_LAMBDA(const std::size_t c, int &update, const bool final_pass) {
        if (final_pass) begin[c] = update;
        update += counts[c];
      };
      Kokkos::parallel_scan("Cabana::LinkedCellList::build::offset_scan", bucketRange, offset_scan);
      Kokkos::fence();

      Kokkos::deep_copy(cellsize, 0ul);

      // compute permutation vector
      Kokkos::View<size_t *> permutes(label + "permutation view", _size);
      auto permutationLambda = KOKKOS_LAMBDA(const size_t &i) {
        size_t cellId = particleBinningLambda(_particleListImp[i]);
//        size_t cellId = 0;
        int c = Kokkos::atomic_fetch_add(&cellsize[i], 1ul);
        permutes[begin[cellId] + c] = i;
      };

      // insert particles into copy
      Kokkos::View<ParticleType *> copyParticleListImpl(label + "intermediate particles copy target view", _size);
      Kokkos::parallel_for(
          particleRange, KOKKOS_LAMBDA(const size_t &i) {
            copyParticleListImpl[i] = _particleListImp[permutes[i]];
          },
          "");

      // copy back
      Kokkos::parallel_for(
          particleRange, KOKKOS_LAMBDA(const size_t &i) {
            _particleListImp[i] = copyParticleListImpl[i];
                                         }, "");
    }
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, std::string label = "") {
    std::array<double, 3> dummy{};
    _forEach<false, false>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy, 0, _size, dummy, dummy, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
               std::string label = "") {
    _forEach<false, true>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy, 0, _size, lowerCorner,
                          higherCorner, label);
  }

  // TODO (lgaertner): temporary solution until container binning is implemented (maybe keep as backup?)
  template <typename Lambda>
  void forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior,
               std::string label = "ParticleView::forEach(behavior)") {
    std::array<double, 3> dummy{};
    _forEach<true, false>(forEachLambda, behavior, 0, _size, dummy, dummy, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
               autopas::IteratorBehavior behavior, std::string label = "") {
    _forEach<true, true>(forEachLambda, behavior, 0, _size, lowerCorner, higherCorner, label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, size_t begin, size_t cellSize, std::string label = "") {
    std::array<double, 3> dummy{};
    _forEach<false, false>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy, begin, cellSize, dummy, dummy,
                           label);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior, size_t begin, size_t cellSize,
               std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label = "") {
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
        label, Kokkos::RangePolicy<>(begin, begin + cellSize), KOKKOS_LAMBDA(const size_t &i) {
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
