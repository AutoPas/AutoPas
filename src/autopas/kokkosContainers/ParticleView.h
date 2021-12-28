/**
 * @file ParticleView.h
 *
 * @date 2 Nov 2021
 * @author lgaertner
 */

#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include "autopas/cells/KokkosParticleCell.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

/**
 * ParticleView class.
 * This class uses a Kokkos::view to store particles and allows access to its references.
 * It also keeps the order of particles based on the currently used container.
 * @tparam ParticleType Type of the Particle that is being stored
 */
template <class ParticleType>
class ParticleView {
 public:
  ParticleView<ParticleType>() : _particleViewImp("ParticleView", _capacity){};

  void addParticle(const ParticleType &p) {
    _particleViewLock.lock();

    if (_size == _capacity) {
      _capacity *= 2;
      Kokkos::resize(_particleViewImp, _capacity);
    }
    deep_copy(subview(_particleViewImp, _size), p);
    _size++;

    _particleViewLock.unlock();
  }

  template <bool deleteDummy, typename Lambda>
  void binParticles(Lambda particleBinningLambda, Kokkos::View<autopas::KokkosParticleCell<ParticleType> *> cells,
                    std::string label = "binParticles") {
    // scan to create buckets
    Kokkos::View<size_t *> begin;
    Kokkos::View<size_t *> cellsize;
    auto nBuckets = cells.size();
    Kokkos::RangePolicy<> bucketRange(0, nBuckets);
    Kokkos::RangePolicy<> particleRange(0, _size);

    Kokkos::View<size_t *> counts("nParticles-per-cell.", nBuckets);

    auto counts_sv = Kokkos::Experimental::create_scatter_view(counts);

    // count number of particles per cell
    Kokkos::parallel_for(
        label + "compute-cellcount", particleRange, KOKKOS_LAMBDA(const size_t &i) {
          auto p = _particleViewImp[i];
          // if dummys should be deleted, check if particle is dummy and then disregard for binning
          if (deleteDummy and p.isDummy()) {
            Kokkos::atomic_sub(&_size, 1ul);
          } else {
            auto counts_data = counts_sv.access();
            counts_data(particleBinningLambda(p)) += 1;
          }
        });
    Kokkos::fence();

    Kokkos::Experimental::contribute(counts, counts_sv);

    // Scan offsets and compute begin cell indices
    Kokkos::parallel_scan(
        "offsetScan", bucketRange, KOKKOS_LAMBDA(const std::size_t c, int &update, const bool final_pass) {
          if (final_pass) cells[c].begin = update;
          update += counts[c];
        });
    Kokkos::fence();

    // Initialize cellsize with 0s
    Kokkos::parallel_for(
        bucketRange, KOKKOS_LAMBDA(const size_t i) { cells[i].cellSize = 0; }, label + "initialize_cellsize");
    Kokkos::fence();

    // compute permutation vector
    Kokkos::View<size_t *> permutes(label + "particle-permutation-view", _size);
    Kokkos::parallel_for(
        "", particleRange, KOKKOS_LAMBDA(const size_t &i) {
          auto p = _particleViewImp[i];
          if (not(deleteDummy and p.isDummy())) {
            size_t cellId = particleBinningLambda(_particleViewImp[i]);
            int c = Kokkos::atomic_fetch_add(&cells[cellId].cellSize, 1ul);
            permutes[cells[cellId].begin + c] = i;
          }
        });
    Kokkos::fence();

    // particle range without dummy particles
    Kokkos::RangePolicy<> newParticleRange(0, _size);

    // insert particles into copy
    Kokkos::View<ParticleType *> copyParticleListImpl(label + "intermediate-particles-sort-target-view", _size);
    Kokkos::parallel_for(
        "insert-into-copyParticleListImpl", newParticleRange,
        KOKKOS_LAMBDA(const size_t &i) { copyParticleListImpl[i] = _particleViewImp[permutes[i]]; });
    Kokkos::fence();

    // copy back
    Kokkos::parallel_for(
        label + "copyBack", newParticleRange,
        KOKKOS_LAMBDA(const size_t &i) { _particleViewImp[i] = copyParticleListImpl[i]; });
    Kokkos::fence();

    // set _particleViewImp pointer in cells
    Kokkos::parallel_for(
        label + "setParticlePtr", bucketRange,
        KOKKOS_LAMBDA(const size_t &i) { cells[i].particlesPtr = &_particleViewImp; });
    Kokkos::fence();
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, std::string label = "") {
    std::array<double, 3> dummy{};
    _forEach<parallel, false, false>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy,
                                     {0, _size}, dummy, dummy, label);
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
               std::string label = "") {
    _forEach<parallel, false, true>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy,
        {0, _size}, lowerCorner, higherCorner, label);
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior,
               std::string label = "ParticleView::forEach(behavior)") {
    std::array<double, 3> dummy{};
    _forEach<parallel, true, false>(forEachLambda, behavior, {0, _size}, dummy, dummy, label);
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
               autopas::IteratorBehavior behavior, std::string label = "") {
    _forEach<parallel, true, true>(forEachLambda, behavior, {0, _size}, lowerCorner, higherCorner,
                                   label);
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, autopas::KokkosParticleCell<ParticleType> cell, std::string label = "") {
    std::array<double, 3> dummy{};
    _forEach<parallel, false, false>(forEachLambda, autopas::IteratorBehavior::ownedOrHaloOrDummy,
                                     cell.getRange(), dummy, dummy, label);
  }

  template <bool parallel, typename Lambda>
  void forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior, autopas::KokkosParticleCell<ParticleType> cell,
               std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label = "") {
    _forEach<parallel, true, true>(forEachLambda, behavior, cell.getRange(), lowerCorner, higherCorner,
                                   label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, std::string label = "") {
    std::array<double, 3> dummy{};
    _reduce<false, false>(reduceLambda, result, autopas::IteratorBehavior::ownedOrHalo, Kokkos::RangePolicy<>(0, _size),
                          dummy, dummy, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
              std::string label = "") {
    _reduce<false, true>(reduceLambda, result, autopas::IteratorBehavior::ownedOrHalo, Kokkos::RangePolicy<>(0, _size),
                         lowerCorner, higherCorner, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, std::string label = "") {
    std::array<double, 3> dummy{};
    _reduce<true, false>(reduceLambda, result, behavior, Kokkos::RangePolicy<>(0, _size), dummy, dummy, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, std::array<double, 3> lowerCorner,
              std::array<double, 3> higherCorner, std::string label = "") {
    _reduce<true, true>(reduceLambda, result, behavior, Kokkos::RangePolicy<>(0, _size), lowerCorner, higherCorner,
                        label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::KokkosParticleCell<ParticleType> cell, std::string label = "") {
    std::array<double, 3> dummy{};
    _reduce<false, false>(reduceLambda, result, autopas::IteratorBehavior::ownedOrHalo, cell.getKokkosRangePolicy(),
                          dummy, dummy, label);
  }

  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior,
              autopas::KokkosParticleCell<ParticleType> cell, std::array<double, 3> lowerCorner,
              std::array<double, 3> higherCorner, std::string label = "") {
    _reduce<true, true>(reduceLambda, result, behavior, cell.getKokkosRangePolicy(), lowerCorner, higherCorner, label);
  }

  // TODO (lgaertner) reduce in region

  void deleteAll() { _size = 0; }

  size_t getSize() const { return _size; }
  size_t getCapacity() const { return _capacity; }
  Kokkos::View<ParticleType *> getParticles() const { return _particleViewImp; };

  void inspect() {
    printf("");
  }

 private:
  template <bool ownershipCheck, bool regionCheck, typename Lambda, typename A>
  void _reduce(Lambda reduceLambda, A &result, autopas::IteratorBehavior behavior, Kokkos::RangePolicy<> rangePolicy,
               std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label) {
    Kokkos::parallel_reduce(
        label, rangePolicy,
        KOKKOS_LAMBDA(const size_t &i, A &a) {
          if ((not ownershipCheck) or behavior.contains(_particleViewImp[i])) {
            if ((not regionCheck) or autopas::utils::inBox(_particleViewImp[i].getR(), lowerCorner, higherCorner)) {
              reduceLambda(_particleViewImp[i], a);
            }
          }
        },
        Kokkos::Sum<A>(result));
    Kokkos::fence();
  }

  template <bool parallel, bool ownershipCheck, bool regionCheck, typename Lambda>
  void _forEach(Lambda forEachLambda, autopas::IteratorBehavior behavior, std::array<size_t, 2> range,
                std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner, std::string label) {
    size_t forBegin;
    size_t forEnd;
    Kokkos::RangePolicy<> rangePolicy;

    if (not parallel) {
      forBegin = range[0];
      forEnd = range[1];
      rangePolicy = Kokkos::RangePolicy<>(0, 1);
    } else {
      forBegin = 0;
      forEnd = 1;
      rangePolicy = Kokkos::RangePolicy<>(range[0], range[1]);
    }

    Kokkos::parallel_for(
        label, rangePolicy, KOKKOS_LAMBDA(const size_t &i) {
          for (size_t j = forBegin; j < forEnd; j++) {
            auto index = parallel ? i : j;
            if ((not ownershipCheck) or behavior.contains(_particleViewImp[index])) {
              if ((not regionCheck) or
                  autopas::utils::inBox(_particleViewImp[index].getR(), lowerCorner, higherCorner)) {
                forEachLambda(_particleViewImp[index]);
              }
            }
          }
        });
    Kokkos::fence();
  }

  autopas::AutoPasLock _particleViewLock;

  size_t _capacity{8};
  size_t _size{0};
  Kokkos::View<ParticleType *> _particleViewImp;
};
