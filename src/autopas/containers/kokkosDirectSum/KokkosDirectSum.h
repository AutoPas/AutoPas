/**
 * @file KokkosDirectSum.h
 * @date 31 October 2025
 * @author Luis Gall
 */

// This file is oriented on the implementation of the CellBasedParticleContainer

#pragma once

#include <Kokkos_Core.hpp>

#include "autopas/containers/ParticleContainerInterface.h"
#include "traversals/KokkosDsNaiveParallelTraversal.h"

namespace autopas {

const size_t INIT_PARTICLES = 0;

template <class Particle_T>
    class KokkosDirectSum : public ParticleContainerInterface<Particle_T> {

        public:
            KokkosDirectSum(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double skin)
                : ParticleContainerInterface<Particle_T>(boxMin, boxMax, skin)
                //_ownedParticles(Kokkos::view_alloc(Kokkos::WithoutInitializing, "owned"), INIT_PARTICLES),
                //_haloParticles(Kokkos::view_alloc(Kokkos::WithoutInitializing, "halos"), INIT_PARTICLES/4)
                {}

            [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; }

            void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
                
                if (numParticles > capacityOwned) {
                    Kokkos::realloc(_ownedParticles, numParticles);
                  capacityOwned = numParticles;
                }

                if (numParticlesHaloEstimate > capacityHalo) {
                    Kokkos::realloc(_haloParticles, numParticlesHaloEstimate);
                  capacityHalo = numParticlesHaloEstimate;
                }
            }

            // If possible, a method for adding a whole vector would make more sense (first touch policy, etc.)
            // TODO: provide an additional function for this

            void addParticleImpl(const Particle_T &p) override {
              std::lock_guard<AutoPasLock> guard (ownedLock);

                if (numberOfOwned < capacityOwned) {
                    _ownedParticles(numberOfOwned++) = p;
                } else {
                    // TODO: resize
                }
            }

            void addHaloParticleImpl(const Particle_T &haloP) override {

              std::lock_guard<AutoPasLock> guard(haloLock);

                if (numberOfHalo < capacityHalo) {
                    _haloParticles(numberOfHalo++) = haloP;
                } else {
                    // TODO: resize
                }
            }

            bool updateHaloParticle(const Particle_T &haloParticle) override {
                // TODO

                return false;
            }

            void deleteHaloParticles() override {
                numberOfHalo = 0;
            }

            void deleteAllParticles() override {
                numberOfHalo = 0;
                numberOfOwned = 0;
            }

            size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
                size_t number = 0;

                return number;
            }

            size_t size() const override {
                return numberOfHalo + numberOfOwned;
            }

            void rebuildNeighborLists(TraversalInterface *traversal) override {
                // No-Op
            }

            void computeInteractions(TraversalInterface *traversal) override {
                prepareTraversal(traversal);

                // TODO: if this is the common structure, why isn't this generalized to some extent?
                traversal->initTraversal();
                traversal->traverseParticles();
                traversal->endTraversal();
            }

            [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {
                //TODO
                return {};
            }

            [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
                // TODO
                return TraversalSelectorInfo();
            }

            /* TODO: Begin of Code Crime (investigate generalizations of this) */
            [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
                IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) override {
                // Copy from DirectSum.h
                return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
                IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors =
                    nullptr) const override {
                // Copy from DirectSum.h
                return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
            }

            [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                typename ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors = nullptr) override {
                // Copy from DirectSum.h
                return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                typename ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors =
                    nullptr) const override {
                // Copy from DirectSum.h
                return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            /* TODO: End of Code Crime */

            template <typename Lambda>
            void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                    const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
                // TODO
            }

            // TODO: maybe there should also be a forEach

            [[nodiscard]] double getCutoff() const final { return 0; }

            void setCutoff(double cutoff) final { /* TODO */ };

            [[nodiscard]] double getInteractionLength() const final { return 0; }

            std::tuple<const Particle_T*, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                    IteratorBehavior iteratorBehavior) const final {

              constexpr std::array<double, 3> boxMin{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                                           std::numeric_limits<double>::lowest()};
              constexpr std::array<double, 3> boxMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                           std::numeric_limits<double>::max()};
              return getParticleImpl<false>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
            }

            std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                                     IteratorBehavior iteratorBehavior,
                                                                     const std::array<double, 3> &boxMin,
                                                                     const std::array<double, 3> &boxMax) const final {
                return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
            }

            bool deleteParticle(Particle_T &particle) final {
                // TODO
                return false;
            }

            bool deleteParticle(size_t cellIndex, size_t particleIndex) final {
                // TODO
                return false;
            }
    
    private:
        
        template <typename Traversal>
        void prepareTraversal(Traversal &traversal) {
            auto* kokkosDsTraversal = dynamic_cast<DSKokkosTraversalInterface<Particle_T> *>(traversal);
            if (kokkosDsTraversal) {
                kokkosDsTraversal->setOwnedToTraverse(_ownedParticles);
                kokkosDsTraversal->setHaloToTraverse(_haloParticles);
            }
            else {
                utils::ExceptionHandler::exception("The selected traversal is not compatible with the KokkosDirectSum container.");
            }
        }

        // Nearly a clone of DirectSum::getParticleImpl to guarantee that we fulfill all the requirements from getParticle in the ParticleContainerInterface's method description
        template <bool regionIter>
        std::tuple<const Particle_T*, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
          IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) const {

          const auto [startCellIndex, endCellIndex] = [&]() -> std::tuple<size_t, size_t> {
            // shortcuts to limit the iterator to only part of the domain.
            if (not(iteratorBehavior & IteratorBehavior::halo)) {
              // only owned region
              return {0, 0};
            }
            if (not(iteratorBehavior & IteratorBehavior::owned)) {
              // only halo region
              return {1, 1};
            }
            if constexpr (regionIter) {
              // if the region lies fully within the container only look at the owned cell
              if (utils::ArrayMath::less(this->getBoxMin(), boxMin) and utils::ArrayMath::less(boxMax, this->getBoxMax())) {
                return {0, 0};
              }
            }
            // all cells
            return {0, 1};
          }();

          if (cellIndex == 0 and particleIndex == 0) {
            cellIndex =
                startCellIndex + ((iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num());
          }

          if (cellIndex >= 2) {
            return {nullptr, 0, 0};
          }

          size_t sizeToCheck = (cellIndex == 0) ? numberOfOwned : numberOfHalo;
          const Kokkos::View<Particle_T*>& viewToCheck = (cellIndex == 0) ? _ownedParticles : _haloParticles;

          if (particleIndex >= sizeToCheck or
              not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            viewToCheck(particleIndex), iteratorBehavior, boxMin, boxMax)) {

            std::tie(cellIndex, particleIndex) =
            advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
          }

          if (cellIndex >= 2) {
            return {nullptr, 0, 0};
          }

          const Kokkos::View<Particle_T*>& viewToExtract = (cellIndex == 0) ? _ownedParticles : _haloParticles;
          const Particle_T * result = &viewToExtract(particleIndex);

          return {result, cellIndex, particleIndex};
        }

  template <bool regionIter>
  std::tuple<size_t, size_t> advanceIteratorIndices(size_t cellIndex, size_t particleIndex,
                                                    IteratorBehavior iteratorBehavior,
                                                    const std::array<double, 3> &boxMin,
                                                    const std::array<double, 3> &boxMax) const {
          // Find the indices for the next particle
          const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();

          do {
            // advance to the next particle
            ++particleIndex;

            // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.
            while (particleIndex >= (cellIndex == 0 ? numberOfOwned : numberOfHalo)) {
              cellIndex += stride;
              particleIndex = 0;
              // If there are no more reasonable cells return invalid indices.
              if (cellIndex > ((not(iteratorBehavior & IteratorBehavior::halo)) ? 0 : 1)) {
                return {std::numeric_limits<decltype(cellIndex)>::max(), std::numeric_limits<decltype(particleIndex)>::max()};
              }
            }
          } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
              ((cellIndex == 0) ? _ownedParticles(particleIndex) : _haloParticles(particleIndex)), iteratorBehavior, boxMin, boxMax));

          // the indices returned at this point should always be valid
          return {cellIndex, particleIndex};
        }

        //TODO: maybe think about a threshold above which particles are stored no longer on the GPU
        // However, for now, initialize and store particles on the host and offload them
        using MemSpace = Kokkos::HostSpace;

        mutable AutoPasLock ownedLock {};

        mutable AutoPasLock haloLock {};

        /* AoS Data structure */
        Kokkos::View<Particle_T*, MemSpace> _ownedParticles {};

        Kokkos::View<Particle_T*, MemSpace> _haloParticles {};

        /* SoA data structure */

        // TODO: implement

        /* AoSoA data structure */

        // TODO: implement

        size_t numberOfOwned {0};

        size_t numberOfHalo {0};

        size_t capacityOwned {INIT_PARTICLES};

        size_t capacityHalo {INIT_PARTICLES/4};

        };
}