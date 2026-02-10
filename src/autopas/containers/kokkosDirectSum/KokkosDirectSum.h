/**
 * @file KokkosDirectSum.h
 * @date 31 October 2025
 * @author Luis Gall
 */

// This file is oriented on the implementation of the CellBasedParticleContainer, just that it does not use any cells to store particles

#pragma once

#include <Kokkos_Core.hpp>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utils/KokkosAoS.h"
#include "autopas/utils/KokkosSoA.h"
#include "autopas/utils/KokkosStorage.h"
#include "traversals/KokkosDsNaiveParallelTraversal.h"

namespace autopas {

template <class Particle_T>
    class KokkosDirectSum : public ParticleContainerInterface<Particle_T> {

        public:
            KokkosDirectSum(DataLayoutOption dataLayout, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double skin)
                : ParticleContainerInterface<Particle_T>(boxMin, boxMax, skin), _dataLayout(dataLayout) {
              if (dataLayout == DataLayoutOption::aos) {
                _aosUpToDate = true;
              }
              else if (dataLayout == DataLayoutOption::soa) {
                _soaUpToDate = true;
              }
              _ownedParticles.setLayout(dataLayout);
              _haloParticles.setLayout(dataLayout);
            }

            [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; }

            bool allowsKokkos() const override { return true; }

            void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
                
                if (numParticles > capacityOwned) {
                  _ownedParticles.resize(numParticles);
                  capacityOwned = numParticles;
                }

                if (numParticlesHaloEstimate > capacityHalo) {
                  _haloParticles.resize(numParticlesHaloEstimate);
                  capacityHalo = numParticlesHaloEstimate;
                }
            }

            // If possible, a method for adding a whole vector would make more sense (first touch policy, etc.)
            // TODO: provide an additional function for this

            void addParticleImpl(const Particle_T &p) override {

              std::lock_guard<AutoPasLock> guard (ownedLock);

              if (numberOfOwned >= capacityOwned) {
                // TODO: resize
                std::cout << "Problem" << std::endl;
              }
              _ownedParticles.addParticle(numberOfOwned++, p);
            }

            void addHaloParticleImpl(const Particle_T &haloP) override {

              std::lock_guard<AutoPasLock> guard(haloLock);

                if (numberOfHalo < capacityHalo) {
                    _haloParticles.addParticle(numberOfHalo++, haloP);
                } else {
                    // TODO: resize
                }
            }

            bool updateHaloParticle(const Particle_T &haloParticle) override {
              _haloParticles.template sync<HostSpace::execution_space>();

              const auto skinHalf = this->getVerletSkin() * 0.5;

              for (int i = 0; i < numberOfHalo; ++i) {
                const auto id = _haloParticles.template operator()<Particle_T::AttributeNames::id, true, true>(i);
                if (id == haloParticle.getID()) {
                  const typename Particle_T::ParticleSoAFloatPrecision x1 = _haloParticles.template operator()<Particle_T::AttributeNames::posX, true, true>(i);
                  const typename Particle_T::ParticleSoAFloatPrecision x2 = _haloParticles.template operator()<Particle_T::AttributeNames::posX, true, true>(i);
                  const typename Particle_T::ParticleSoAFloatPrecision x3 = _haloParticles.template operator()<Particle_T::AttributeNames::posX, true, true>(i);

                  const typename Particle_T::ParticleSoAFloatPrecision dX = x1 - haloParticle.getR().at(0);
                  const typename Particle_T::ParticleSoAFloatPrecision dY = x2 - haloParticle.getR().at(1);
                  const typename Particle_T::ParticleSoAFloatPrecision dZ = x3 - haloParticle.getR().at(2);

                  const typename Particle_T::ParticleSoAFloatPrecision distanceSqr = dX*dX + dY*dY + dZ*dZ;

                  if (distanceSqr < skinHalf*skinHalf) {
                    _haloParticles.addParticle(i, haloParticle);
                    _haloParticles.template markModified<HostSpace::execution_space>();
                    return true;
                  }
                }
              }

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
              if (behavior & 0b1) {
                number += numberOfOwned;
              }
              if (behavior & 0b10) {
                number += numberOfHalo;
              }
              // TODO: other behaviors (Actually, dummies can be somewhere in both? lists)
              // maybe find dummies in owned list with the help of reduceKokkos()
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

                // TODO: if this is the common structure, why isn't this generalized and called in a higher level of the hierarchy?
                traversal->initTraversal();
                traversal->traverseParticles();
                traversal->endTraversal();

                finishTraversal(traversal);
            }

            [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {
              if (keepNeighborListsValid) {
                // TODO: collect particles and mark non owned as dummy
                // i.e.: those particles which are outside of the box shall be returned, halo particles should be marked as dummies
                return {};
              }
              deleteHaloParticles();
              // TODO: delete dummy particles
              // TODO: determine those particles which are not inside of the box and return them

              return {};
            }

            [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
                // TODO
                return TraversalSelectorInfo{};
            }

            /* TODO: Begin of Code Crime (investigate generalizations of this -> every container does EXACTLY the same) */
            [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
                IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) override {
                // Copy from DirectSum.h
                convertToAoS();
                _soaUpToDate = false;
                return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
                IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors =
                    nullptr) const override {
                // Copy from DirectSum.h
                // TODO: think about how to handle case when particles are stored in SoA Format, or mark this as discouraged
                return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
            }

            [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                typename ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors = nullptr) override {
                // Copy from DirectSum.h
                convertToAoS();
                _soaUpToDate = false;
                return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                typename ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors =
                    nullptr) const override {
                // Copy from DirectSum.h
                // TODO: think about how to handle case when particles are stored in SoA Format, or mark this as discouraged
                return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            /* TODO: End of Code Crime */

            template <typename Lambda>
            void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                    const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
                // TODO
            }

            template <typename Lambda>
            void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
              // TODO
            }

            template <class ExecSpace, typename Lambda>
            void forEachKokkos(Lambda& forEachLambda, IteratorBehavior behavior) {

              if (_dataLayout == DataLayoutOption::aos) {
                convertToAoS();
                _soaUpToDate = false;
              }
              else if (_dataLayout == DataLayoutOption::soa) {
                convertToSoA();
                _ownedParticles.template sync<ExecSpace>();
                _haloParticles.template sync<ExecSpace>();
                _aosUpToDate = false;
              }

              // TODO: make sure that no AoS runs on Cuda
              if (behavior & 0b1) {
                auto& owned = _ownedParticles;
                Kokkos::parallel_for("forEachKokkosOwned", Kokkos::RangePolicy<ExecSpace>(0, numberOfOwned), KOKKOS_LAMBDA(int i)  {
                  forEachLambda(i, owned);
                });
              }
              if (behavior & 0b10) {
                auto& halo = _haloParticles;
                Kokkos::parallel_for("forEachKokkosHalo", Kokkos::RangePolicy<ExecSpace>(0, numberOfHalo), KOKKOS_LAMBDA(int i)  {
                  forEachLambda(i, halo);
                });
              }
              // TODO: consider other behavior such as dummies, container only, ...

              if (_dataLayout == DataLayoutOption::soa) {
                _ownedParticles.template markModified<ExecSpace>();
                _haloParticles.template markModified<ExecSpace>();
              }
            }

            template<class ExecSpace, typename Result, typename Reduction, typename Lambda>
            void reduceKokkos(Lambda reduceLambda, Result& result, IteratorBehavior behavior) {

              if (_dataLayout == DataLayoutOption::aos) {
                convertToAoS();
              }
              else if (_dataLayout == DataLayoutOption::soa) {
                convertToSoA();
                _ownedParticles.template sync<ExecSpace>();
                _haloParticles.template sync<ExecSpace>();
              }

              if (behavior & 0b1) {
                auto& owned = _ownedParticles;
                Kokkos::parallel_reduce("reduceKokkosOwned", Kokkos::RangePolicy<ExecSpace>(0, numberOfOwned), KOKKOS_LAMBDA(int i, Result& localResult)  {
                  reduceLambda(i, owned, localResult);
                }, Reduction(result));
              }
              if (behavior & 0b10) {
                auto& halo = _haloParticles;
                Kokkos::parallel_reduce("reduceKokkosHalo", Kokkos::RangePolicy<ExecSpace>(0, numberOfHalo), KOKKOS_LAMBDA(int i, Result& localResult)  {
                  reduceLambda(i, halo, localResult);
                }, Reduction(result));
              }
              // TODO: consider other behavior such as dummies, container only, ...
            }

            [[nodiscard]] double getCutoff() const final { return 0; }

            void setCutoff(double) final {};

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
              if (cellIndex == 0 && particleIndex > numberOfOwned) {
                // Particle does not exist in owned list
                return false;
              } else if (cellIndex == 1 && particleIndex > numberOfHalo) {
                // Particle does not exist in halo list
                return false;
              }
              auto& storage = (cellIndex == 0) ? _ownedParticles : _haloParticles;
              storage.template operator()<Particle_T::AttributeNames::ownershipState, true, true>(particleIndex) = OwnershipState::dummy;
              return true;
            }
    
    private:
        
        template <typename Traversal>
        void prepareTraversal(Traversal &traversal) {
            auto* kokkosDsTraversal = dynamic_cast<DSKokkosTraversalInterface<Particle_T> *>(traversal);
            if (kokkosDsTraversal) {

              auto targetLayout = traversal->getDataLayout();

              if (targetLayout != _dataLayout) {
                convertTo(targetLayout);
              }

              kokkosDsTraversal->setOwnedToTraverse(_ownedParticles);
              kokkosDsTraversal->setHaloToTraverse(_haloParticles);
            }
            else {
                utils::ExceptionHandler::exception("The selected traversal is not compatible with the KokkosDirectSum container.");
            }
        }

        template <typename Traversal>
        void finishTraversal(Traversal traversal) {
          auto* kokkosDsTraversal = dynamic_cast<DSKokkosTraversalInterface<Particle_T> *>(traversal);

          if (kokkosDsTraversal) {

            auto targetLayout = traversal->getDataLayout();

            kokkosDsTraversal->retrieveOwned(_ownedParticles);

            if (targetLayout == DataLayoutOption::aos) {
              _soaUpToDate = false;
              _aosUpToDate = true;
              if (_dataLayout != targetLayout) {
                convertToSoA();
              }
            }
            else if (targetLayout == DataLayoutOption::soa) {
              _aosUpToDate = false;
              _soaUpToDate = true;
              if (_dataLayout != targetLayout) {
                convertToAoS();
              }
            }
          }
          else {
            utils::ExceptionHandler::exception("The selected traversal is not compatible with the KokkosDirectSum container.");
          }
        }

        void convertTo(DataLayoutOption layout) {
          if (layout == DataLayoutOption::soa) {
            convertToSoA();
          }
          else if (layout == DataLayoutOption::aos) {
            convertToAoS();
          }
        }

        void convertToSoA() {

          if (_soaUpToDate) {
            return;
          }

          _ownedParticles.convertToSoA(numberOfOwned);
          _haloParticles.convertToSoA(numberOfHalo);

          _soaUpToDate = true;
        }

        void convertToAoS() {

          if (_aosUpToDate) {
            return;
          }

          _ownedParticles.convertToAoS(numberOfOwned);
          _haloParticles.convertToAoS(numberOfHalo);

          _aosUpToDate = true;
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
          const auto& viewToCheck = (cellIndex == 0) ? _ownedParticles.getAoS() : _haloParticles.getAoS();

          if (particleIndex >= sizeToCheck or
              not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            viewToCheck.getParticle(particleIndex), iteratorBehavior, boxMin, boxMax)) {

            std::tie(cellIndex, particleIndex) =
            advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
          }

          if (cellIndex >= 2) {
            return {nullptr, 0, 0};
          }

          const auto& viewToExtract = (cellIndex == 0) ? _ownedParticles.getAoS() : _haloParticles.getAoS();
          const Particle_T * result = &viewToExtract.getParticle(particleIndex);

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
              ((cellIndex == 0) ? _ownedParticles.getAoS().getParticle(particleIndex) : _haloParticles.getAoS().getParticle(particleIndex)), iteratorBehavior, boxMin, boxMax));

          // the indices returned at this point should always be valid
          return {cellIndex, particleIndex};
        }

#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif

  using HostSpace = Kokkos::HostSpace;

        mutable AutoPasLock ownedLock {};

        mutable AutoPasLock haloLock {};

        DataLayoutOption _dataLayout {};

        bool _aosUpToDate = false;

        bool _soaUpToDate = false;

        utils::KokkosStorage<Particle_T> _ownedParticles {};

        utils::KokkosStorage<Particle_T> _haloParticles {};

        size_t numberOfOwned {0};

        size_t numberOfHalo {0};

        size_t capacityOwned {0};

        size_t capacityHalo {0};

        };
}