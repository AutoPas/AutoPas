/**
 * @file KokkosDirectSum.h
 * @date 31 October 2025
 * @author Luis Gall
 */

// This file is oriented on the implementation of the CellBasedParticleContainer, just that it does not use any cells to store particles

#pragma once

#include <Kokkos_Core.hpp>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utilsKokkos/KokkosStorage.h"
#include "autopas/containers/kokkosDirectSum/traversals/DSKokkosTraversalInterface.h"

namespace autopas {

/**
 * This class stores all owned particles in a single list (meaning no cells or spatial decomposition)
 * The particle interactions are calculated directly, such that each particle interacts with every other particle
 * As a consequence, we discourage the use of this container for larger numbers of particles
 *
 * So far, the surrounding volume is not decomposed and also stored in a single list of particles (this can be subject
 * to change in future implementations)
 *
 * The most interesting feature of this class is that it allows to store particle data also in the SoA format and not
 * only AoS. For further information, we refer to the KokkosStorage class
 *
 * @tparam Particle_T AoS Particle type that is used with this container
 */
template <class Particle_T>
    class KokkosDirectSum : public ParticleContainerInterface<Particle_T> {

        public:
            KokkosDirectSum(DataLayoutOption dataLayout, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double skin)
                : ParticleContainerInterface<Particle_T>(boxMin, boxMax, skin), _dataLayout(dataLayout) {
              _ownedParticles.setLayout(dataLayout);
              _haloParticles.setLayout(dataLayout);
            }

            [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; }

            bool allowsKokkos() const override { return true; }

            void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {

              _ownedParticles.resize(numParticles);
              _haloParticles.resize(numParticlesHaloEstimate);
            }

            // If possible, a method for adding a whole vector would make more sense (first touch policy, etc.)
            // TODO: provide an additional function for this

            bool updateHaloParticle(const Particle_T &haloParticle) override {
              _haloParticles.template sync<HostSpace::execution_space>();

              const auto skinHalf = this->getVerletSkin() * 0.5;

              for (int i = 0; i < _haloParticles.size(); ++i) {
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
                    this->addHaloParticleImpl(haloParticle);
                    _haloParticles.template markModified<HostSpace::execution_space>();
                    return true;
                  }
                }
              }

              return false;
            }

            void deleteHaloParticles() override {
              _haloParticles.clear();
            }

            void deleteAllParticles() override {
              _ownedParticles.clear();
              _haloParticles.clear();
            }

            size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
              // TODO: this should maybe better just count the number of particles in both lists that fulfill the behavior requirement
              size_t number = 0;
              return number;
            }

            size_t size() const override {
              // TODO: double-check with the other containers how they handle that
              return -1;
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

            /* TODO: the current approach of massive contention locks is certainly not optimal but we can avoid costly data transfers to the CPU by doing everything on the GPU
                For the future, it might make sense to improve here if it turns out that also for the other containers (!), this becomes a crucial bottleneck */
            [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {

              /* Prepare data structures */
              utils::KokkosStorage<Particle_T> migrants {_ownedParticles.getLayout(), _ownedParticles.size()};

              auto& boxMin = ParticleContainerInterface<Particle_T>::_boxMin;
              auto& boxMax = ParticleContainerInterface<Particle_T>::_boxMax;

              const auto& boxMinKokkos = Kokkos::Array<double, 3> {boxMin.at(0), boxMin.at(1), boxMin.at(2)};
              const auto& boxMaxKokkos = Kokkos::Array<double, 3> {boxMax.at(0), boxMax.at(1), boxMax.at(2)};

              auto& owned = _ownedParticles;

              Kokkos::View<int*, DeviceSpace> migrantCounter {"migrantCounter", 1};

              if (keepNeighborListsValid) {
                // i.e.: those particles which are outside of the box shall be returned, halo particles should be marked as dummies
                // TODO: think about the validity of this approach because in theory, the particles should still exist for neighbor lists being valid (but this is not the case for DirectSum as there are no neighbor lists)
                _haloParticles.clear();

                // TODO: make sure that owned is synced to the right memory space
                Kokkos::parallel_for("collectMigrants", Kokkos::RangePolicy<DeviceSpace::execution_space>(0, _ownedParticles.size()), KOKKOS_LAMBDA(int i) {

                  // TODO: make sure this logic is working, i.e. the correct particles are found, i.e. the owned particles outside of the box
                  if ((not owned.template fulfillsIteratorRequirements<true, useHostView>(i, IteratorBehavior::ownedOrHaloOrDummy, boxMinKokkos, boxMaxKokkos))
                     and owned.template fulfillsIteratorRequirements<false, useHostView>(i, IteratorBehavior::owned, boxMinKokkos, boxMaxKokkos)) {

                    int migrantIndex = Kokkos::atomic_fetch_inc(&migrantCounter(0));
                    auto id = owned.template operator()<Particle_T::AttributeNames::id, true, useHostView>(i);
                    migrants.template copyParticle<useHostView>(migrantIndex, owned, i);

                    owned.template operator()<Particle_T::AttributeNames::ownershipState, true, useHostView>(i) = OwnershipState::dummy;
                  }
                });
              } else {
                deleteHaloParticles();

                utils::KokkosStorage<Particle_T> survivors {_ownedParticles.getLayout(), _ownedParticles.size()};

                Kokkos::View<int*, DeviceSpace> survivorCounter {"survivorCounter", 1};

                Kokkos::parallel_for("collectMigrantsAndReinsertOwned", Kokkos::RangePolicy<DeviceSpace::execution_space>(0, _ownedParticles.size()), KOKKOS_LAMBDA(int i) {

                  if (owned.template fulfillsIteratorRequirements<false, useHostView>(i, IteratorBehavior::owned, boxMinKokkos, boxMaxKokkos)) {
                    /* Particle is owned */

                    if (owned.template fulfillsIteratorRequirements<true, useHostView>(i, IteratorBehavior::owned, boxMinKokkos, boxMaxKokkos)) {
                      /* Particle is within the container */

                      // TODO: this can lead to heavy contention on the locks... -> maybe think of a better approach
                      int survivorIndex = Kokkos::atomic_fetch_inc(&survivorCounter(0));
                      survivors.template copyParticle<useHostView>(survivorIndex, owned, i);
                    } else {
                      /* Particle is outside the container boundary */

                      int migrantIndex = Kokkos::atomic_fetch_inc(&migrantCounter(0));
                      migrants.template copyParticle<useHostView>(migrantIndex, owned, i);
                    }
                  }
                });
                migrants.template markModified<DeviceSpace::execution_space>();
                survivors.template markModified<DeviceSpace::execution_space>();

                auto surviorCounterMirror = Kokkos::create_mirror_view(survivorCounter);
                Kokkos::deep_copy(surviorCounterMirror, survivorCounter);

                int numSurvivors = surviorCounterMirror(0);
                survivors.resize(numSurvivors);
                survivors.overrideSize(numSurvivors);

                _ownedParticles = survivors;
              }

              auto migrantCounterMirror = Kokkos::create_mirror_view(migrantCounter);
              Kokkos::deep_copy(migrantCounterMirror, migrantCounter);

              int numMigrants = migrantCounterMirror(0);
              migrants.resize(numMigrants);

              _ownedParticles.template markModified<DeviceSpace::execution_space>();
              _ownedParticles.markLayoutModified(owned.getLayout());

              migrants.template sync<HostSpace::execution_space>();
              migrants.syncSoAToAoS();
              std::vector<Particle_T> migrantVector {};
              migrantVector.reserve(numMigrants);

              for (int i = 0; i < numMigrants; ++i) {
                Particle_T migrant = migrants.getAoS().getParticle(i);
                migrantVector.push_back(migrant);
              }

              return migrantVector;
            }

            [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
                // TODO
                return TraversalSelectorInfo{};
            }

            /* TODO: Begin of Code Crime (investigate generalizations of this -> every container does EXACTLY the same) */
            [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
                IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                utils::optRef<typename ContainerIterator<Particle_T, true, false>::ParticleVecType> additionalVectors =
                  std::nullopt) override {
                // Copy from DirectSum.h
                convertToAoS();
                return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
                IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                utils::optRef<typename ContainerIterator<Particle_T, false, false>::ParticleVecType> additionalVectors =
                    std::nullopt) const override {
                // Copy from DirectSum.h
                // TODO: think about how to handle case when particles are stored in SoA Format, or mark this as discouraged (Problem: const qualifier prohibits necessary data conversions/copies...)
                return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
            }

            [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                utils::optRef<typename ContainerIterator<Particle_T, true, true>::ParticleVecType> additionalVectors =
                    std::nullopt) override {
                // Copy from DirectSum.h
                convertToAoS();
                return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                utils::optRef<typename ContainerIterator<Particle_T, false, true>::ParticleVecType> additionalVectors =
                    std::nullopt) const override {
                // Copy from DirectSum.h
                // TODO: think about how to handle case when particles are stored in SoA Format, or mark this as discouraged (Problem: const qualifier prohibits necessary data conversions/copies...)
                return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            /* TODO: End of Code Crime */

            template <typename Lambda>
            void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
              // TODO: decide if we really want to use this/if we can merge it with the Kokkos functions
            }

            template <typename Lambda>
            void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                    const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
              // TODO: decide if we really want to use this/if we can merge it with the Kokkos functions
            }

            template <class ExecSpace, typename Lambda>
            void forEachKokkos(Lambda& forEachLambda, IteratorBehavior behavior, const std::string& label = "forEachKokkos") {

              auto& boxMin = ParticleContainerInterface<Particle_T>::_boxMin;
              auto& boxMax = ParticleContainerInterface<Particle_T>::_boxMax;

              forEachInRegionKokkos<ExecSpace, false>(forEachLambda, behavior, boxMin, boxMax, label);
            }

              template <class ExecSpace, bool regionIter, typename Lambda>
            void forEachInRegionKokkos(Lambda forEachLambda, IteratorBehavior behavior, const std::array<double, 3> &lowerCorner, const std::array<double, 3>& higherCorner, const std::string& label = "forEachInRegionKokkos") {

              const auto& lowerCornerKokkos = Kokkos::Array{lowerCorner.at(0), lowerCorner.at(1), lowerCorner.at(2)};
              const auto& higherCornerKokkos = Kokkos::Array{higherCorner.at(0), higherCorner.at(1), higherCorner.at(2)};

              if (_dataLayout == DataLayoutOption::aos) {
                convertToAoS(); /* this guarantees syncing to the Host (AoS must not run on CUDA - so far) */
                _ownedParticles.markLayoutModified(DataLayoutOption::aos);
              }
              else if (_dataLayout == DataLayoutOption::soa) {
                convertToSoA();
                _ownedParticles.template sync<ExecSpace>();
                _haloParticles.template sync<ExecSpace>();
              }

              // TODO: make sure that no AoS runs on Cuda
              // TODO: decide on how to deduce this template parameters based on the ExecSpace
              constexpr bool host = false;

              /* owned or dummies; this basically filters out only halo, sequential, container only */
              if (behavior & 0b101) {
                auto& owned = _ownedParticles;
                Kokkos::parallel_for(label + "_owned", Kokkos::RangePolicy<ExecSpace>(0, _ownedParticles.size()), KOKKOS_LAMBDA(int i)  {

                  if (owned.template fulfillsIteratorRequirements<regionIter, host>(i, behavior, lowerCornerKokkos, higherCornerKokkos)) {
                    forEachLambda(i, owned);
                  }
                });
              }
              /* halo or dummies; this basically filters out only owned, sequential, container only */
              if (behavior & 0b110)  {
                auto& halo = _haloParticles;
                Kokkos::parallel_for(label + "_halo", Kokkos::RangePolicy<ExecSpace>(0, _haloParticles.size()), KOKKOS_LAMBDA(int i)  {
                  if (halo.template fulfillsIteratorRequirements<regionIter, host>(i, behavior, lowerCornerKokkos, higherCornerKokkos)) {
                    forEachLambda(i, halo);
                  }
                });
              }
              /* force sequential (sequential yes, but which particles? all? only owned? */
              if (behavior & 0b1000) {
                // TODO: implement
              }
              /* container only (whatever that means...) */
              if (behavior & 0b10000) {
                // TODO: implement
              }
              // TODO: also consider particles in the additionalStorage (they come from the LogicHandler's additional buffer)

              if (_dataLayout == DataLayoutOption::soa) {
                _ownedParticles.template markModified<ExecSpace>();
                _haloParticles.template markModified<ExecSpace>();
                _ownedParticles.markLayoutModified(DataLayoutOption::soa);
              }
            }

            template<class ExecSpace, typename Result, typename Reduction, typename Lambda>
            void reduceKokkos(Lambda reduceLambda, Result& result, IteratorBehavior behavior, const std::string& label = "reduceKokkos") {

              auto& boxMin = ParticleContainerInterface<Particle_T>::_boxMin;
              auto& boxMax = ParticleContainerInterface<Particle_T>::_boxMax;

              reduceInRegionKokkos<ExecSpace, false, Result, Reduction>(reduceLambda, result, behavior, boxMin, boxMax, label);
            }

            // TODO: declare that changes to particles is actually undefined behavior (as no memory syncs and data conversions are issued)
            template<class ExecSpace, bool regionIter, typename Result, typename Reduction, typename Lambda>
            void reduceInRegionKokkos(Lambda reduceLambda, Result& result, IteratorBehavior behavior, const std::array<double, 3>& lowerCorner, const std::array<double, 3>& higherCorner, const std::string& label = "reduceInRegionKokkos") {
              if (_dataLayout == DataLayoutOption::aos) {
                convertToAoS();
              }
              else if (_dataLayout == DataLayoutOption::soa) {
                convertToSoA();
                _ownedParticles.template sync<ExecSpace>();
                _haloParticles.template sync<ExecSpace>();
              }
              // TODO: make sure that no AoS runs on Cuda
              // TODO: decide on how to deduce this template parameters based on the ExecSpace
              constexpr bool host = false;

              const auto& lowerCornerKokkos = Kokkos::Array{lowerCorner.at(0), lowerCorner.at(1), lowerCorner.at(2)};
              const auto& higherCornerKokkos = Kokkos::Array{higherCorner.at(0), higherCorner.at(1), higherCorner.at(2)};

              /* owned or dummies; this basically filters out only halo, sequential, container only */
              if (behavior & 0b101) {
                auto& owned = _ownedParticles;
                Kokkos::parallel_reduce(label + "_owned", Kokkos::RangePolicy<ExecSpace>(0, _ownedParticles.size()), KOKKOS_LAMBDA(int i, Result& localResult)  {
                  if (owned.template fulfillsIteratorRequirements<regionIter, host>(i, behavior, lowerCornerKokkos, higherCornerKokkos)) {
                    reduceLambda(i, owned, localResult);
                  }
                }, Reduction(result));
              }
              /* halo or dummies; this basically filters out only owned, sequential, container only */
              if (behavior & 0b110) {
                auto& halo = _haloParticles;
                Kokkos::parallel_reduce(label + "_halo", Kokkos::RangePolicy<ExecSpace>(0, _haloParticles.size()), KOKKOS_LAMBDA(int i, Result& localResult)  {
                  if (halo.template fulfillsIteratorRequirements<regionIter, host>(i, behavior, lowerCornerKokkos, higherCornerKokkos)) {
                    reduceLambda(i, halo, localResult);
                  }
                }, Reduction(result));
              }
              /* force sequential (sequential yes, but which particles? all? only owned? */
              if (behavior & 0b1000) {
                // TODO: implement
              }
              /* container only (whatever that means...) */
              if (behavior & 0b10000) {
                // TODO: implement
              }
              // TODO: consider other behavior such as dummies, container only, ...

              // TODO: also consider particles in the additionalStorage (they come from the LogicHandler's additional buffer)
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
              if (cellIndex == 0 && particleIndex > _ownedParticles.size()) {
                // Particle does not exist in owned list
                return false;
              } else if (cellIndex == 1 && particleIndex > _haloParticles.size()) {
                // Particle does not exist in halo list
                return false;
              }
              auto& storage = (cellIndex == 0) ? _ownedParticles : _haloParticles;
              storage.template operator()<Particle_T::AttributeNames::ownershipState, true, true>(particleIndex) = OwnershipState::dummy;
              return true;
            }

protected:

  void addParticleImpl(const Particle_T &p) override {

    // TODO: check whether the target data layout is up to date

    std::lock_guard<AutoPasLock> guard (ownedLock);
    _ownedParticles.addParticle(p);
  }

  void addHaloParticleImpl(const Particle_T &haloP) override {

    // TODO: check whether the target layout is up to date

    std::lock_guard<AutoPasLock> guard(haloLock);
    _haloParticles.addParticle(haloP);
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

              kokkosDsTraversal->setOwnedToTraverse(_ownedParticles, targetLayout);
              kokkosDsTraversal->setHaloToTraverse(_haloParticles, targetLayout);
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
            _ownedParticles.markLayoutModified(targetLayout);
            kokkosDsTraversal->retrieveOwned(_ownedParticles);

            if (targetLayout == DataLayoutOption::aos) {
              if (_dataLayout != targetLayout) {
                convertToSoA();
              }
            }
            else if (targetLayout == DataLayoutOption::soa) {
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
          _ownedParticles.syncAoSToSoA();
          _haloParticles.syncAoSToSoA();
        }

        void convertToAoS() {
          _ownedParticles.syncSoAToAoS();
          _haloParticles.syncSoAToAoS();
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

          size_t sizeToCheck = (cellIndex == 0) ? _ownedParticles.size() : _haloParticles.size();
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
            while (particleIndex >= (cellIndex == 0 ? _ownedParticles.size() : _haloParticles.size())) {
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
  constexpr static bool useHostView = false;
#else
  using DeviceSpace = Kokkos::HostSpace;
  constexpr static bool useHostView = true;
#endif

  using HostSpace = Kokkos::HostSpace;

        mutable AutoPasLock ownedLock {};

        mutable AutoPasLock haloLock {};

        DataLayoutOption _dataLayout {};

        // TODO: find a better solution than having the flags in the container (maybe it makes sense to have them in the kokkosStorage
        // bool _aosUpToDate = false;

        // bool _soaUpToDate = false;

        utils::KokkosStorage<Particle_T> _ownedParticles {};

        utils::KokkosStorage<Particle_T> _haloParticles {};

        };
}