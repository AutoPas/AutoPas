/**
 * @file VerletListsKokkos.h
 * @date 28.05.2026
 * @author Franziska Duhr
 */
#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/utils/KokkosAoS.h"
#include "autopas/utils/KokkosSoA.h"
#include "autopas/utils/KokkosStorage.h"
#include "traversals/VerletListsKokkosTraversalInterface.h"

namespace autopas {

/**
 * @class VerletListsKokkos
 * @brief A container for managing particles using Kokkos for parallel computation.
 * @tparam Particle_T AoS Particle type that is used with this container
 */
template <class Particle_T>
class VerletListsKokkos : public ParticleContainerInterface<Particle_T> {
    public:
    VerletListsKokkos(DataLayoutOption dataLayout, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double skin, const double cutoff)
        : ParticleContainerInterface<Particle_T>(boxMin, boxMax, skin), _dataLayout(dataLayout), _cutoff(cutoff) {
            if (dataLayout == DataLayoutOption::aos) {
                _aosUpToDate = true;
            }
            else if (dataLayout == DataLayoutOption::soa) {
                _soaUpToDate = true;
            }
            _ownedParticles.setLayout(dataLayout);
            _haloParticles.setLayout(dataLayout); 
    }
    
    [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletListsKokkos; }

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
            capacityOwned = capacityOwned == 0 ? 1 : capacityOwned * 2;
            _ownedParticles.resize(capacityOwned);
        }
        _ownedParticles.addParticle(numberOfOwned++, p);
        _neighborListValid = false;
    }

    void addHaloParticleImpl(const Particle_T &haloP) override {

        std::lock_guard<AutoPasLock> guard(haloLock);

        if (numberOfHalo >= capacityHalo) {
            capacityHalo = capacityHalo == 0 ? 1 : capacityHalo * 2;
            _haloParticles.resize(capacityHalo);
        }
        _haloParticles.addParticle(numberOfHalo++, haloP);
        _neighborListValid = false;
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
        _neighborListValid = false;
    }

    void deleteAllParticles() override {
        numberOfHalo = 0;
        numberOfOwned = 0;
        _neighborListValid = false;
    }

    size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
        size_t number = 0;
        if (behavior & 0b1) {
            number += numberOfOwned;
        }
        if (behavior & 0b10) {
            number += numberOfHalo;
        }
        // TODO: other behaviors (Actually, dummies can be somewhere in both(?) lists)
        // maybe find dummies in owned list with the help of reduceKokkos()
        return number;
    }

    size_t size() const override {
        return numberOfHalo + numberOfOwned;
    }

    // TODO: move to a faster (cell-binned) algorithm.
    void rebuildNeighborLists([[maybe_unused]]TraversalInterface *traversal) override {
        spdlog::debug("Rebuilding Verlet Lists. Number of owned particles: {}, number of halo particles: {}", numberOfOwned, numberOfHalo);
        if (_neighborListValid) {
            return;
        }
        spdlog::info("Rebuilding Verlet Lists with cutoff {} and skin {}", _cutoff, this->getVerletSkin());
        convertToAoS();

        const double interactionLength = _cutoff + this->getVerletSkin();
        const double interactionLengthSqr = interactionLength * interactionLength;

        Kokkos::resize(_neighborListOffsets, numberOfOwned + 1);
        _neighborListOffsets.clear_sync_state();
        _neighborListOffsets.modify_host();
        auto h_offsets = _neighborListOffsets.h_view;

        // Full neighbor list: every pair {i,j} is stored under both i and j (j != i).
        h_offsets(0) = 0;
        for (size_t i = 0; i < numberOfOwned; ++i) {
            const auto &ri = _ownedParticles.getAoS().getParticle(i).getR();
            size_t count = 0;
            for (size_t j = 0; j < numberOfOwned; ++j) {
                if (j == i) {
                    continue;
                }
                const auto &rj = _ownedParticles.getAoS().getParticle(j).getR();
                const double dx = ri[0] - rj[0];
                const double dy = ri[1] - rj[1];
                const double dz = ri[2] - rj[2];
                if (dx * dx + dy * dy + dz * dz < interactionLengthSqr) {
                    ++count;
                }
            }
            h_offsets(i + 1) = count;
        }
        for (size_t i = 0; i < numberOfOwned; ++i) {
            h_offsets(i + 1) += h_offsets(i);
        }

        const size_t totalNeighbors = numberOfOwned == 0 ? 0 : h_offsets(numberOfOwned);
        Kokkos::resize(_neighborListEntries, totalNeighbors);
        _neighborListEntries.clear_sync_state();
        _neighborListEntries.modify_host();
        auto h_entries = _neighborListEntries.h_view;

        for (size_t i = 0; i < numberOfOwned; ++i) {
            const auto &ri = _ownedParticles.getAoS().getParticle(i).getR();
            size_t slot = h_offsets(i);
            for (size_t j = 0; j < numberOfOwned; ++j) {
                if (j == i) {
                    continue;
                }
                const auto &rj = _ownedParticles.getAoS().getParticle(j).getR();
                const double dx = ri[0] - rj[0];
                const double dy = ri[1] - rj[1];
                const double dz = ri[2] - rj[2];
                if (dx * dx + dy * dy + dz * dz < interactionLengthSqr) {
                    h_entries(slot++) = j;
                }
            }
        }

        _neighborListOffsets.template sync<typename DeviceSpace::execution_space>();
        _neighborListEntries.template sync<typename DeviceSpace::execution_space>();

        // owned-halo neighbor list: for each owned particle, the halo particles within range.
        // Halo particles are ghosts and never receive force, so this list is one-directional
        // (each owned i gathers from its halo neighbors) and needs no newton3.
        Kokkos::resize(_haloNeighborListOffsets, numberOfOwned + 1);
        _haloNeighborListOffsets.clear_sync_state();
        _haloNeighborListOffsets.modify_host();
        auto h_haloOffsets = _haloNeighborListOffsets.h_view;

        h_haloOffsets(0) = 0;
        for (size_t i = 0; i < numberOfOwned; ++i) {
            const auto &ri = _ownedParticles.getAoS().getParticle(i).getR();
            size_t count = 0;
            for (size_t j = 0; j < numberOfHalo; ++j) {
                const auto &rj = _haloParticles.getAoS().getParticle(j).getR();
                const double dx = ri[0] - rj[0];
                const double dy = ri[1] - rj[1];
                const double dz = ri[2] - rj[2];
                if (dx * dx + dy * dy + dz * dz < interactionLengthSqr) {
                    ++count;
                }
            }
            h_haloOffsets(i + 1) = count;
        }
        for (size_t i = 0; i < numberOfOwned; ++i) {
            h_haloOffsets(i + 1) += h_haloOffsets(i);
        }

        const size_t totalHaloNeighbors = numberOfOwned == 0 ? 0 : h_haloOffsets(numberOfOwned);
        Kokkos::resize(_haloNeighborListEntries, totalHaloNeighbors);
        _haloNeighborListEntries.clear_sync_state();
        _haloNeighborListEntries.modify_host();
        auto h_haloEntries = _haloNeighborListEntries.h_view;

        for (size_t i = 0; i < numberOfOwned; ++i) {
            const auto &ri = _ownedParticles.getAoS().getParticle(i).getR();
            size_t slot = h_haloOffsets(i);
            for (size_t j = 0; j < numberOfHalo; ++j) {
                const auto &rj = _haloParticles.getAoS().getParticle(j).getR();
                const double dx = ri[0] - rj[0];
                const double dy = ri[1] - rj[1];
                const double dz = ri[2] - rj[2];
                if (dx * dx + dy * dy + dz * dz < interactionLengthSqr) {
                    h_haloEntries(slot++) = j;
                }
            }
        }

        _haloNeighborListOffsets.template sync<typename DeviceSpace::execution_space>();
        _haloNeighborListEntries.template sync<typename DeviceSpace::execution_space>();

        _neighborListValid = true;
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
            // Caller asserts positions haven't moved enough to invalidate the list.
            // TODO: collect particles and mark non owned as dummy
            // i.e.: those particles which are outside of the box shall be returned, halo particles should be marked as dummies
            return {};
        }
        deleteHaloParticles();
        _neighborListValid = false;
        // TODO: delete dummy particles
        // TODO: determine those particles which are not inside of the box and return them

        return {};
    }

    [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
        using namespace autopas::utils::ArrayMath::literals;
        // No spatial decomposition (yet) - treat the whole box as a single cell, like DirectSum.
        return TraversalSelectorInfo({1, 1, 1}, getInteractionLength(),
                                     this->getBoxMax() - this->getBoxMin(), 0);
    }

    /* TODO: Begin of Code Crime (investigate generalizations of this -> every container does EXACTLY the same) */
    [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
        IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
        utils::optRef<typename ContainerIterator<Particle_T, true, false>::ParticleVecType> additionalVectors =
            std::nullopt) override {
        // Copy from DirectSum.h
        convertToAoS();
        _soaUpToDate = false;
        return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
    }

    [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
        IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
        utils::optRef<typename ContainerIterator<Particle_T, false, false>::ParticleVecType> additionalVectors =
            std::nullopt) const override {
        // Copy from DirectSum.h
        // TODO: think about how to handle case when particles are stored in SoA Format, or mark this as discouraged
        return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
    }

    [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
        const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
        utils::optRef<typename ContainerIterator<Particle_T, true, true>::ParticleVecType> additionalVectors =
            std::nullopt) override {
        // Copy from DirectSum.h
        convertToAoS();
        _soaUpToDate = false;
        return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
    }

    [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
        const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
        utils::optRef<typename ContainerIterator<Particle_T, false, true>::ParticleVecType> additionalVectors =
            std::nullopt) const override {
        // Copy from DirectSum.h
        // TODO: think about how to handle case when particles are stored in SoA Format, or mark this as discouraged (Problem: const qualifier prohibits data conversions...)
        return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
    }

    /* TODO: End of Code Crime */

    template <typename Lambda>
    void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
            const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
        // TODO
    }

    template <typename Lambda>
    void forEachInRegionKokkos(Lambda forEachLambda, const std::array<double, 3> &lowerCorner, const std::array<double, 3>& higherCorner, IteratorBehavior behavior) {
        // TODO: this will still require to pass the all particle's indices to the lambda clause with a preceeding check if the current particle is within the requested region in addition to the test for the ownershipstate
    }


    template <typename Lambda>
    void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
        // TODO
    }

    template <class ExecSpace, typename Lambda>
    void forEachKokkos(Lambda forEachLambda, IteratorBehavior behavior) {

        convertToSoA();
        _ownedParticles.template sync<ExecSpace>();
        _haloParticles.template sync<ExecSpace>();
        _aosUpToDate = false;

        // TODO: make sure that no AoS runs on Cuda
        if (behavior & 0b1) {
        if (numberOfOwned > 0 and not _ownedParticles.deviceViewsAllocated()) {
            utils::ExceptionHandler::exception("VerletListsKokkos::forEachKokkos: owned device views are not allocated.");
        }
        const auto owned = _ownedParticles.deviceView();
        const auto lambda = forEachLambda;
        Kokkos::parallel_for("forEachKokkosOwned", Kokkos::RangePolicy<ExecSpace>(0, numberOfOwned), KOKKOS_LAMBDA(int i)  {
            // TODO: here, also a dummy check is required
            lambda(i, owned);
        });
        }
        if (behavior & 0b10) {
        if (numberOfHalo > 0 and not _haloParticles.deviceViewsAllocated()) {
            utils::ExceptionHandler::exception("VerletListsKokkos::forEachKokkos: halo device views are not allocated.");
        }
        const auto halo = _haloParticles.deviceView();
        const auto lambda = forEachLambda;
        Kokkos::parallel_for("forEachKokkosHalo", Kokkos::RangePolicy<ExecSpace>(0, numberOfHalo), KOKKOS_LAMBDA(int i)  {
            // TODO: here, also a dummy check is required
            lambda(i, halo);
        });
        }
        // TODO: consider other behavior such as dummies, container only, ...

        _ownedParticles.template markModified<ExecSpace>();
        _haloParticles.template markModified<ExecSpace>();
        _soaUpToDate = true;
        _aosUpToDate = false;
    }

    template<class ExecSpace, typename Result, typename Reduction, typename Lambda>
    void reduceKokkos(Lambda reduceLambda, Result& result, IteratorBehavior behavior) {

        convertToSoA();
        _ownedParticles.template sync<ExecSpace>();
        _haloParticles.template sync<ExecSpace>();

        if (behavior & 0b1) {
        if (numberOfOwned > 0 and not _ownedParticles.deviceViewsAllocated()) {
            utils::ExceptionHandler::exception("VerletListsKokkos::reduceKokkos: owned device views are not allocated.");
        }
        const auto owned = _ownedParticles.deviceView();
        const auto lambda = reduceLambda;
        Kokkos::parallel_reduce("reduceKokkosOwned", Kokkos::RangePolicy<ExecSpace>(0, numberOfOwned), KOKKOS_LAMBDA(int i, Result& localResult)  {
            lambda(i, owned, localResult);
        }, Reduction(result));
        }
        if (behavior & 0b10) {
        if (numberOfHalo > 0 and not _haloParticles.deviceViewsAllocated()) {
            utils::ExceptionHandler::exception("VerletListsKokkos::reduceKokkos: halo device views are not allocated.");
        }
        const auto halo = _haloParticles.deviceView();
        const auto lambda = reduceLambda;
        Kokkos::parallel_reduce("reduceKokkosHalo", Kokkos::RangePolicy<ExecSpace>(0, numberOfHalo), KOKKOS_LAMBDA(int i, Result& localResult)  {
            lambda(i, halo, localResult);
        }, Reduction(result));
        }
        // TODO: consider other behavior such as dummies, container only, ...
    }

    [[nodiscard]] double getCutoff() const final { return _cutoff; }

    void setCutoff(double cutoff) final { _cutoff = cutoff; }

    [[nodiscard]] double getInteractionLength() const final { return _cutoff + this->getVerletSkin(); }

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
        if (cellIndex == 0 && particleIndex >= numberOfOwned) {
            // Particle does not exist in owned list
            return false;
        } else if (cellIndex == 1 && particleIndex >= numberOfHalo) {
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
            auto* kokkosVLTraversal = dynamic_cast<VerletListsKokkosTraversalInterface<Particle_T> *>(traversal);
            if (kokkosVLTraversal) {

              auto targetLayout = traversal->getDataLayout();

              if (targetLayout != _dataLayout) {
                convertTo(targetLayout);
              }

              kokkosVLTraversal->setOwnedToTraverse(_ownedParticles);
              kokkosVLTraversal->setHaloToTraverse(_haloParticles);
              kokkosVLTraversal->setNeighborList(_neighborListOffsets.d_view, _neighborListEntries.d_view);
              kokkosVLTraversal->setHaloNeighborList(_haloNeighborListOffsets.d_view, _haloNeighborListEntries.d_view);
            }
            else {
                utils::ExceptionHandler::exception("The selected traversal is not compatible with the VerletListsKokkos container.");
            }
        }

        template <typename Traversal>
        void finishTraversal(Traversal traversal) {
          auto* kokkosVLTraversal = dynamic_cast<VerletListsKokkosTraversalInterface<Particle_T> *>(traversal);

          if (kokkosVLTraversal) {

            auto targetLayout = traversal->getDataLayout();

            kokkosVLTraversal->retrieveOwned(_ownedParticles);

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
            utils::ExceptionHandler::exception("The selected traversal is not compatible with the VerletListsKokkos container.");
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
#elif defined(KOKKOS_ENABLE_HIP)
    using DeviceSpace = Kokkos::HIPSpace;
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

    double _cutoff {0.0};

    // h_view is written by rebuildNeighborLists; d_view is read by the traversal.
    // owned-owned neighbor list (entries index into the owned particles)
    Kokkos::DualView<size_t*> _neighborListOffsets {"vl_neighborListOffsets", 0};
    Kokkos::DualView<size_t*> _neighborListEntries {"vl_neighborListEntries", 0};
    // owned-halo neighbor list (offsets per owned particle, entries index into the halo particles)
    Kokkos::DualView<size_t*> _haloNeighborListOffsets {"vl_haloNeighborListOffsets", 0};
    Kokkos::DualView<size_t*> _haloNeighborListEntries {"vl_haloNeighborListEntries", 0};
    bool _neighborListValid {false};

};

} 
