/**
 * @file KokkosDirectSum.h
 * @date 31 October 2025
 * @author Luis Gall
 */

// This file is oriented on the implementation of the CellBasedParticleContainer

#pragma once

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/TraversalInterface.h"

#include <Kokkos_Core.hpp>

namespace autopas {

template <class Particle_T>
    class KokkosDirectSum : public ParticleContainerInterface<Particle_T> {

        public:
            KokkosDirectSum(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double skin)
                : ParticleContainerInterface<Particle_T>(boxMin, boxMax, skin)
                {}

            [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::kokkosDirectSum; }

            void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
                // TODO
            }

            void addParticleImpl(const Particle_T &p) override {
                // TODO
            }

            void addHaloParticleImpl(const Particle_T &haloP) override {
                // TODO
            }

            bool updateHaloParticle(const Particle_T &haloParticle) override {
                // TODO

                return false;
            }

            void deleteHaloParticles() override {
                // TODO
            }

            void deleteAllParticles() override {
                // TODO
            }

            size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const override {
                // TODO
                return 0;
            }

            size_t size() const override {
                // TODO
                return 0;
            }

            void rebuildNeighborLists(TraversalInterface *traversal) override {
                // No-Op
            }

            void computeInteractions(TraversalInterface *traversal) override {
                // TODO
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
                return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
                const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
                typename ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors =
                    nullptr) const override {
                return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
            }

            /* TODO: End of Code Crime */

            template <typename Lambda>
            void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                    const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
                // TODO
            }

            [[nodiscard]] double getCutoff() const final { return 0; }

            void setCutoff(double cutoff) final { /* No-Op */ };

            [[nodiscard]] double getInteractionLength() const final { return 0; }

            std::tuple<const Particle_T*, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                    IteratorBehavior iteratorBehavior) const final {
                // TODO
                return {};
            }

            std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                                     IteratorBehavior iteratorBehavior,
                                                                     const std::array<double, 3> &boxMin,
                                                                     const std::array<double, 3> &boxMax) const final {
                // TODO
                return {};
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

    #if defined(Kokkos_ENABLE_CUDA)
        using MemSpace = Kokkos::CudaSpace;

    #elif defined(Kokkos_ENABLE_OPENMP)
        using MemSpace = Kokkos::OpenMP;
    #endif
        
    
        /* AoS Data structure */
        Kokkos::View<Particle_T*> _ownedParticles {};

        Kokkos::View<Particle_T*> _haloParticles {};

        /* SoA data structure */

        // TODO: implement

        /* AoSoA data structure */

        // TODO: implement

        };
}