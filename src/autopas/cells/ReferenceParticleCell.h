/**
 * @file FullParticleCell.h
 * @date 18.01.2018
 * @author seckler
 */

#pragma once

#include <mutex>
#include <vector>

#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/utils/CudaSoA.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
    template <class Particle, class SoAArraysType = typename Particle::SoAArraysType>
    class ReferenceParticleCell : public ParticleCell<Particle> {
    public:
        /**
         * Constructs a new FullParticleCell.
         */
        ReferenceParticleCell()
                : _cellLength({std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max()}) {}

        /**
         * Constructs a new FullParticleCell with the given cell side length.
         * @param cellLength cell side length
         */
        explicit ReferenceParticleCell(const std::array<double, 3> &cellLength) : _cellLength(cellLength) {}

        /**
         * @copydoc ParticleCell::addParticle()
         */
        void addParticle(const Particle &p) override {
            particlesLock.lock();
            _particles.emplace_back(std::make_unique<Particle>(p));
            particlesLock.unlock();
        }

        SingleCellIteratorWrapper<Particle, true> begin() override {
            return SingleCellIteratorWrapper<Particle, true>(new iterator_t(this));
        }

        SingleCellIteratorWrapper<Particle, false> begin() const override {
            return SingleCellIteratorWrapper<Particle, false>(new const_iterator_t(this));
        }

        unsigned long numParticles() const override { return _particles.size(); }

        /**
         * Returns a reference to the element at position n in the cell.
         * @param n Position of an element in the container
         * @return Reference to the element
         */
        Particle &operator[](size_t n) { return *(_particles[n]); }

        /**
         * Returns a const reference to the element at position n in the cell.
         * @param n Position of an element in the container
         * @return Reference to the element
         */
        const Particle &operator[](size_t n) const { return *(_particles[n]); }

        /**
         * Returns the particle at position index. Needed by SingleCellIterator.
         * @param index the position of the particle to return.
         * @return the particle at position index.
         */
        Particle &at(size_t index) { return *(_particles.at(index)); }

        /**
         * Returns the const particle at position index. Needed by SingleCellIterator.
         * @param index the position of the particle to return.
         * @return the particle at position index.
         */
        const Particle &at(size_t index) const { return *(_particles.at(index)); }

        bool isNotEmpty() const override { return numParticles() > 0; }

        void clear() override { _particles.clear(); }

        void deleteByIndex(size_t index) override {
            std::lock_guard<AutoPasLock> lock(particlesLock);
            if (index >= numParticles()) {
                utils::ExceptionHandler::exception("Index out of range (range: [0, {}[, index: {})", numParticles(), index);
            }

            if (index < numParticles() - 1) {
                _particles[index].swap(_particles[numParticles() - 1]);
            }
            _particles.pop_back();
        }

        void setCellLength(std::array<double, 3> &cellLength) override { _cellLength = cellLength; }

        std::array<double, 3> getCellLength() const override { return _cellLength; }

        /**
         * Resizes the container so that it contains n elements.
         * @param n New container size
         * @param toInsert Particle to insert. This is needed to allow for non-default-constructible particles.
         */
        void resize(size_t n, const Particle &toInsert) { _particles.resize(n, std::unique_ptr<Particle>(toInsert)); }

        /**
         * Sort the particles in the cell by a dimension.
         * @param dim dimension to sort
         */
        void sortByDim(const size_t dim) {
            std::sort(_particles.begin(), _particles.end(),
                      [dim](const std::unique_ptr<Particle> a, const std::unique_ptr<Particle> b) -> bool { return a->getR()[dim] < b->getR()[dim]; });
        }

        /**
         * Requests that the vector capacity be at least enough to contain n elements.
         * @param n Minimum capacity for the vector.
         */
        void reserve(size_t n) { _particles.reserve(n); }

        /**
         * Storage of the molecules of the cell.
         */
        std::vector<std::unique_ptr<Particle>> _particles;

        /**
         * SoA buffer of this cell.
         */
        SoA<SoAArraysType> _particleSoABuffer;

        /**
         * Device particle SoABuffer.
         */
        CudaSoA<typename Particle::CudaDeviceArraysType> _particleSoABufferDevice;

        /**
         * Type of the internal iterator.
         */
        using iterator_t = internal::SingleCellIterator<Particle, ReferenceParticleCell<Particle, SoAArraysType>, true>;

        /**
         * Type of the internal const iterator.
         */
        using const_iterator_t = internal::SingleCellIterator<Particle, ReferenceParticleCell<Particle, SoAArraysType>, false>;

    private:
        AutoPasLock particlesLock;
        std::array<double, 3> _cellLength;
    };
}  // namespace autopas
