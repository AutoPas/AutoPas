/**
 * @file SortedCellView.h
 * @date 18.01.2018
 * @author C. Menges
 */

#pragma once

#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class defines a sorted view on a given ParticleCell. Particles are sorted along the normalized vector r.
 * \image html SortingPrinciple.png "Projection of particles in 2D"
 * The projected position of a particle as well as a pointer to the underlying particle are stored in a vector, which is
 * sorted by the projected positions.
 * @note Insertion, deletion and change of position of particles invalidates the view.
 * @tparam Particle
 */
template <class Particle, class ParticleCellType>
class SortedCellView : public ParticleCell<Particle> {
 public:
  /**
   * Type that holds or refers to the actual particles.
   */
  using StorageType = std::vector<Particle *>;
  /**
   * Constructs a FullSortedParticleCell.
   * @param cell Cell whose particles are sorted.
   * @param r Normalized vector along particles are sorted.
   */
  SortedCellView(ParticleCellType &cell, const std::array<double, 3> &r) : _cell(&cell) {
    _particles.reserve(cell.numParticles());
    for (auto &p : cell) {
      _particles.push_back(std::make_pair(utils::ArrayMath::dot(p.getR(), r), &p));
    }
    std::sort(_particles.begin(), _particles.end(),
              [](const auto &a, const auto &b) -> bool { return a.first < b.first; });
  }

  CellType getParticleCellTypeAsEnum() override { return CellType::SortedCellView; }

  /**
   * @copydoc ParticleCell::addParticle()
   */
  void addParticle(const Particle &p) override {}

  /**
   * @copydoc autopas::FullParticleCell::begin()
   * @note: iterating via begin / iterators does not yield the sorted order.
   */
  CellIterator<typename ParticleCellType::StorageType, true> begin() {
    return CellIterator<typename ParticleCellType::StorageType, true>(_cell->begin());
  }

  /**
   * @copydoc autopas::CellIterator<typename ParticleCellType::StorageType, false>::begin()
   * @note const version
   */
  CellIterator<typename ParticleCellType::StorageType, false> begin() const {
    return CellIterator<typename ParticleCellType::StorageType, false>(_cell->begin());
  }

  /**
   * @copydoc autopas::FullParticleCell::end()
   */
  CellIterator<typename ParticleCellType::StorageType, true> end() {
    return CellIterator<typename ParticleCellType::StorageType, true>(_cell->end());
  }

  /**
   * @copydoc autopas::FullParticleCell::end()
   * @note const version
   */
  CellIterator<typename ParticleCellType::StorageType, false> end() const {
    return CellIterator<typename ParticleCellType::StorageType, false>(_cell->end());
  }

  unsigned long numParticles() const override { return _particles.size(); }

  bool isEmpty() const override { return numParticles() == 0; }

  void clear() override { _particles.clear(); }

  void deleteDummyParticles() override {
    _particles.erase(std::remove_if(_particles.begin(), _particles.end(),
                                    [](const auto &particlePosPair) { return particlePosPair.second->isDummy(); }),
                     _particles.end());
  }

  void deleteByIndex(size_t index) override {
    if (index >= numParticles()) {
      AutoPasLog(ERROR, "Index out of range");
      utils::ExceptionHandler::exception("Error: Index out of range");
    }

    if (index < numParticles() - 1) {
      std::swap(_particles[index], _particles[numParticles() - 1]);
    }
    _particles.pop_back();
  }

  void setCellLength(std::array<double, 3> &cellLength) override { _cell->setCellLength(cellLength); }

  std::array<double, 3> getCellLength() const override { return _cell->getCellLength(); }

  /**
   * Returns the particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  Particle &at(size_t index) { return _particles.at(index); }

  /**
   * Returns the const particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  const Particle &at(size_t index) const { return _particles.at(index); }

  /**
   * Sorted vector of projected positions and particle pointers.
   */
  std::vector<std::pair<double, Particle *>> _particles;

  /**
   * Underlying cell.
   */
  ParticleCellType *_cell;

 private:
};
}  // namespace autopas
