/**
 * @file SortedCellView.h
 * @date 18.01.2018
 * @author C. Menges
 */

#pragma once

#include <vector>
#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/CudaSoA.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class defines a sorted view on a given ParticleCell. Particles are sorted along the normalized vector r. The
 * projected position of a particle as well as a pointer to the underlying particle are stored in a vector, which is
 * sorted by the projected positions.
 * @note Insertion, deletion and change of position of particles invalidates the view.
 * @tparam Particle
 */
template <class Particle, class ParticleCellType>
class SortedCellView : public ParticleCell<Particle> {
 public:
  SortedCellView() = default;
  /**
   * Constructs a FullSortedParticleCell.
   * @param cell Cell whose particles are sorted.
   * @param r Normalized vector along particles are sorted.
   */
  SortedCellView(ParticleCellType &cell, const std::array<double, 3> &r) : _cell(&cell) { rebuild(r); }

  void addParticle(const Particle &p) override { _cell->addParticle(p); }

  /**
   * Rebuilds the internal sorted data structure.
   * @param r Normalized vector along particles are sorted.
   */
  void rebuild(const std::array<double, 3> &r) {
    _particles.clear();
    _particles.reserve(_cell->numParticles());
    for (auto p = getStaticCellIter(*_cell); p.isValid(); ++p) {
      _particles.push_back(std::make_pair(ArrayMath::dot(p->getR(), r), &(*p)));
    }
    std::sort(_particles.begin(), _particles.end(),
              [](const auto &a, const auto &b) -> bool { return a.first < b.first; });
  }

  SingleCellIteratorWrapper<Particle> begin() override { return _cell->begin(); }

  unsigned long numParticles() const override { return _particles.size(); }

  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particles.clear(); }

  void deleteByIndex(size_t index) override {
    if (index >= numParticles()) {
      AutoPasLog(error, "Index out of range");
      utils::ExceptionHandler::exception("Error: Index out of range");
    }

    if (index < numParticles() - 1) {
      std::swap(_particles[index], _particles[numParticles() - 1]);
    }
    _particles.pop_back();
  }

  void setCellLength(const std::array<double, 3> &cellLength) override { _cell->setCellLength(cellLength); }

  std::array<double, 3> getCellLength() const override { return _cell->getCellLength(); }

  /**
   * Type of the internal iterator.
   */
  typedef internal::SingleCellIterator<Particle, SortedCellView<Particle, ParticleCellType>> iterator_t;

  /**
   * Sorted vector of projected positions and particle pointers.
   */
  std::vector<std::pair<double, Particle *>> _particles;

  /**
   * Underlying cell.
   */
  ParticleCellType *_cell;
};
}  // namespace autopas
