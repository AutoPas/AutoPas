/**
 * @file LCC01Traversal3B.h
 * @author muehlhaeusser
 * @date 05.09.2023
 */

#pragma once

#include "LCTriTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c01 traversal for 3-body interactions.
 *
 * The traversal uses a c01 base step performed on every single cell.
 * This includes every unique combination of cells that contain particle triplets within the cutoff.
 *
 * newton3 cannot be applied!
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor type that defines the interaction of three particles.
 */
template <class ParticleCell, class Functor>
class LCC01Traversal3B : public C01BasedTraversal<ParticleCell, Functor, 3>,
                         public LCTriTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param functor The functor that defines the interaction of three particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC01Traversal3B(const std::array<unsigned long, 3> &dims, Functor *functor, const double interactionLength,
                            const std::array<double, 3> &cellLength, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        C01BasedTraversal<ParticleCell, Functor, 3>(dims, functor, interactionLength, cellLength, dataLayout,
                                                    useNewton3),
        _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/,
                     dataLayout, useNewton3),
        _functor(functor) {
    computeOffsets();
  }

  /**
   * Computes all cell triplets used in processBaseCell()
   */
  void computeOffsets();

  void traverseParticleTriplets() override;

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   *
   * This is because the cell functor in the c01 traversal is hardcoded to not allow newton 3 even if only one thread is
   * used.
   *
   * @return
   */
  [[nodiscard]] bool isApplicable() const override { return not this->_useNewton3; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c01_3b; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellFunctor.setSortingThreshold(sortingThreshold); }

 private:
  /**
   * Computes all interactions between the base
   * cell and adjacent cells.
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  inline void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  /**
   * Appends all needed Attributes to the SoA buffer in cell.
   * @tparam I
   * @param cell
   * @param appendCell
   */
  template <std::size_t... I>
  inline constexpr void appendNeeded(ParticleCell &cell, ParticleCell &appendCell, std::index_sequence<I...>) {
    cell._particleSoABuffer.template append<std::get<I>(Functor::getNeededAttr(std::false_type()))...>(
        appendCell._particleSoABuffer);
  }

  /**
   * Combination of triplets for processBaseCell().
   */
  std::vector<std::tuple<long, long, std::array<double, 3>>> _cellOffsets;

  /**
   * CellFunctor to be used for the traversal defining the interaction between three cells.
   */
  internal::CellFunctor3B<ParticleCell, Functor, false> _cellFunctor;

  Functor *_functor;
};

template <class ParticleCell, class Functor>
inline void LCC01Traversal3B<ParticleCell, Functor>::computeOffsets() {
  using namespace utils::ArrayMath::literals;

  // Helper function to get minimal distance between two cells
  auto cellDistance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    return std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * this->_cellLength[0],
                                 std::max(0l, (std::abs(y1 - y2) - 1l)) * this->_cellLength[1],
                                 std::max(0l, (std::abs(z1 - z2) - 1l)) * this->_cellLength[2]};
  };

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);
  _cellOffsets.emplace_back(0, 0, std::array<double, 3>{1., 1., 1.});

  // offsets for the first cell
  for (long x1 = -this->_overlap[0]; x1 <= static_cast<long>(this->_overlap[0]); ++x1) {
    for (long y1 = -this->_overlap[1]; y1 <= static_cast<long>(this->_overlap[1]); ++y1) {
      for (long z1 = -this->_overlap[2]; z1 <= static_cast<long>(this->_overlap[2]); ++z1) {
        // check distance between base cell and cell 1
        const auto dist01 = cellDistance(0l, 0l, 0l, x1, y1, z1);

        const double distSquare = utils::ArrayMath::dot(dist01, dist01);
        if (distSquare > interactionLengthSquare) continue;

        // offsets for cell 2
        for (long x2 = -this->_overlap[0]; x2 <= static_cast<long>(this->_overlap[0]); ++x2) {
          for (long y2 = -this->_overlap[1]; y2 <= static_cast<long>(this->_overlap[1]); ++y2) {
            for (long z2 = -this->_overlap[2]; z2 <= static_cast<long>(this->_overlap[2]); ++z2) {
              // check distance between cell 1 and cell 2
              const auto dist12 = cellDistance(x1, y1, z1, x2, y2, z2);

              const double dist12Squared = utils::ArrayMath::dot(dist12, dist12);
              if (dist12Squared > interactionLengthSquare) continue;

              // check distance between base cell and cell 2
              const auto dist02 = cellDistance(0l, 0l, 0l, x2, y2, z2);

              const double dist02Squared = utils::ArrayMath::dot(dist02, dist02);
              if (dist02Squared > interactionLengthSquare) continue;

              const long offset1 = utils::ThreeDimensionalMapping::threeToOneD(
                  x1, y1, z1, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

              const long offset2 = utils::ThreeDimensionalMapping::threeToOneD(
                  x2, y2, z2, utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension));

              if (offset2 <= offset1) continue;
              // sorting direction towards middle of cell 1 and cell 2
              const std::array<double, 3> sortDirection = {
                  (x1 + x2) * this->_cellLength[0], (y1 + y2) * this->_cellLength[1], (z1 + z2) * this->_cellLength[2]};
              _cellOffsets.emplace_back(offset1, offset2, utils::ArrayMath::normalize(sortDirection));
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class Functor>
inline void LCC01Traversal3B<ParticleCell, Functor>::processBaseCell(std::vector<ParticleCell> &cells, unsigned long x,
                                                                     unsigned long y, unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];

  for (auto const &[offset1, offset2, r] : _cellOffsets) {
    const unsigned long otherIndex1 = baseIndex + offset1;
    const unsigned long otherIndex2 = baseIndex + offset2;
    ParticleCell &otherCell1 = cells[otherIndex1];
    ParticleCell &otherCell2 = cells[otherIndex2];

    if (baseIndex == otherIndex1 && baseIndex == otherIndex2) {
      this->_cellFunctor.processCell(baseCell);
    } else if (baseIndex == otherIndex1 && baseIndex != otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell2);
    } else if (baseIndex != otherIndex1 && baseIndex == otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell1);
    } else if (baseIndex != otherIndex1 && otherIndex1 == otherIndex2) {
      this->_cellFunctor.processCellPair(baseCell, otherCell1);
    } else {
      this->_cellFunctor.processCellTriple(baseCell, otherCell1, otherCell2);
    }
  }
}

template <class ParticleCell, class Functor>
inline void LCC01Traversal3B<ParticleCell, Functor>::traverseParticleTriplets() {
  auto &cells = *(this->_cells);
  if (not this->isApplicable()) {
    utils::ExceptionHandler::exception("The C01 traversal cannot work with enabled newton3!");
  }

  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
