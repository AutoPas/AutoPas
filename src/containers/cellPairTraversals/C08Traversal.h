/**
 * @file C08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include <utils/WrapOpenMP.h>
#include "C08BasedTraversal.h"

namespace autopas {

/**
 * This class provides the c08 traversal.
 *
 * The traversal uses the c08 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eight colors is applied.
 *
 * @tparam ParticleCell the type of cells
 * @tparam CellFunctor the cell functor that defines the interaction of the
 * particles of two specific cells
 */
template <class ParticleCell, class CellFunctor>
class C08Traversal : public C08BasedTraversal<ParticleCell, CellFunctor> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param cellfunctor The cell functor that defines the interaction of
   * particles between two different cells.
   */
  explicit C08Traversal(const std::array<unsigned long, 3> &dims, CellFunctor *cellfunctor)
      : C08BasedTraversal<ParticleCell, CellFunctor>(dims, cellfunctor) {}
  // documentation in base class
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;
  bool isApplicable() override;
};

template <class ParticleCell, class CellFunctor>
inline bool C08Traversal<ParticleCell, CellFunctor>::isApplicable() {
  return true;
}

template <class ParticleCell, class CellFunctor>
inline void C08Traversal<ParticleCell, CellFunctor>::traverseCellPairs(std::vector<ParticleCell> &cells) {
  using std::array;
  const array<unsigned long, 3> stride = {2, 2, 2};
  array<unsigned long, 3> end;
  for (int d = 0; d < 3; ++d) {
    end[d] = this->_cellsPerDimension[d] - 1;
  }

#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
    for (unsigned long col = 0; col < 8; ++col) {
      std::array<unsigned long, 3> start = ThreeDimensionalMapping::oneToThreeD(col, stride);

      // intel compiler demands following:
      const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
      const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
      const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3)
#endif
      for (unsigned long z = start_z; z < end_z; z += stride_z) {
        for (unsigned long y = start_y; y < end_y; y += stride_y) {
          for (unsigned long x = start_x; x < end_x; x += stride_x) {
            unsigned long baseIndex = ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
            this->processBaseCell(cells, baseIndex);
          }
        }
      }
    }
  }
};

}  // namespace autopas