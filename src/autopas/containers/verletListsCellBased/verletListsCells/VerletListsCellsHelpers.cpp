/**
 * @file VerletListsCellsHelpers.cpp
 * @author F.Gratl
 * @date 29.05.24
 */

#include "VerletListsCellsHelpers.h"

#include <algorithm>
#include <cmath>

#include "utils/ThreeDimensionalMapping.h"

namespace autopas::VerletListsCellsHelpers {

size_t estimateListLength(size_t numParticles, const std::array<double, 3> &boxSize, double interactionLength,
                          double correctionFactor) {
  const auto boxVolume = boxSize[0] * boxSize[1] * boxSize[2];
  constexpr double sphereConstant = 4. / 3. * M_PI;
  const auto listVolume = sphereConstant * interactionLength * interactionLength * interactionLength;
  const auto volumeFraction = listVolume / boxVolume;
  // ceil because if we need space for e.g. 3.5 particles better reserve for 4.
  return static_cast<size_t>(std::ceil(numParticles * volumeFraction * correctionFactor));
}

std::vector<BaseStepOffsets> buildC08BaseStep(const std::array<int, 3> &cellsPerDim) {
  // cellOffset1, cellOffset2, list estimation factor
  std::vector<BaseStepOffsets> offsets{};
  offsets.reserve(14);
  // This is currently guaranteed by the VerletListsCells constructor
  constexpr int interactionCellsPerDim = 1;
  // We estimate factors for fraction of particles needing a neighbor in the base cell.
  // Four cases: self (cell interacts with itself), and cells which share a face, edge, or corner with the base cell.
  // The factors are very conservative and express which portion of the cell could possibly have an interaction partner
  // in the other cell. Considerations for each factor are:
  // [0] self: Every particle could find a partner within the cell => 1
  // [1] face: Only case where a particle could not find a partner is if they are on the opposite cell edge,
  // which is unlikely => 1
  // [2] edge: Particles further away from the shared edge than r_i can't find a partner in the other cell.
  // The volume of the possible region is a quarter of a cylinder of height and radius r_i.
  // [3] corner: Particles further away from the shared corner than r_i can't find a partner in the other cell.
  // The volume of the possible region is an eighth of a sphere of radius r_i.
  constexpr std::array<double, 4> estimatorFactors{1., 1., 1. / 4. * M_PI, 1. / 6. * M_PI};
  // Idea: go over all 27 partner cells, create a pair of offsets from the base cell
  // and then shift the interaction into the 2x2 block in the positive axis directions.
  for (int z = -interactionCellsPerDim; z <= interactionCellsPerDim; ++z) {
    for (int y = -interactionCellsPerDim; y <= interactionCellsPerDim; ++y) {
      for (int x = -interactionCellsPerDim; x <= interactionCellsPerDim; ++x) {
        // Find the relative cell indices for the cell-cell interaction and shift them inside the 2x2x2 box
        int baseCell = 0;
        int partnerCell = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, cellsPerDim);
        if (x < 0) {
          baseCell -= x;
          partnerCell -= x;
        }
        if (y < 0) {
          baseCell -= y * cellsPerDim[0];
          partnerCell -= y * cellsPerDim[0];
        }
        if (z < 0) {
          baseCell -= z * cellsPerDim[0] * cellsPerDim[1];
          partnerCell -= z * cellsPerDim[0] * cellsPerDim[1];
        }
        // Count number of non-aligned dimensions
        const auto factor = estimatorFactors[std::abs(x) + std::abs(y) + std::abs(z)];
        const auto &[smallerIndex, biggerIndex] = std::minmax(baseCell, partnerCell);
        // Check if this offset tuple is already in the offsets list and if not add it.
        if (auto tuple = BaseStepOffsets{smallerIndex, biggerIndex, factor};
            std::find(offsets.begin(), offsets.end(), tuple) == offsets.end()) {
          offsets.emplace_back(tuple);
        }
      }
    }
  }
  // Sort offsets to group processing of the same cells (-> slightly better cache re-usage)
  std::sort(offsets.begin(), offsets.end(), [](const auto &pair1, const auto &pair2) {
    const auto &[a1, b1, f1] = pair1;
    const auto &[a2, b2, f2] = pair2;

    if (a1 == a2) {
      return b1 < b2;
    } else {
      return a1 < a2;
    }
  });
  return offsets;
}

bool BaseStepOffsets::operator==(const BaseStepOffsets &rhs) const {
  return offset1 == rhs.offset1 && offset2 == rhs.offset2 && listSizeEstimateFactor == rhs.listSizeEstimateFactor;
}

bool BaseStepOffsets::operator!=(const BaseStepOffsets &rhs) const { return !(rhs == *this); }

}  // namespace autopas::VerletListsCellsHelpers