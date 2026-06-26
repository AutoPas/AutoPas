/**
 * @file TwoCellsInteractionHitrateGeneratorTest.cpp
 *
 * @date 26.06.2026
 * @author hmeyran
 */

#include "TwoCellsInteractionHitrateGeneratorTest.h"

#include <array>
#include <cmath>
#include <set>
#include <vector>

#include "autopas/utils/ExceptionHandler.h"
#include "autopasTools/PseudoContainer.h"
#include "autopasTools/generators/TwoCellsInteractionHitrateGenerator.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

using MoleculeType = mdLib::MoleculeLJ;
using CellType = std::vector<MoleculeType>;

/**
 * Cutoff radius used across all tests.
 */
constexpr double kCutoff = 3.0;
/**
 * Shared x-coordinate of the boundary plane between the two cells.
 */
constexpr double kBoundary = 10.0;
/**
 * x-extent of each cell; must exceed 1.2 * kCutoff so the far zone fits.
 */
constexpr double kExtent = 8.0;

/**
 * Bounding box of cell 1: occupies [kBoundary - kExtent, kBoundary] in x.
 */
const std::array<double, 3> kBoxMin1 = {kBoundary - kExtent, 0., 0.};
const std::array<double, 3> kBoxMax1 = {kBoundary, 8., 8.};
/**
 * Bounding box of cell 2: occupies [kBoundary, kBoundary + kExtent] in x.
 */
const std::array<double, 3> kBoxMin2 = {kBoundary, 0., 0.};
const std::array<double, 3> kBoxMax2 = {kBoundary + kExtent, 8., 8.};

/**
 * Computes the Euclidean distance between two particles.
 * @param a First particle.
 * @param b Second particle.
 * @return Distance between a and b.
 */
double dist3D(const MoleculeType &a, const MoleculeType &b) {
  double d = 0.;
  for (int dim = 0; dim < 3; ++dim) {
    const double delta = a.getR()[dim] - b.getR()[dim];
    d += delta * delta;
  }
  return std::sqrt(d);
}

/**
 * Fills two cells with n particles at the requested hitrate using the module-level box definitions.
 * @param v1 Target vector for cell 1 particles.
 * @param v2 Target vector for cell 2 particles.
 * @param n Number of particles per cell.
 * @param hitrate Fraction of particle pairs that are within cutoff distance.
 * @param seed RNG seed for deterministic placement.
 */
void fill(CellType &v1, CellType &v2, std::size_t n, double hitrate, unsigned int seed = 42) {
  autopasTools::PseudoContainer c1(v1), c2(v2);
  autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      c1, c2, kBoxMin1, kBoxMax1, kBoxMin2, kBoxMax2, n, hitrate, kCutoff, MoleculeType{}, seed);
}

/**
 * Both cells receive exactly n particles for every (n, hitrate) combination.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testParticleCount) {
  for (std::size_t n : {0ul, 1ul, 10ul, 50ul}) {
    for (double hitrate : {0.0, 0.5, 1.0}) {
      CellType v1, v2;
      fill(v1, v2, n, hitrate);
      EXPECT_EQ(v1.size(), n) << "n=" << n << " hitrate=" << hitrate;
      EXPECT_EQ(v2.size(), n) << "n=" << n << " hitrate=" << hitrate;
    }
  }
}

/**
 * Every particle lies strictly within its own cell's bounding box.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testPositionsInBounds) {
  CellType v1, v2;
  fill(v1, v2, 20, 0.5);
  for (std::size_t d = 0; d < 3; ++d) {
    for (const auto &p : v1) {
      EXPECT_GE(p.getR()[d], kBoxMin1[d]) << "cell1 below min, dim=" << d;
      EXPECT_LT(p.getR()[d], kBoxMax1[d]) << "cell1 above max, dim=" << d;
    }
    for (const auto &p : v2) {
      EXPECT_GE(p.getR()[d], kBoxMin2[d]) << "cell2 below min, dim=" << d;
      EXPECT_LT(p.getR()[d], kBoxMax2[d]) << "cell2 above max, dim=" << d;
    }
  }
}

/**
 * With hitrate=0 every particle is in the far zone: no cross-cell pair is within cutoff.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testInteractionCountHitrateZero) {
  CellType v1, v2;
  fill(v1, v2, 30, 0.0);
  for (const auto &p1 : v1) {
    for (const auto &p2 : v2) {
      EXPECT_GE(dist3D(p1, p2), kCutoff);
    }
  }
}

/**
 * At least k cross-cell pairs are within cutoff after filling.
 * Particles are shuffled so we count all in-range pairs rather than checking by index.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testPairedParticlesInRange) {
  for (std::size_t n : {1ul, 10ul, 21ul}) {
    for (double hitrate : {0.5, 1.0}) {
      CellType v1, v2;
      fill(v1, v2, n, hitrate);
      const auto k = static_cast<std::size_t>(std::round(static_cast<double>(n) * hitrate));
      std::size_t inRangeCount = 0;
      for (const auto &p1 : v1) {
        for (const auto &p2 : v2) {
          if (dist3D(p1, p2) < kCutoff) {
            ++inRangeCount;
          }
        }
      }
      EXPECT_GE(inRangeCount, k) << "n=" << n << " hitrate=" << hitrate;
    }
  }
}

/**
 * Far particles are out of range of every particle in the other cell.
 * Far particles are identified geometrically: cell1 far particles have x < boundary - cutoff,
 * cell2 far particles have x > boundary + cutoff.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testFarParticlesOutOfRange) {
  constexpr std::size_t n = 20;
  constexpr double hitrate = 0.5;
  CellType v1, v2;
  fill(v1, v2, n, hitrate);

  for (const auto &p1 : v1) {
    if (p1.getR()[0] < kBoundary - kCutoff) {
      for (const auto &p2 : v2) {
        EXPECT_GE(dist3D(p1, p2), kCutoff)
            << "far cell1 particle at x=" << p1.getR()[0] << " within cutoff of a cell2 particle";
      }
    }
  }
  for (const auto &p2 : v2) {
    if (p2.getR()[0] > kBoundary + kCutoff) {
      for (const auto &p1 : v1) {
        EXPECT_GE(dist3D(p1, p2), kCutoff)
            << "far cell2 particle at x=" << p2.getR()[0] << " within cutoff of a cell1 particle";
      }
    }
  }
}

/**
 * No two particles across both cells share an ID.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testUniqueIDs) {
  CellType v1, v2;
  fill(v1, v2, 25, 0.6);
  std::set<std::size_t> ids;
  for (const auto &p : v1) {
    EXPECT_EQ(ids.count(p.getID()), 0u) << "duplicate ID " << p.getID() << " in cell1";
    ids.insert(p.getID());
  }
  for (const auto &p : v2) {
    EXPECT_EQ(ids.count(p.getID()), 0u) << "duplicate ID " << p.getID() << " in cell2";
    ids.insert(p.getID());
  }
}

/**
 * The same seed produces identical particle positions across two independent calls.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testDeterminism) {
  CellType a1, a2, b1, b2;
  fill(a1, a2, 15, 0.4, /*seed=*/123);
  fill(b1, b2, 15, 0.4, /*seed=*/123);
  for (std::size_t i = 0; i < a1.size(); ++i) {
    EXPECT_EQ(a1[i].getR(), b1[i].getR()) << "cell1 particle " << i << " differs";
  }
  for (std::size_t i = 0; i < a2.size(); ++i) {
    EXPECT_EQ(a2[i].getR(), b2[i].getR()) << "cell2 particle " << i << " differs";
  }
}

/**
 * Non-adjacent cells (gap between boxMax1[0] and boxMin2[0]) must throw.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testPreconditionNonAdjacentCells) {
  CellType v1, v2;
  autopasTools::PseudoContainer c1(v1), c2(v2);
  const std::array<double, 3> gappedMin2 = {kBoundary + 1.0, 0., 0.};
  auto call = [&]() {
    autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
        c1, c2, kBoxMin1, kBoxMax1, gappedMin2, kBoxMax2, 10ul, 0.5, kCutoff, MoleculeType{});
  };
  EXPECT_THROW(call(), autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Hitrate outside [0, 1] must throw.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testPreconditionHitrateOutOfRange) {
  CellType v1, v2;
  autopasTools::PseudoContainer c1(v1), c2(v2);
  auto call = [&](double hitrate) {
    autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
        c1, c2, kBoxMin1, kBoxMax1, kBoxMin2, kBoxMax2, 10ul, hitrate, kCutoff, MoleculeType{});
  };
  EXPECT_THROW(call(1.5), autopas::utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(call(-0.1), autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * A cell whose x-extent is <= 1.5 * cutoff (insufficient for the far zone) must throw.
 */
TEST_F(TwoCellsInteractionHitrateGeneratorTest, testPreconditionCellTooNarrow) {
  CellType v1, v2;
  autopasTools::PseudoContainer c1(v1), c2(v2);
  // x-extent = kCutoff exactly, which is <= 1.5 * kCutoff
  const std::array<double, 3> tooNarrowMin1 = {kBoundary - kCutoff, 0., 0.};
  auto call = [&]() {
    autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
        c1, c2, tooNarrowMin1, kBoxMax1, kBoxMin2, kBoxMax2, 10ul, 0.5, kCutoff, MoleculeType{});
  };
  EXPECT_THROW(call(), autopas::utils::ExceptionHandler::AutoPasException);
}
