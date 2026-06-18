/**
 * @file TwoCellsInteractionRatioGeneratorTest.cpp
 * @author H. Meyran
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <set>
#include <vector>

#include "autopas/utils/ExceptionHandler.h"
#include "autopasTools/PseudoContainer.h"
#include "autopasTools/generators/TwoCellsInteractionHitrateGenerator.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

using Molecule = mdLib::MoleculeLJ;
using Cell = std::vector<Molecule>;

namespace {

constexpr double kCutoff = 3.0;
constexpr double kBoundary = 10.0;
constexpr double kExtent = 8.0;  // > 1.5 * kCutoff = 4.5

const std::array<double, 3> kBoxMin1 = {kBoundary - kExtent, 0., 0.};
const std::array<double, 3> kBoxMax1 = {kBoundary, 8., 8.};
const std::array<double, 3> kBoxMin2 = {kBoundary, 0., 0.};
const std::array<double, 3> kBoxMax2 = {kBoundary + kExtent, 8., 8.};

double dist3D(const Molecule &a, const Molecule &b) {
  double d = 0.;
  for (int dim = 0; dim < 3; ++dim) {
    const double delta = a.getR()[dim] - b.getR()[dim];
    d += delta * delta;
  }
  return std::sqrt(d);
}

void fill(Cell &v1, Cell &v2, std::size_t n, double ratio, unsigned int seed = 42) {
  autopasTools::PseudoContainer c1(v1), c2(v2);
  autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      c1, c2, kBoxMin1, kBoxMax1, kBoxMin2, kBoxMax2, n, ratio, kCutoff, Molecule{}, seed);
}

}  // namespace

// Both cells receive exactly n particles for every (n, ratio) combination.
TEST(TwoCellsInteractionRatioGeneratorTest, ParticleCount) {
  for (std::size_t n : {0ul, 1ul, 10ul, 50ul}) {
    for (double ratio : {0.0, 0.5, 1.0}) {
      Cell v1, v2;
      fill(v1, v2, n, ratio);
      EXPECT_EQ(v1.size(), n) << "n=" << n << " ratio=" << ratio;
      EXPECT_EQ(v2.size(), n) << "n=" << n << " ratio=" << ratio;
    }
  }
}

// Every particle lies strictly within its own cell's bounding box.
TEST(TwoCellsInteractionRatioGeneratorTest, PositionsInBounds) {
  Cell v1, v2;
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

// With ratio=0 every particle is in the far zone: no cross-cell pair is within cutoff.
TEST(TwoCellsInteractionRatioGeneratorTest, InteractionCountRatioZero) {
  Cell v1, v2;
  fill(v1, v2, 30, 0.0);
  for (const auto &p1 : v1) {
    for (const auto &p2 : v2) {
      EXPECT_GE(dist3D(p1, p2), kCutoff);
    }
  }
}

// The k explicitly paired partners (cell1[i], cell2[i]) are within cutoff.
TEST(TwoCellsInteractionRatioGeneratorTest, PairedParticlesInRange) {
  for (std::size_t n : {1ul, 10ul, 21ul}) {
    for (double ratio : {0.5, 1.0}) {
      Cell v1, v2;
      fill(v1, v2, n, ratio);
      const auto k = static_cast<std::size_t>(std::round(static_cast<double>(n) * ratio));
      for (std::size_t i = 0; i < k; ++i) {
        EXPECT_LT(dist3D(v1[i], v2[i]), kCutoff)
            << "pair " << i << " out of range (n=" << n << " ratio=" << ratio << ")";
      }
    }
  }
}

// Far particles (indices k..n-1) are out of range of every particle in the other cell.
// This holds because x-distance alone exceeds cutoff for all far↔any-other-cell pairs.
TEST(TwoCellsInteractionRatioGeneratorTest, FarParticlesOutOfRange) {
  constexpr std::size_t n = 20;
  constexpr double ratio = 0.5;
  Cell v1, v2;
  fill(v1, v2, n, ratio);
  const auto k = static_cast<std::size_t>(std::round(static_cast<double>(n) * ratio));

  for (std::size_t i = k; i < n; ++i) {
    for (const auto &p2 : v2) {
      EXPECT_GE(dist3D(v1[i], p2), kCutoff) << "far cell1[" << i << "] within cutoff of a cell2 particle";
    }
  }
  for (std::size_t i = k; i < n; ++i) {
    for (const auto &p1 : v1) {
      EXPECT_GE(dist3D(v2[i], p1), kCutoff) << "far cell2[" << i << "] within cutoff of a cell1 particle";
    }
  }
}

// No two particles across both cells share an ID.
TEST(TwoCellsInteractionRatioGeneratorTest, UniqueIDs) {
  Cell v1, v2;
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

// The same seed produces identical particle positions across two independent calls.
TEST(TwoCellsInteractionRatioGeneratorTest, Determinism) {
  Cell a1, a2, b1, b2;
  fill(a1, a2, 15, 0.4, /*seed=*/123);
  fill(b1, b2, 15, 0.4, /*seed=*/123);
  for (std::size_t i = 0; i < a1.size(); ++i) {
    EXPECT_EQ(a1[i].getR(), b1[i].getR()) << "cell1 particle " << i << " differs";
  }
  for (std::size_t i = 0; i < a2.size(); ++i) {
    EXPECT_EQ(a2[i].getR(), b2[i].getR()) << "cell2 particle " << i << " differs";
  }
}

// Non-adjacent cells (gap between boxMax1[0] and boxMin2[0]) must throw.
TEST(TwoCellsInteractionRatioGeneratorTest, PreconditionNonAdjacentCells) {
  Cell v1, v2;
  autopasTools::PseudoContainer c1(v1), c2(v2);
  const std::array<double, 3> gappedMin2 = {kBoundary + 1.0, 0., 0.};
  auto call = [&]() {
    autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
        c1, c2, kBoxMin1, kBoxMax1, gappedMin2, kBoxMax2, 10ul, 0.5, kCutoff, Molecule{});
  };
  EXPECT_THROW(call(), autopas::utils::ExceptionHandler::AutoPasException);
}

// Ratio outside [0, 1] must throw.
TEST(TwoCellsInteractionRatioGeneratorTest, PreconditionRatioOutOfRange) {
  Cell v1, v2;
  autopasTools::PseudoContainer c1(v1), c2(v2);
  auto call = [&](double ratio) {
    autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
        c1, c2, kBoxMin1, kBoxMax1, kBoxMin2, kBoxMax2, 10ul, ratio, kCutoff, Molecule{});
  };
  EXPECT_THROW(call(1.5), autopas::utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(call(-0.1), autopas::utils::ExceptionHandler::AutoPasException);
}

// A cell whose x-extent is <= 1.5 * cutoff (insufficient for the far zone) must throw.
TEST(TwoCellsInteractionRatioGeneratorTest, PreconditionCellTooNarrow) {
  Cell v1, v2;
  autopasTools::PseudoContainer c1(v1), c2(v2);
  // x-extent = kCutoff exactly, which is <= 1.5 * kCutoff
  const std::array<double, 3> tooNarrowMin1 = {kBoundary - kCutoff, 0., 0.};
  auto call = [&]() {
    autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
        c1, c2, tooNarrowMin1, kBoxMax1, kBoxMin2, kBoxMax2, 10ul, 0.5, kCutoff, Molecule{});
  };
  EXPECT_THROW(call(), autopas::utils::ExceptionHandler::AutoPasException);
}
