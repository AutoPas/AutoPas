/**
 * @file FmmOperators.h
 * @author Joachim Marin
 * @date 13.11.2019
 */

#pragma once

#include <complex>
#include <vector>
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/fastMultipoleMethod/FmmOperatorInterface.h"
#include "autopas/fastMultipoleMethod/FmmTree.h"
#include "autopas/fastMultipoleMethod/FmmTreeNode.h"
#include "autopas/utils/FmmMath.h"

namespace autopas::fmm {

template <class ParticleCell>
class PotentialOperators : public FmmOperatorInterface<ParticleCell> {
  using Complex = std::complex<double>;

 private:
  long orderOfExpansion;
  // i^(|k-m|-|k|-|m|)
  std::vector<std::vector<Complex>> powerM2L;
  autopas::utils::FmmMath<double, long> fmmMath;

 public:
  explicit PotentialOperators(long orderOfExpansion) : orderOfExpansion(orderOfExpansion) {
    this->fmmMath = autopas::utils::FmmMath<double, long>();

    // The power i^(|k-m|-|k|-|m|) is needed frequently in the M2L operator.
    this->powerM2L =
        std::vector<std::vector<Complex>>(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion * 2 + 1));
    for (long k = -orderOfExpansion; k <= orderOfExpansion; ++k) {
      for (long m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
        using namespace std::complex_literals;
        this->powerM2L[k + orderOfExpansion][m + orderOfExpansion] =
            autopas::utils::FmmMath<double, long>::powI(std::abs(k - m) - std::abs(k) - std::abs(m));
      }
    }
  }

  void P2M(FmmTreeNode<ParticleCell> &leaf) override {
    // Loop order changed, so the spherical harmonics cache is built less frequently.

    for (auto iter = leaf.getTree().getContainer()->getRegionIterator(leaf.getBoxMin(), leaf.getBoxMax());
         iter.isValid(); ++iter) {
      auto spherical = autopas::utils::FmmMath<double, long>::toSpherical(
          autopas::utils::ArrayMath::sub(iter->getR(), leaf.getBoxCenter()));
      for (long m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
        fmmMath.sphericalHarmonicsBuildCache(-m, orderOfExpansion, spherical[1], spherical[2]);
        for (long n = std::abs(m); n <= orderOfExpansion; n++) {
          auto add = iter->charge * std::pow(spherical[0], n) *
                     fmmMath.sphericalHarmonicsCached(-m, n, spherical[1], spherical[2]);
          leaf.setM(m, n, leaf.getM(m, n) + add);
        }
      }
    }
  }

  void M2M(FmmTreeNode<ParticleCell> &parent) override {
    if (parent.isLeaf()) {
      return;
    }
    for (std::size_t c = 0; c < 8; ++c) {
      auto child = parent.getChild(c);
      /*if (child->isZeroM()) {
        continue;
      }*/
      auto spherical = autopas::utils::FmmMath<double, long>::toSpherical(
          autopas::utils::ArrayMath::sub(child.getBoxCenter(), parent.getBoxCenter()));
      double r = spherical[0];
      double theta = spherical[1];
      double phi = spherical[2];
      for (long j = 0; j <= orderOfExpansion; ++j) {
        for (long k = -j; k <= j; ++k) {
          Complex mmn = 0;
          for (long n = 0; n <= j; ++n) {
            for (long m = -n; m <= n; ++m) {
              auto childM = child.getM(k - m, j - n);

              using namespace std::complex_literals;
              auto complex = std::pow(1i, std::abs(k) - std::abs(m) - std::abs(k - m));
              Complex product = childM * complex * std::pow(r, n);

              if (product != 0.0) {
                auto harmonics = fmmMath.sphericalHarmonics(-m, n, theta, phi);
                product *= harmonics * autopas::utils::FmmMath<double, long>::getA(m, n) *
                           autopas::utils::FmmMath<double, long>::getA(k - m, j - n) /
                           autopas::utils::FmmMath<double, long>::getA(k, j);
              }

              mmn += product;
              assert(!__isnan(mmn.real()) && !__isnan(mmn.imag()));
            }
          }
          // Add the matrix defined in 5.22
          parent.setM(k, j, parent.getM(k, j) + mmn);
        }
      }
    }
  }

  // Most performance critical function. Over 90% of fmm is spent here.
  void M2L(FmmTreeNode<ParticleCell> &node) override {
    for (auto inter : node.getInteractionList()) {
      auto spherical = autopas::utils::FmmMath<double, long>::toSpherical(
          autopas::utils::ArrayMath::sub(inter->getBoxCenter(), node.getBoxCenter()));
      double rho = spherical[0];
      double theta = spherical[1];
      double phi = spherical[2];

      // Loop order changed, so the spherical harmonics cache is built less frequently.
      for (long k = -orderOfExpansion; k <= orderOfExpansion; ++k) {
        long absK = std::abs(k);
        for (long m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
          // rhoPower1 * rhoPower2 = std::pow(rho, n+j+1)
          // Avoids std::pow in the innermost loop.
          double rhoPower1 = std::pow(rho, absK + 1);
          long absM = std::abs(m);
          fmmMath.sphericalHarmonicsBuildCache(m - k, 2 * orderOfExpansion, theta, phi);
          for (long j = absK; j <= orderOfExpansion; ++j) {
            // sign = (-1)^n
            // sign has to be reset before after every n-loop, so it cannot be put outside the j-loop.
            auto sign = absM % 2 == 1 ? 1 : -1;

            double rhoPower2 = std::pow(rho, absM) * rhoPower1;

            // complex = i^(|k-m|-|k|-|m|)
            auto complex = this->powerM2L[k + orderOfExpansion][m + orderOfExpansion];

            auto aKJ = autopas::utils::FmmMath<double, long>::getA(k, j);

            Complex sum = 0;
            for (long n = absM; n <= orderOfExpansion; ++n) {
              sign = -sign;
              auto product = inter->getM(m, n);
              if (product != 0.0) {
                product *= complex * autopas::utils::FmmMath<double, long>::getA(m, n) * aKJ *
                           fmmMath.sphericalHarmonicsCached(m - k, j + n, theta, phi);
                product /= sign * autopas::utils::FmmMath<double, long>::getA(m - k, j + n) * rhoPower2;
              }

              sum += product;
              // assert(!__isnan(sum.real()) && !__isnan(sum.imag()));

              rhoPower2 *= rho;
            }

            node.setL(k, j, node.getL(k, j) + sum);
            rhoPower1 *= rho;
          }
        }
      }
    }
  }

  void L2L(FmmTreeNode<ParticleCell> &node) override {
    if (node.getDepth() > 0) {
      auto parent = node.getOctreeParent();
      auto cartesian = autopas::utils::ArrayMath::sub(parent.getBoxCenter(), node.getBoxCenter());
      auto spherical = autopas::utils::FmmMath<double, long>::toSpherical(cartesian);
      double rho = spherical[0];
      double theta = spherical[1];
      double phi = spherical[2];
      for (long j = 0; j <= orderOfExpansion; ++j) {
        for (long k = -j; k <= j; ++k) {
          Complex lmn = 0;
          for (long n = j; n <= orderOfExpansion; ++n) {
            for (long m = -n; m <= n; ++m) {
              using namespace std::complex_literals;

              auto parentL = parent.getL(m, n);

              auto complex = std::pow(1i, std::abs(m) - std::abs(m - k) - std::abs(k));

              Complex product = parentL * complex * autopas::utils::FmmMath<double, long>::getA(m - k, n - j);
              if (product != 0.0) {
                auto harmonics = fmmMath.sphericalHarmonics(m - k, n - j, theta, phi);
                product *= autopas::utils::FmmMath<double, long>::getA(k, j) * harmonics * std::pow(rho, n - j) /
                           (std::pow(-1, n + j) * autopas::utils::FmmMath<double, long>::getA(m, n));
              }

              lmn += product;
              assert(!__isnan(lmn.real()) && !__isnan(lmn.imag()));
            }
          }
          // Add the matrix defined in 5.30
          node.setL(k, j, node.getL(k, j) + lmn);
        }
      }
    }
  }

  void L2P(FmmTreeNode<ParticleCell> &leaf) override {
    for (auto iter = leaf.getTree().getContainer()->getRegionIterator(leaf.getBoxMin(), leaf.getBoxMax());
         iter.isValid(); ++iter) {
      auto sphericalPos = autopas::utils::FmmMath<double, long>::toSpherical(
          autopas::utils::ArrayMath::sub(iter->getR(), leaf.getBoxCenter()));

      double rho = sphericalPos[0];
      double theta = sphericalPos[1];
      double phi = sphericalPos[2];

      // evaluate
      Complex potential = 0;
      for (long n = 0; n <= orderOfExpansion; ++n) {
        for (long m = -n; m <= n; m++) {
          Complex product = std::pow(rho, n);

          if (product != 0.0) {
            product *= fmmMath.sphericalHarmonics(m, n, theta, phi);
            if (product != 0.0) {
              product *= leaf.getL(m, n);
            }
          }
          potential += product;
          assert(!__isnan(potential.real()) && !__isnan(potential.imag()));
        }
      }
      if (std::abs(potential.imag()) > 0.00001) {
        std::cerr << "Potential has non-zero imaginary part: " << potential << std::endl;
      }
      iter->longRange = potential.real();
      iter->resultFMM = potential.real();
    }
  }
};
}  // namespace autopas::fmm
