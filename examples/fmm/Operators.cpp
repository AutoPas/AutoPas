/**
 * @file Operators.cpp
 * @date 15.09.19
 * @author Joachim Marin
 */

#include "Operators.h"

Operators::Operators(int orderOfExpansion) {
  this->orderOfExpansion = orderOfExpansion;
  this->mathOpt = Math3D::MathOpt();

  // The power i^(|k-m|-|k|-|m|) is needed frequently in the M2L operator.
  this->powerM2L =
      std::vector<std::vector<Complex>>(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion * 2 + 1));
  for (int k = -orderOfExpansion; k <= orderOfExpansion; ++k) {
    for (int m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
      using namespace std::complex_literals;
      this->powerM2L[k + orderOfExpansion][m + orderOfExpansion] =
          Math3D::MathOpt::powI(std::abs(k - m) - std::abs(k) - std::abs(m));
    }
  }
}

void Operators::P2M(AdaptiveOctreeNode &leaf) {
  // Loop order changed, so the spherical harmonics cache is built less frequently.
  for (auto iter = leaf.getTree()->getDomain()->getRegionIterator(leaf.getNodeMinCorner(), leaf.getNodeMaxCorner());
       iter.isValid(); ++iter) {
    auto spherical = Math3D::toSpherical(autopas::ArrayMath::sub(iter->getR(), leaf.getNodeCenter()));
    for (int m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
      mathOpt.sphericalHarmonicsBuildCache(-m, orderOfExpansion, spherical[1], spherical[2]);
      for (int n = std::abs(m); n <= orderOfExpansion; n++) {
        auto add = iter->charge * std::pow(spherical[0], n) *
            mathOpt.sphericalHarmonicsCached(-m, n, spherical[1], spherical[2]);
        leaf.setM(m, n, leaf.getM(m, n) + add);
      }
    }
  }
}

void Operators::M2M(AdaptiveOctreeNode &parent) {
  if (parent.isLeaf()) {
    return;
  }
  for (int c = 0; c < 8; ++c) {
    auto child = parent.getChild(c);
    if (child->isZeroM()) {
      continue;
    }
    auto spherical = Math3D::toSpherical(autopas::ArrayMath::sub(child->getNodeCenter(), parent.getNodeCenter()));
    double r = spherical[0];
    double theta = spherical[1];
    double phi = spherical[2];
    for (int j = 0; j <= orderOfExpansion; ++j) {
      for (int k = -j; k <= j; ++k) {
        Complex mmn = 0;
        for (int n = 0; n <= j; ++n) {
          for (int m = -n; m <= n; ++m) {
            auto childM = child->getM(k - m, j - n);

            using namespace std::complex_literals;
            auto complex = std::pow(1i, std::abs(k) - std::abs(m) - std::abs(k - m));
            Complex product = childM * complex * std::pow(r, n);

            if (product != 0.0) {
              auto harmonics = mathOpt.sphericalHarmonics(-m, n, theta, phi);
              product *= harmonics * Math3D::MathOpt::getA(m, n) * Math3D::MathOpt::getA(k - m, j - n) / Math3D::MathOpt::getA(k, j);
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
void Operators::M2L(AdaptiveOctreeNode &node) {
  
  for (auto inter : node.getInteractionList()) {
    if (inter->isZeroM()) {
      continue;
    }
    auto spherical = Math3D::toSpherical(autopas::ArrayMath::sub(inter->getNodeCenter(), node.getNodeCenter()));
    double rho = spherical[0];
    double theta = spherical[1];
    double phi = spherical[2];

    // Loop order changed, so the spherical harmonics cache is built less frequently.
    for (int k = -orderOfExpansion; k <= orderOfExpansion; ++k) {
      int absK = std::abs(k);
      for (int m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
        // rhoPower1 * rhoPower2 = std::pow(rho, n+j+1)
        // Avoids std::pow in the innermost loop.
        double rhoPower1 = std::pow(rho, absK + 1);
        int absM = std::abs(m);
        mathOpt.sphericalHarmonicsBuildCache(m - k, 2 * orderOfExpansion, theta, phi);
        for (int j = absK; j <= orderOfExpansion; ++j) {
          // sign = (-1)^n
          // sign has to be reset before after every n-loop, so it cannot be put outside the j-loop.
          auto sign = absM % 2 == 1 ? 1 : -1;

          double rhoPower2 = std::pow(rho, absM) * rhoPower1;

          // complex = i^(|k-m|-|k|-|m|)
          auto complex = this->powerM2L[k + orderOfExpansion][m + orderOfExpansion];

          auto aKJ = Math3D::MathOpt::getA(k, j);

          Complex sum = 0;
          for (int n = absM; n <= orderOfExpansion; ++n) {
            sign = -sign;
            auto product = inter->getM(m, n);
            if (product != 0.0) {
              product *=
                  complex * Math3D::MathOpt::getA(m, n) * aKJ * mathOpt.sphericalHarmonicsCached(m - k, j + n, theta, phi);
              product /= sign * Math3D::MathOpt::getA(m - k, j + n) * rhoPower2;
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

void Operators::L2L(AdaptiveOctreeNode &node) {
  auto parent = node.getParent();
  if (parent != nullptr) {
    if (parent->isZeroL()) {
      return;
    }
    auto cartesian = autopas::ArrayMath::sub(parent->getNodeCenter(), node.getNodeCenter());
    auto spherical = Math3D::toSpherical(cartesian);
    double rho = spherical[0];
    double theta = spherical[1];
    double phi = spherical[2];
    for (int j = 0; j <= orderOfExpansion; ++j) {
      for (int k = -j; k <= j; ++k) {
        Complex lmn = 0;
        for (int n = j; n <= orderOfExpansion; ++n) {
          for (int m = -n; m <= n; ++m) {
            using namespace std::complex_literals;

            auto parentL = parent->getL(m, n);

            auto complex = std::pow(1i, std::abs(m) - std::abs(m - k) - std::abs(k));

            Complex product = parentL * complex * Math3D::MathOpt::getA(m - k, n - j);
            if (product != 0.0) {
              auto harmonics = mathOpt.sphericalHarmonics(m - k, n - j, theta, phi);
              product *=
                  Math3D::MathOpt::getA(k, j) * harmonics * std::pow(rho, n - j) / (std::pow(-1, n + j) * Math3D::MathOpt::getA(m, n));
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

void Operators::L2P(AdaptiveOctreeNode &leaf) {
  for (auto iter = leaf.getTree()->getDomain()->getRegionIterator(leaf.getNodeMinCorner(), leaf.getNodeMaxCorner());
       iter.isValid(); ++iter) {
    auto sphericalPos = Math3D::toSpherical(autopas::ArrayMath::sub(iter->getR(), leaf.getNodeCenter()));

    double rho = sphericalPos[0];
    double theta = sphericalPos[1];
    double phi = sphericalPos[2];

    // evaluate
    Complex potential = 0;
    for (int n = 0; n <= orderOfExpansion; ++n) {
      for (int m = -n; m <= n; m++) {
        Complex product = std::pow(rho, n);

        if (product != 0.0) {
          product *= mathOpt.sphericalHarmonics(m, n, theta, phi);
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

// compare old:

/**
 * @file Operators.cpp
 * @date 15.09.19
 * @author Joachim Marin
 */

void Operators::P2M_Old(OctreeNode &leaf) {
  leaf.fmmM = ComplexMatrix(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion + 1, 0));

  // Loop order changed, so the spherical harmonics cache is built less frequently.
  for (auto iter = leaf.getContainer()->getRegionIterator(leaf.getLowCorner(), leaf.getHighCorner()); iter.isValid();
       ++iter) {
    auto spherical = Math3D::toSpherical(autopas::ArrayMath::sub(iter->getR(), leaf.getCenter()));
    for (int m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
      mathOpt.sphericalHarmonicsBuildCache(-m, orderOfExpansion, spherical[1], spherical[2]);
      for (int n = std::abs(m); n <= orderOfExpansion; n++) {
        auto add = iter->charge * std::pow(spherical[0], n) *
            mathOpt.sphericalHarmonicsCached(-m, n, spherical[1], spherical[2]);
        leaf.setM(m, n, leaf.getM(m, n) + add);
      }
    }
  }
}

void Operators::M2M_Old(OctreeNode &parent) {
  parent.fmmM = ComplexMatrix(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion + 1, 0));

  for (int c = 0; c < 8; ++c) {
    auto child = parent.getChild(c);
    if (child->getIsZeroM()) {
      continue;
    }
    auto spherical = Math3D::toSpherical(autopas::ArrayMath::sub(child->getCenter(), parent.getCenter()));
    double r = spherical[0];
    double theta = spherical[1];
    double phi = spherical[2];
    for (int j = 0; j <= orderOfExpansion; ++j) {
      for (int k = -j; k <= j; ++k) {
        Complex mmn = 0;
        for (int n = 0; n <= j; ++n) {
          for (int m = -n; m <= n; ++m) {
            auto childM = child->getM(k - m, j - n);

            using namespace std::complex_literals;
            auto complex = std::pow(1i, std::abs(k) - std::abs(m) - std::abs(k - m));
            Complex product = childM * complex * std::pow(r, n);

            if (product != 0.0) {
              auto harmonics = mathOpt.sphericalHarmonics(-m, n, theta, phi);
              product *= harmonics * Math3D::MathOpt::getA(m, n) * Math3D::MathOpt::getA(k - m, j - n) / Math3D::MathOpt::getA(k, j);
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
void Operators::M2L_Old(OctreeNode &node) {
  node.fmmL = ComplexMatrix(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion + 1, 0));

  // No far field interactions.
  if (node.getDepth() < 2) {
    return;
  }

  for (auto inter : *node.getInteractionList()) {
    if (inter->getIsZeroM()) {
      continue;
    }
    auto spherical = Math3D::toSpherical(autopas::ArrayMath::sub(inter->getCenter(), node.getCenter()));
    double rho = spherical[0];
    double theta = spherical[1];
    double phi = spherical[2];

    // Loop order changed, so the spherical harmonics cache is built less frequently.
    for (int k = -orderOfExpansion; k <= orderOfExpansion; ++k) {
      int absK = std::abs(k);
      for (int m = -orderOfExpansion; m <= orderOfExpansion; ++m) {
        // rhoPower1 * rhoPower2 = std::pow(rho, n+j+1)
        // Avoids std::pow in the innermost loop.
        double rhoPower1 = std::pow(rho, absK + 1);
        int absM = std::abs(m);
        mathOpt.sphericalHarmonicsBuildCache(m - k, 2 * orderOfExpansion, theta, phi);
        for (int j = absK; j <= orderOfExpansion; ++j) {
          // sign = (-1)^n
          // sign has to be reset before after every n-loop, so it cannot be put outside the j-loop.
          auto sign = absM % 2 == 1 ? 1 : -1;

          double rhoPower2 = std::pow(rho, absM) * rhoPower1;

          // complex = i^(|k-m|-|k|-|m|)
          auto complex = this->powerM2L[k + orderOfExpansion][m + orderOfExpansion];

          auto aKJ = Math3D::MathOpt::getA(k, j);

          Complex sum = 0;
          for (int n = absM; n <= orderOfExpansion; ++n) {
            sign = -sign;
            auto product = inter->getM(m, n);
            if (product != 0.0) {
              product *=
                  complex * Math3D::MathOpt::getA(m, n) * aKJ * mathOpt.sphericalHarmonicsCached(m - k, j + n, theta, phi);
              product /= sign * Math3D::MathOpt::getA(m - k, j + n) * rhoPower2;
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

void Operators::L2L_Old(OctreeNode &node) {
  auto parent = node.getParent();
  if (parent != nullptr) {
    if (parent->getIsZeroL()) {
      return;
    }
    auto cartesian = autopas::ArrayMath::sub(parent->getCenter(), node.getCenter());
    auto spherical = Math3D::toSpherical(cartesian);
    double rho = spherical[0];
    double theta = spherical[1];
    double phi = spherical[2];
    for (int j = 0; j <= orderOfExpansion; ++j) {
      for (int k = -j; k <= j; ++k) {
        Complex lmn = 0;
        for (int n = j; n <= orderOfExpansion; ++n) {
          for (int m = -n; m <= n; ++m) {
            using namespace std::complex_literals;

            auto parentL = parent->getL(m, n);

            auto complex = std::pow(1i, std::abs(m) - std::abs(m - k) - std::abs(k));

            Complex product = parentL * complex * Math3D::MathOpt::getA(m - k, n - j);
            if (product != 0.0) {
              auto harmonics = mathOpt.sphericalHarmonics(m - k, n - j, theta, phi);
              product *=
                  Math3D::MathOpt::getA(k, j) * harmonics * std::pow(rho, n - j) / (std::pow(-1, n + j) * Math3D::MathOpt::getA(m, n));
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

void Operators::L2P_Old(OctreeNode &leaf) {
  for (auto iter = leaf.getContainer()->getRegionIterator(leaf.getLowCorner(), leaf.getHighCorner()); iter.isValid();
       ++iter) {
    auto sphericalPos = Math3D::toSpherical(autopas::ArrayMath::sub(iter->getR(), leaf.getCenter()));

    double rho = sphericalPos[0];
    double theta = sphericalPos[1];
    double phi = sphericalPos[2];

    // evaluate
    Complex potential = 0;
    for (int n = 0; n <= orderOfExpansion; ++n) {
      for (int m = -n; m <= n; m++) {
        Complex product = std::pow(rho, n);

        if (product != 0.0) {
          product *= mathOpt.sphericalHarmonics(m, n, theta, phi);
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

    iter->longRange_Old = potential.real();
    iter->resultFMM_Old = potential.real();
  }
}