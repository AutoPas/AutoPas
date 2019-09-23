/**
 * @file Operators.cpp
 * @date 15.09.19
 * @author Joachim Marin
 */

#include "Operators.h"
#include "Math3D.h"

void printVec(const std::array<double, 3> &vec) {
  std::cout << "(" << vec[0] << "," << vec[1] << "," << vec[2] << ")";
}

double getA(int m, int n) {

  if (n - m < 0 || n + m < 0) {
    std::cerr << "getA(" << m << "," << n << ") is not defined for n - m < 0 or n + m < 0" << std::endl;
  }
  return std::pow(-1, n) / std::sqrt(factorial(n - m) * factorial(n + m));
}

void Operators::P2M(OctreeNode *leaf) {

  leaf->fmmM = ComplexMatrix(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion + 1, 0));

  for (int n = 0; n <= orderOfExpansion; ++n) {
    for (int m = -n; m <= n; ++m) {
      Complex mmn = 0;
      for (auto iter = leaf->getContainer()->begin(); iter.isValid(); ++iter) {
        auto sphericalPos = toSpherical(center(iter->getR(), leaf->getCenter()));


        mmn += iter->charge * std::pow(sphericalPos[0], n) *
               sphericalHarmonics(-m, n, sphericalPos[1], sphericalPos[2]);
      }

      // Save the matrix defined in 5.16
      leaf->setM(m, n, mmn);
    }
  }
}


void Operators::M2M(OctreeNode *parent) {

  parent->fmmM = ComplexMatrix(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion + 1, 0));

  for (int j = 0; j <= orderOfExpansion; ++j) {
    for (int k = -j; k <= j; ++k) {
      parent->setM(k, j, 0);
    }
  }

  for (int c = 0; c < 8; ++c) {
    auto child = parent->getChild(c);
    if(child->getIsZeroM())
    {
      /*std::cout << "skip child with depth = " << child->getDepth() << " and center = ";
      printVec(child->getCenter());
      std::cout << std::endl;*/
      continue;
    }
    auto cartesian = center(child->getCenter(), parent->getCenter());
    auto spherical = toSpherical(cartesian);
    double r = spherical[0];
    double alpha = spherical[1];
    double beta = spherical[2];
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
              auto harmonics = sphericalHarmonics(-m, n, alpha, beta);
              product *= harmonics * getA(m, n) * getA(k - m, j - n) / getA(k, j);
            }

            mmn += product;
          }
        }
        // Add the matrix defined in 5.22
        parent->setM(k, j, parent->getM(k, j) + mmn);
      }
    }
  }
}

void Operators::M2L(OctreeNode *node) {

  node->fmmL = ComplexMatrix(orderOfExpansion * 2 + 1, std::vector<Complex>(orderOfExpansion + 1, 0));

  // No far field interactions.
  if (node->getDepth() < 2) {
    return;
  }

  for (auto inter : *node->getInteractionList()) {
    if(inter->getIsZeroM()) {
      /*std::cout << "skip interactor with depth = " << inter->getDepth() << " and center = ";
      printVec(inter->getCenter());
      std::cout << std::endl;*/
      continue;
    }
    auto cartesian = center(inter->getCenter(), node->getCenter());
    auto spherical = toSpherical(cartesian);
    double r = spherical[0];
    double alpha = spherical[1];
    double beta = spherical[2];
    for (int j = 0; j <= orderOfExpansion; ++j) {
      for (int k = -j; k <= j; ++k) {
        Complex lmn = 0;
        for (int n = 0; n <= orderOfExpansion; ++n) {
          for (int m = -n; m <= n; ++m) {
            using namespace std::complex_literals;
            auto interM = inter->getM(m, n);

            Complex product = interM;
            if (product != 0.0) {
              auto complex = std::pow(1i, std::abs(k - m) - std::abs(k) - std::abs(m));
              product *= complex * getA(m, n) * getA(k, j);
              if (product != 0.0) {
                auto harmonics = sphericalHarmonics(m - k, j + n, alpha, beta);
                product *= harmonics;
                product /= (std::pow(-1, n) * getA(m - k, j + n) * std::pow(r, j + n + 1));
              }
            }
            lmn += product;
          }
        }

        // Add the matrix defined in 5.26
        node->setL(k, j, node->getL(k, j) + lmn);
      }
    }
  }

}

void Operators::L2L(OctreeNode *node) {
  auto parent = node->getParent();
  if (parent != nullptr) {
    if(parent->getIsZeroL()){
      /*std::cout << "skip parent with depth = " << parent->getDepth() << " and center = ";
      printVec(parent->getCenter());
      std::cout << std::endl;*/
     return;
    }
    auto cartesian = center(parent->getCenter(), node->getCenter());
    auto spherical = toSpherical(cartesian);
    double r = spherical[0];
    double alpha = spherical[1];
    double beta = spherical[2];
    for (int j = 0; j <= orderOfExpansion; ++j) {
      for (int k = -j; k <= j; ++k) {
        Complex lmn = 0;
        for (int n = j; n <= orderOfExpansion; ++n) {
          for (int m = -n; m <= n; ++m) {
            using namespace std::complex_literals;

            auto parentL = parent->getL(m, n);

            auto complex = std::pow(1i, std::abs(m) - std::abs(m - k) - std::abs(k));

            Complex product = parentL * complex;
            if (product != 0.0) {

              auto harmonics = sphericalHarmonics(m - k, n - j, alpha, beta);
              product *= getA(m - k, n - j) * getA(k, j) * harmonics * std::pow(r, n - j) /
                         (std::pow(-1, n + j) * getA(m, n));
            }

            lmn += product;

          }
        }
        // Add the matrix defined in 5.30
        node->setL(k, j, node->getL(k, j) + lmn);

      }
    }
  }
}

void Operators::L2P(OctreeNode *leaf) {
  for (auto iter = leaf->getContainer()->begin(); iter.isValid(); ++iter) {
    auto sphericalPos = toSpherical(center(iter->getR(), leaf->getCenter()));

    double r = sphericalPos[0];
    double alpha = sphericalPos[1];
    double beta = sphericalPos[2];

    /*printVec(iter->getR());
    std::cout << ", ";
    printVec(leaf->getCenter());
    std::cout << ", ";
    printVec(center(iter->getR(), leaf->getCenter()));
    std::cout << ", ";
    printVec(sphericalPos);
    std::cout << std::endl;*/

    // evaluate
    Complex potential = 0;
    for (int n = 0; n <= orderOfExpansion; ++n) {
      for (int m = -n; m <= n; m++) {
        Complex product = std::pow(r, n);

        if (product != 0.0) {
          product *= sphericalHarmonics(m, n, alpha, beta);
          if (product != 0.0) {
            product *= leaf->getL(m, n);
          }
        }
        //std::cout << potential << "+" << product << "=" << potential + product << std::endl;
        potential += product;
      }
    }
    if (std::abs(potential.imag()) > 0.00001) {
      std::cerr << "Potential has non-zero imaginary part: " << potential << std::endl;
    }

    iter->resultFMM = potential.real();

    std::cout << "result = " << iter->resultFMM << std::endl;
  }
}


