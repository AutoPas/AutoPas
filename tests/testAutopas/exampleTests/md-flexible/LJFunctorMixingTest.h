#pragma once
#include <gtest/gtest.h>
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "AutoPasTestBase.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class LJFunctorMixingTest : public AutoPasTestBase {
 public:
  LJFunctorMixingTest()
      : AutoPasTestBase(),
        epsilon1{1.0},
        sigma1{1.0},
        cutoff{1.},
        boxmin{{0., 0., 0.}},
        boxmax{{5., 5., 5.}},
        PCL{ParticleClassLibrary(epsilonMapForTwoParticles, sigmaMapForTwoParticles, massMapForTwoParticles)},
        functor{autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>,
                                   autopas::FunctorN3Modes::Both, true>(cutoff, PCL, 0.0)} {}

  static double L2Norm(std::array<double, 3> array) {
    double square_sum = 0;
    for (double i : array) {
      square_sum += (i * i);
    }
    return sqrt(square_sum);
  }

  /**Berechnet lennordForce zwischen zwei particles i und j mit positionen xi und xj, gibt modifiziertes xi zurück
   * */
  std::array<double, 3> lennardForceCalculation(std::array<double, 3> x1, std::array<double, 3> x2);

  void testAoSNoGlobals(bool newton3);

 protected:
  double epsilon1 = 1.0;
  double sigma1 = 1.0;
  double epsilon2 = 2.0;
  double sigma2 = 2.0;
  double shift = 0.1;
  double absDelta = 1e-7;
  std::array<double, 3> expectedForce = {835415983.76769447326660156, 1670831967.53538894653320312,
                                         2506247951.30308341979980469};
  std::map<unsigned long, double> epsilonMapForTwoParticles = {{0, epsilon1}, {1, epsilon2}};
  std::map<unsigned long, double> sigmaMapForTwoParticles = {{0, sigma1}, {1, sigma2}};
  std::map<unsigned long, double> massMapForTwoParticles = {{0, 1.0}, {1, 1.0}};
  double cutoff;
  std::array<double, 3> boxmin;
  std::array<double, 3> boxmax;
  ParticleClassLibrary PCL;
  autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both, true>
      functor;
};