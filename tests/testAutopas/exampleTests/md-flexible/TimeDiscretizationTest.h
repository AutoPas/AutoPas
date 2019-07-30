#include <gtest/gtest.h>
#include <vector>
#include "../../../../examples/md-flexible/PrintableMolecule.h"
#include "../../../../examples/md-flexible/Simulation.h"
#include "../../../../examples/md-flexible/TimeDiscretization.h"
#include "../../../../src/autopas/utils/ArrayMath.h"
#include "../../testingHelpers/GridGenerator.h"
#include "../../testingHelpers/RandomGenerator.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

class TimeDiscretizationTest : public AutoPasTestBase {
 public:
  TimeDiscretizationTest()
      : AutoPasTestBase(),
        epsilon{1.0},
        sigma{1.0},
        cutoff{1.},
        boxmin{{0., 0., 0.}},
        boxmax{{5., 5., 5.}},
        PCL{ParticleClassLibrary(epsilon, sigma, 1.0},
        functor{autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>,
                                   autopas::FunctorN3Modes::Both, true>(cutoff, PCL, 0.0)} {}

  void globalForceTest(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &auto1,
                       autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &auto2,
                       int iterations);
  void initFillWithParticles(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
                             std::array<unsigned long, 3> particlesPerDim);

  static std::array<double, 3> nextPosition(std::array<double, 3> position, std::array<double, 3> force,
                                            std::array<double, 3> velocity, double particle_delta_t);
  static std::array<double, 3> nextVelocity(std::array<double, 3> velocity, std::array<double, 3> force,
                                            std::array<double, 3> oldf, double particle_delta_t);

  void Pos_and_Velo_Test(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autopas,
                         size_t numberOfParticles, int iterations);

 protected:
  double epsilon;
  double sigma;
  double cutoff;
  std::array<double, 3> boxmin;
  std::array<double, 3> boxmax;
  ParticleClassLibrary PCL;
  autopas::LJFunctor<PrintableMolecule, autopas::ParticleCell<PrintableMolecule>, autopas::FunctorN3Modes::Both, true>
      functor;
};