/**
 * @file TraversalRaceConditionTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include <vector>
#include "../../examples/md/mdutils.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/containers/cellPairTraversals/SlicedTraversal.h"

class TraversalRaceConditionTest : public AutoPasTestBase {
 public:
  TraversalRaceConditionTest() = default;

  ~TraversalRaceConditionTest() override = default;

  void fillWithParticles(autopas::AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autoPas,
                         std::array<size_t, 3> particlesPerDim);

  /*
   * Simple AoS only functor which repulses paritcles from each other with a
   * constant force of 1.
   */
  class SimpleFunctor : public autopas::Functor<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> {
   public:
    typedef PrintableMolecule Particle;
    typedef PrintableMolecule::SoAArraysType SoAArraysType;
    typedef autopas::FullParticleCell<PrintableMolecule> ParticleCell;

    void AoSFunctor(PrintableMolecule &i, PrintableMolecule &j, bool newton3) override {
      auto coordsI = i.getR();
      auto coordsJ = j.getR();

      std::array<double, 3> dr = ArrayMath::sub(coordsI, coordsJ);
      double dr2 = ArrayMath::dot(dr, dr);

      if (dr2 > CUTOFFSQUARE) return;

      std::array<double, 3> f;

      for (int dim = 0; dim < 3; ++dim) {
        if (coordsI[dim] < coordsJ[dim]) {
          f[dim] = -1;
        } else if (coordsI[dim] > coordsJ[dim]) {
          f[dim] = 1;
        } else {
          f[dim] = 0;
        }
      }

      i.addF(f);
      j.subF(f);
    }

    AUTOPAS_FUNCTOR_SOAEXTRACTOR(, , , );

    AUTOPAS_FUNCTOR_SOALOADER(, , , );

   private:
    // in a grid with separation 1 this includes all neighbors with a Chebyshev
    // distance of 1
    constexpr static double CUTOFFSQUARE = 3;
  };
};
