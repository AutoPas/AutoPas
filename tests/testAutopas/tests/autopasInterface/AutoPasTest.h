/**
 * @file AutoPasTest.h
 * @author seckler
 * @date 29.05.18
 */

#pragma once

#include <gtest/gtest.h>

#include "NumThreadGuard.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;

class AutoPasTest : public testing::Test {
 public:
  AutoPasTest() {
    autoPas.setBoxMin({0., 0., 0.});
    autoPas.setBoxMax({10., 10., 10.});
    autoPas.setCutoff(1.);
    autoPas.init();
  }

 protected:
  NumThreadGuard _numThreadGuard{autopas::autopas_get_max_threads(), TUNED_THREADS};
  void expectedParticles(size_t expectedOwned, size_t expectedHalo);

  autopas::AutoPas<Molecule> autoPas;
};
