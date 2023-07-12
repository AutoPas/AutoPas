/**
 * @file ParticleSerializationToolsTest.h
 * @author J. Körner
 * @date 19/05/21
 */
#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "src/ParticleSerializationTools.h"
#include "src/TypeDefinitions.h"

/**
 * Test class for the ParticleSerializationTools.
 */
class ParticleSerializationToolsTest : public AutoPasTestBase {
 public:
  /**
   * Constructor.
   */
  ParticleSerializationToolsTest();

 protected:
  /**
   * Molecule used for testing.
   */
  ParticleType _molecule;
};
