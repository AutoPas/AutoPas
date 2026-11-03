/**
 * @file OpenMPConfiguratorTest.h
 * @author S. J. Newcome
 * @date 28/09/2026
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/OpenMPKindOption.h"
#include "autopas/tuning/Configuration.h"

/**
 * Tests that the OpenMP schedule kinds Configuration::hasCompatibleValues() accepts can actually be set with the
 * compiler and OpenMP runtime AutoPas is built with, and that the others fail when trying to set them.
 *
 * Parameterized over all OpenMPKindOptions.
 */
class OpenMPConfiguratorTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::OpenMPKindOption> {
 public:
  /**
   * The configuration all tests use, which only varies in the OpenMP schedule kind (and in the chunk size where
   * relevant).
   * @param ompKind The OpenMP schedule kind.
   * @return The configuration.
   */
  static autopas::Configuration getConfiguration(autopas::OpenMPKindOption ompKind);
};