/**
 * @file MulticenteredLJFunctorTest.cpp
 * @author S. Newcome
 * @date 16/05/2022
 */

#pragma once

#include <gtest/gtest.h>

#include "MulticenteredLJFunctorTest.h"

void MulticenteredLJFunctorTest::testForceCalculation(multisiteMolecule molA, multisiteMolecule molB, double cutoff) {
  // check size of site property vectors match
  EXPECT_EQ(molA.sitePositions.size(),molA.siteEpsilons.size()) << "molA: Number of epsilons does not match number of site positions";
  EXPECT_EQ(molA.sitePositions.size(),molA.siteSigmas.size()) << "molA: Number of sigmas does not match number of site positions";
  EXPECT_EQ(molB.sitePositions.size(),molB.siteEpsilons.size()) << "molB: Number of epsilons does not match number of site positions";
  EXPECT_EQ(molB.sitePositions.size(),molB.siteSigmas.size()) << "molB: Number of sigmas does not match number of site positions";
  // create + fill PPL
  ParticlePropertiesLibrary PPL(cutoff);
  std::vector<size_t> molASiteIds;
  for (int i = 0; i < molA.sitePositions.size(); ++i) {
    PPL.addSiteType(i,molA.siteEpsilons[i],molA.siteSigmas[i],1);
    molASiteIds[i] = i;
  }
  PPL.addMolType(0, molASiteIds, molA.sitePositions);

  std::vector<size_t> molBSiteIds;
  for (int i = 0; i < molB.sitePositions.size(); ++i) {
    PPL.addSiteType(i+molA.sitePositions.size(),molB.siteEpsilons[i],molB.siteSigmas[i],1);
    molBSiteIds[i] = i+molA.sitePositions.size();
  }
  PPL.addMolType(1, molBSiteIds, molB.sitePositions);

  // create molecules
  autopas::MulticenteredMoleculeLJ molAParticle, molBParticle;
  molAParticle.setR()

}