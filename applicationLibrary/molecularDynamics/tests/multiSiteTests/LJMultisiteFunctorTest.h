/**
 * @file LJMultisiteFunctorTest.h
 * @author S. Newcome
 * @date 12/05/2022
 */

#pragma once

#include <gtest/gtest.h>

#include "LJParentMultiSiteFunctorTest.h"

#include "autopas/AutoPasDecl.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

/**
 * Test class for LJMultisiteFunctor
 */
class LJMultisiteFunctorTest : public LJParentMultiSiteFunctorTest<mdLib::MultisiteMoleculeLJ> {};
