/**
 * @file IteratorBehaviorTest.cpp
 * @author F. Gratl
 * @date 13.12.22
 */

#include "autopas/options/IteratorBehavior.h"
#include "autopas/particles/OwnershipState.h"

static_assert(static_cast<unsigned int>(autopas::IteratorBehavior::owned) ==
                  static_cast<unsigned int>(autopas::OwnershipState::owned),
              "IteratorBehavior::owned and OwnershipState::owned should map to the same value");
static_assert(static_cast<unsigned int>(autopas::IteratorBehavior::halo) ==
                  static_cast<unsigned int>(autopas::OwnershipState::halo),
              "IteratorBehavior::halo and OwnershipState::halo should map to the same value");
