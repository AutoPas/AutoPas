/**
 * @file ParticleAttributes.h
 * @author J. Körner
 * @date 13.05.2021
 */
#pragma once

#include <array>

#include "autopas/particles/OwnershipState.h"

/*
 * A struct containing all properties of autopas::MoleculeLJ<double>.
 * This can be used to align the attributes of a particle in memory to make serialization and deserialization easier.
 */
struct ParticleAttributes {
  /**
   * The position attribute
   */
  std::array<double, 3> position;

  /**
   * The velocity attribute
   */
  std::array<double, 3> velocity;

  /**
   * The force attribute
   */
  std::array<double, 3> force;

  /**
   * The id attribute
   */
  unsigned long id;

  /**
   * The ownershipState attribute
   */
  autopas::OwnershipState ownershipState;

  /**
   * The type attribute
   */
  size_t typeId;

  /**
   * The oldForce attribute
   */
  std::array<double, 3> oldForce;
};
