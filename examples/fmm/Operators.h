/**
 * @file Operators.h
 * @date 15.09.19
 * @author Joachim Marin
 */
#pragma once

#include <cmath>
#include <complex>
#include <iostream>
#include "FmmParticle.h"
#include "Math3D.h"
#include "Octree.h"

class Operators {
 public:
  Operators() = default;

  static void P2M(OctreeNode *leaf);

  static void M2M(OctreeNode *parent);

  static void M2L(OctreeNode *node);

  static void L2L(OctreeNode *parent);

  static void L2P(OctreeNode *leaf);
};
