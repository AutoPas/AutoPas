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
#include "AdaptiveOctree.h"

class Operators {
 public:
  explicit Operators(int orderOfExpansion);

  /**
   * Calculates the coefficients M(m,n) for the leaf based on the particles it contains.
   * See formula 5.16.
   * @param The cell whose M(m,n) coefficients are calculated.
   */
  void P2M(AdaptiveOctreeNode &leaf);

  /**
   * Calculates the coefficients M(m,n) for the node based on its children's coefficients.
   * See formula 5.22.
   * @param The cell whose M(m,n) coefficients are calculated.
   */
  void M2M(AdaptiveOctreeNode &parent);

  /**
   * Calculates the coefficients L(m,n) for the node based on the M(m,n) coefficients of the cells in the node's
   * interaction list.
   * See formula 5.26.
   * @param node The cell whose coefficients are calculated.
   */
  void M2L(AdaptiveOctreeNode &node);

  /**
   * Calculates the coefficients L(m,n) for the node's children based on the coefficients of this cell.
   * @param parent The cell whose children's coefficients are calculated.
   * See formula 5.30.
   */
  void L2L(AdaptiveOctreeNode &parent);

  /**
   * Evaluates the expansion for every particle and stores the result in the particle.
   * @param leaf The cell containing the particles that are evaluated.
   * See formula 5.28.
   */
  void L2P(AdaptiveOctreeNode &leaf);

 private:
  int orderOfExpansion;
  // i^(|k-m|-|k|-|m|)
  std::vector<std::vector<Complex>> powerM2L;
};
