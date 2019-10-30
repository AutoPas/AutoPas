/**
 * @file Octree.h
 * @date 14.09.19
 * @author Joachim Marin
 */

#pragma once

#include <complex>
#include <memory>
#include <string>
#include <vector>
#include "FmmParticle.h"
#include "OctreeNode.h"
#include "autopas/AutoPas.h"

using CellMatrix2D = std::vector<std::vector<OctreeNode *>>;
using CellMatrix3D = std::vector<CellMatrix2D>;
using CellMatrix4D = std::vector<CellMatrix3D>;

class Octree {
 public:
  explicit Octree(int size, double cellSize, int contSize);

  std::shared_ptr<AutoPasCont> getAutoPasCont(int x, int y, int z) {
    return contMatrix[x / contSize][y / contSize][z / contSize];
  }

  OctreeNode *getCell(int depth, int x, int y, int z);

  void setCell(int depth, int x, int y, int z, OctreeNode *leaf);

  OctreeNode *getRoot() { return &(*root); }

  int getHeight() { return height; }
  int getSize() { return size; }
  double getCellSize() { return cellSize; }

 private:
  std::unique_ptr<OctreeNode> root;
  CellMatrix4D cellMatrix;
  std::vector<std::vector<std::vector<std::shared_ptr<AutoPasCont>>>> contMatrix;
  int size;
  int height;
  int contSize;
  double cellSize;
};
