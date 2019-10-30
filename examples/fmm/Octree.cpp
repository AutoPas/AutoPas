/**
 * @file Octree.cpp
 * @date 14.09.19
 * @author Joachim Marin
 */

#include "Octree.h"

Octree::Octree(int size, double cellSize, int contSize) {
  this->size = size;
  this->cellSize = cellSize;
  this->height = int(std::lround(std::log2(size)));
  this->contSize = contSize;

  int numContainers = size / contSize;

  // Init cellMatrix.
  cellMatrix =
      CellMatrix4D(height + 1, CellMatrix3D(size, CellMatrix2D(size, std::vector<OctreeNode *>(size, nullptr))));

  contMatrix = std::vector<std::vector<std::vector<std::shared_ptr<AutoPasCont>>>>(
      numContainers, std::vector<std::vector<std::shared_ptr<AutoPasCont>>>(
                         numContainers, std::vector<std::shared_ptr<AutoPasCont>>(numContainers)));

  double contCellSize = cellSize * contSize;
  for (int i = 0; i < numContainers; ++i) {
    for (int j = 0; j < numContainers; ++j) {
      for (int k = 0; k < numContainers; ++k) {
        contMatrix[i][j][k] = std::make_shared<AutoPasCont>();

        contMatrix[i][j][k]->setAllowedContainers(
            std::set<autopas::ContainerOption>{autopas::ContainerOption::directSum});

        contMatrix[i][j][k]->setBoxMin({i * contCellSize, j * contCellSize, k * contCellSize});
        contMatrix[i][j][k]->setBoxMax({(1.0 + i) * contCellSize, (1.0 + j) * contCellSize, (1.0 + k) * contCellSize});

        contMatrix[i][j][k]->init();
      }
    }
  }

  // Create tree.
  root = std::make_unique<OctreeNode>(0, 0, 0, size, cellSize, this, nullptr);
  root->initNeighbourList();
  root->initInteractionList();
}

OctreeNode *Octree::getCell(int depth, int x, int y, int z) { return cellMatrix[depth][x][y][z]; }

void Octree::setCell(int depth, int x, int y, int z, OctreeNode *leaf) { cellMatrix[depth][x][y][z] = leaf; }

