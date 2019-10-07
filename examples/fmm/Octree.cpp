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

OctreeNode::OctreeNode(int x, int y, int z, int size, double cellSize, Octree *tree, OctreeNode *parent) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->tree = tree;
  this->size = size;
  this->cellSize = cellSize;
  this->parent = parent;
  if (parent == nullptr) {
    this->depth = 0;
    this->totalSize = size;
  } else {
    this->depth = parent->depth + 1;
    this->totalSize = parent->totalSize;
  }

  child = std::vector<std::unique_ptr<OctreeNode>>(8);

  if (size > 1) {
    int newSize = size / 2;

    // Create 8 children.
    child[0] = std::make_unique<OctreeNode>(x, y, z, newSize, cellSize, tree, this);
    child[1] = std::make_unique<OctreeNode>(x + newSize, y, z, newSize, cellSize, tree, this);
    child[2] = std::make_unique<OctreeNode>(x, y + newSize, z, newSize, cellSize, tree, this);
    child[3] = std::make_unique<OctreeNode>(x + newSize, y + newSize, z, newSize, cellSize, tree, this);
    child[4] = std::make_unique<OctreeNode>(x, y, z + newSize, newSize, cellSize, tree, this);
    child[5] = std::make_unique<OctreeNode>(x + newSize, y, z + newSize, newSize, cellSize, tree, this);
    child[6] = std::make_unique<OctreeNode>(x, y + newSize, z + newSize, newSize, cellSize, tree, this);
    child[7] = std::make_unique<OctreeNode>(x + newSize, y + newSize, z + newSize, newSize, cellSize, tree, this);

  } else {
    // cont = tree->getAutoPasCont(x,y,z);
    this->lowCorner = std::array<double, 3>({x * cellSize, y * cellSize, z * cellSize});
    this->highCorner = std::array<double, 3>({(x + 1.0) * cellSize, (y + 1.0) * cellSize, (z + 1.0) * cellSize});
    this->haloLowCorner = std::array<double, 3>({(x - 1.0) * cellSize, (y - 1.0) * cellSize, (z - 1.0) * cellSize});
    this->haloHighCorner = std::array<double, 3>({(x + 2.0) * cellSize, (y + 2.0) * cellSize, (z + 2.0) * cellSize});
    // Create AutoPas object at leaf.
    cont = std::make_unique<AutoPasCont>();

    cont->setAllowedContainers(std::set<autopas::ContainerOption>{autopas::ContainerOption::directSum});

    cont->setBoxMin(lowCorner);
    cont->setBoxMax(highCorner);

    cont->init();
  }
  tree->setCell(this->depth, x, y, z, this);
}

void OctreeNode::initNeighbourList() {
  neighbourList = std::set<OctreeNode *>();

  for (int i = std::max(0, x - size); i < std::min(totalSize, x + size + 1); i += size) {
    for (int j = std::max(0, y - size); j < std::min(totalSize, y + size + 1); j += size) {
      for (int k = std::max(0, z - size); k < std::min(totalSize, z + size + 1); k += size) {
        neighbourList.insert(tree->getCell(depth, i, j, k));
      }
    }
  }

  if (size != 1) {
    for (int c = 0; c < 8; ++c) {
      child[c]->initNeighbourList();
    }
  }
}

void OctreeNode::initInteractionList() {
  // interactionList = parent->neighbourList / this->neighbourList
  interactionList = std::set<OctreeNode *>();
  if (parent != nullptr) {
    for (OctreeNode *parentNeighbour : *parent->getNeighbourList()) {
      for (int c = 0; c < 8; ++c) {
        auto childOfParent = parentNeighbour->getChild(c);
        int nX = childOfParent->x;
        int nY = childOfParent->y;
        int nZ = childOfParent->z;

        if (neighbourList.find(tree->getCell(depth, nX, nY, nZ)) == neighbourList.end()) {
          interactionList.insert(childOfParent);
        }
      }
    }
  }

  if (size != 1) {
    for (int c = 0; c < 8; ++c) {
      child[c]->initInteractionList();
    }
  }
}