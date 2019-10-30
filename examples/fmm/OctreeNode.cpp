/**
 * @file OctreeNode.cpp
 * @date 30.10.19
 * @author Joachim Marin
 */

#include "OctreeNode.h"
#include "Octree.h"


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