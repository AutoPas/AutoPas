/**
 * @file main.cpp
 * @date 14.09.19
 * @author Joachim Marin
 */

#include <iostream>
#include <Octree.h>
#include <Operators.h>
#include <FmmParticle.h>


void upwardPassRec(OctreeNode *node) {
  if (node->isLeaf()) {
    Operators::P2M(node);
  } else {
    for (int i = 0; i < 8; ++i) {
      auto child = node->getChild(i);
      upwardPassRec(child);
    }
    Operators::M2M(node);
  }
}

void upwardPass(Octree *tree) {
  upwardPassRec(tree->getRoot());
}

void downwardPassRec1(OctreeNode *node) {
  Operators::M2L(node);
  if (!node->isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node->getChild(i);
      downwardPassRec1(child);
    }
  }
}

void downwardPassRec2(OctreeNode *node) {
  Operators::L2L(node);
  if (!node->isLeaf()) {
    for (int i = 0; i < 8; ++i) {
      auto child = node->getChild(i);
      downwardPassRec2(child);
    }
  } else {
    Operators::L2P(node);
  }
}

void downwardPass(Octree *tree) {
  downwardPassRec1(tree->getRoot());
  downwardPassRec2(tree->getRoot());
}

int main(int argc, char **argv) {

  std::cout << "test" << std::endl;

  Octree tree(4, 2.0);

  auto corner1 = tree.getCell(tree.getHeight(), 0, 0, 0)->getContainer();
  auto cell3 = tree.getCell(tree.getHeight(), 1, 1, 1)->getContainer();

  auto corner2 = tree.getCell(tree.getHeight(), 3, 3, 3)->getContainer();

  FmmParticle particle({1, 1, 1}, {0, 0, 0}, 0, 17.0);
  corner1->addParticle(particle);

  particle = FmmParticle({0.5, 0.5, 0.5}, {0, 0, 0}, 0, 25.0);
  corner1->addParticle(particle);
  particle = FmmParticle({1.5, 1.5, 1.5}, {0, 0, 0}, 0, 25.0);
  corner1->addParticle(particle);


  particle = FmmParticle({7, 7, 7}, {0, 0, 0}, 1, 42.0);
  corner2->addParticle(particle);


  particle = FmmParticle({3, 3, 3}, {0, 0, 0}, 1, 1000);
  cell3->addParticle(particle);


  upwardPass(&tree);
  downwardPass(&tree);






  // Printed Result:
  // result = 4.04145
  // result = 4.15036
  // result = 5.06055
  // result = 6.06218
  // result = 150.369

  // Expected Result:
  // (1/1/1):
  // p = 4.04       good
  //
  // (0.5/0.5/0.5)
  // p = 3.73       off by 0.42
  //
  // (1.5/1.5/1.5)
  // p = 4.41       off by 0.65
  //
  // (3,3,3)
  // p = 6.06       good
  //
  // (7,7,7)
  // p += 1.64
  // p += 2.22
  // p += 2.62
  // p += 144.30
  // p = 150.78     off by 0.41

  return EXIT_SUCCESS;
}