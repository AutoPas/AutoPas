/**
 * @file OctreeTest.h
 * @author C.Menges
 * @date 11.07.2019
 */

#include "OctreeTest.h"
#include <gmock/gmock-generated-matchers.h>

TEST_F(OctreeTest, testToString) {
  std::vector<FPCell> cells;
  cells.resize(8);
  std::array<double, 3> boxMin({0.0, 0.0, 0.0});
  std::array<double, 3> boxMax({8.0, 8.0, 8.0});
  autopas::Octree<Particle, FPCell> octree(cells, boxMin, boxMax);
  octree.init({2ul, 2ul, 2ul});
  // EXPECT_EQ("0", static_cast<std::string>(octree));
}

TEST_F(OctreeTest, testSplit) {
  std::vector<FPCell> cells;
  cells.resize(64);
  std::array<double, 3> boxMin({0.0, 0.0, 0.0});
  std::array<double, 3> boxMax({8.0, 8.0, 8.0});
  autopas::Octree<Particle, FPCell> octree(cells, boxMin, boxMax);
  octree.init({4ul, 4ul, 4ul});
  octree.setMaxElements(1);

  auto &cell = octree.getContainingCell({1.0, 1.0, 1.0});
  EXPECT_EQ(&cells[0], &cell);
  Particle defaultParticle({1.0, 4.0, 1.0}, {1.0, 1.0, 1.0}, 0);
  cell.addParticle(defaultParticle);
  for (int i = 1; i < 4; i++) {
    defaultParticle.setID(i);
    cell.addParticle(defaultParticle);
  }
  octree.update();
  /*octree.apply([](autopas::internal::OctreeNode<Particle, FPCell> &node){
      std::cout << "Index: " << dynamic_cast<autopas::internal::OctreeExternalNode<Particle, FPCell>&>(node).getIndex()
     << std::endl; for(auto n : dynamic_cast<autopas::internal::OctreeExternalNode<Particle, FPCell>&>(node)._neighbors)
          std::cout << n.first << std::endl;}, autopas::internal::ExecutionPolicy::seq);*/
  EXPECT_EQ(octree.getSize(), 4);

  octree.update();
  /* octree.apply([](autopas::internal::OctreeNode<Particle, FPCell> &node){
       std::cout << "Index: " << dynamic_cast<autopas::internal::OctreeExternalNode<Particle, FPCell>&>(node).getIndex()
     << std::endl; for(auto n : dynamic_cast<autopas::internal::OctreeExternalNode<Particle, FPCell>&>(node)._neighbors)
           std::cout << n.first << std::endl;}, autopas::internal::ExecutionPolicy::seq);*/
  auto lastUpdate = static_cast<std::string>(octree);
  octree.update();
  auto update = static_cast<std::string>(octree);
  EXPECT_EQ(octree.getSize(), 4);
  EXPECT_EQ(lastUpdate, update);
  std::cout << static_cast<std::string>(octree);
}