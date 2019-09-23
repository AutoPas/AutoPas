/**
 * @file Octree.h
 * @date 14.09.19
 * @author Joachim Marin
 */

#pragma once

#include <vector>
#include <autopas/AutoPas.h>
#include <memory>
#include <string>
#include <FmmParticle.h>
#include <complex>

constexpr int orderOfExpansion = 6;

using AutoPasCont = autopas::AutoPas<FmmParticle, autopas::FullParticleCell<FmmParticle>>;
using Complex = std::complex<double>;
using ComplexMatrix = std::vector<std::vector<Complex>>;


class OctreeNode;

using CellMatrix2D = std::vector<std::vector<OctreeNode *>>;
using CellMatrix3D = std::vector<CellMatrix2D>;
using CellMatrix4D = std::vector<CellMatrix3D>;

class Octree {
public:
    explicit Octree(int size, double cellSize);

    OctreeNode *getCell(int depth, int x, int y, int z);

    void setCell(int depth, int x, int y, int z, OctreeNode *leaf);

    OctreeNode *getRoot() {
      return &(*root);
    }

    int getHeight() {
      return height;
    }

private:
    std::unique_ptr<OctreeNode> root;
    CellMatrix4D cellMatrix;
    int size;
    int height;

};

class OctreeNode {
public:
    explicit OctreeNode(int x, int y, int z, int size, double cellSize, Octree *tree, OctreeNode *parent);

    bool isLeaf() {
      return size == 1;
    }

    OctreeNode *getChild(int i) {
      return child[i].get();
    }

    AutoPasCont *getContainer() {
      return &(*cont);
    }

    std::array<double, 3> getCenter() {
      double cX, cY, cZ;
      cX = (this->x + 0.5 * this->size) * this->cellSize;
      cY = (this->y + 0.5 * this->size) * this->cellSize;
      cZ = (this->z + 0.5 * this->size) * this->cellSize;
      return std::array<double, 3>({cX, cY, cZ});
    }

    ComplexMatrix fmmM;
    ComplexMatrix fmmL;


    Complex getM(int m, int n) {
      int offset = m <= 0 ? 0 : -1;
      int indexM = 2 * std::abs(m) + offset;
      if (indexM < 0 || indexM >= static_cast<int>(fmmM.size())) {
        return 0.0;
      }
      if (n < 0 || n >= static_cast<int> (fmmM[0].size())) {
        return 0.0;
      }
      return fmmM.at(indexM).at(n);
    }

    void setM(int m, int n, Complex value) {
      int offset = m <= 0 ? 0 : -1;
      fmmM.at(2 * std::abs(m) + offset).at(n) = value;

      if (value != 0.0) {
        isZeroM = false;
      }
    }

    Complex getL(int m, int n) {
      int offset = m <= 0 ? 0 : -1;
      int indexM = 2 * std::abs(m) + offset;
      if (indexM < 0 || indexM >= static_cast<int>(fmmL.size())) {
        return 0.0;
      }
      if (n < 0 || n >= static_cast<int> (fmmL[0].size())) {
        return 0.0;
      }
      return fmmL.at(indexM).at(n);
    }

    void setL(int m, int n, Complex value) {
      int offset = m <= 0 ? 0 : -1;
      fmmL.at(2 * std::abs(m) + offset).at(n) = value;

      if (value != 0.0) {
        isZeroL = false;
      }
    }

    void initNeighbourList();

    void initInteractionList();

    std::set<OctreeNode *> *getNeighbourList() {
      return &neighbourList;
    }

    std::set<OctreeNode *> *getInteractionList() {
      return &interactionList;
    }

    int getDepth() {
      return depth;
    }

    OctreeNode *getParent() {
      return parent;
    }

    bool getIsZeroL() {
      return isZeroL;
    }

    bool getIsZeroM() {
      return isZeroM;
    }

private:
    Octree *tree;
    std::vector<std::unique_ptr<OctreeNode>> child;
    OctreeNode *parent;
    std::unique_ptr<AutoPasCont> cont;
    int x, y, z;
    double cellSize;
    int size;
    int depth;
    int totalSize;
    std::set<OctreeNode *> neighbourList;
    std::set<OctreeNode *> interactionList;
    bool isZeroL = true;
    bool isZeroM = true;

};

