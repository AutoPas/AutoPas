//
// Created by TinaVl on 1/11/2021.
//

#pragma once
namespace autopas
{
template <class Particle>
class VLCCellPairTraversalInterface
{
 public:
  void setVerletList(VLCCellPairNeighborList<Particle> &verlet) { _cellPairVerletList = &verlet; }

 protected:
  VLCCellPairNeighborList<Particle> *_cellPairVerletList;
};
}
