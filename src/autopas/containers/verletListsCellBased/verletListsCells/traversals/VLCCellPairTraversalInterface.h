/**
 * @file VLCCellPairTraversalInterface.h
 * @author tirgendetwas
 * @date 11.01.21
 */

#pragma once
namespace autopas {
template <class Particle>
/**
 * TODO
 * @tparam Particle
 */
class VLCCellPairTraversalInterface {
 public:
  void setVerletList(VLCCellPairNeighborList<Particle> &verlet) { _cellPairVerletList = &verlet; }

 protected:
  VLCCellPairNeighborList<Particle> *_cellPairVerletList;
};
}  // namespace autopas
