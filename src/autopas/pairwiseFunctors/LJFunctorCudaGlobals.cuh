/**
 * @file LJFunctorCuda.cuh
 *
 * @date 26.4.2019
 * @author jspahl
 */

#pragma once

#include <array>
#include "autopas/pairwiseFunctors/FunctorCuda.cuh"
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#include "cuda_runtime.h"

namespace autopas {

/**
 * Stores all constants needed for the calculation
 * @tparam floatType of constants
 */
template <typename floatType>
class LJFunctorGlobalsConstants : public FunctorCudaConstants<floatType> {
 public:
  LJFunctorGlobalsConstants() {}
  LJFunctorGlobalsConstants(floatType csq, floatType ep24, floatType sqs, floatType sh6,
                            std::array<floatType, 3>& lowCorner, std::array<floatType, 3>& highCorner)
      : cutoffsquare(csq), epsilon24(ep24), sigmasquare(sqs), shift6(sh6) {
    boxMin.x = lowCorner[0];
    boxMin.y = lowCorner[1];
    boxMin.z = lowCorner[2];

    boxMax.x = highCorner[0];
    boxMax.y = highCorner[1];
    boxMax.z = highCorner[2];
  }
  floatType cutoffsquare;
  floatType epsilon24;
  floatType sigmasquare;
  floatType shift6;
  typename vec3<floatType>::Type boxMin;
  typename vec3<floatType>::Type boxMax;
};

/**
 * Stores all pointers to the device Memory SoAs as needed by the LJ Functor
 * @tparam floatType of all vectors
 */
template <typename floatType>
class LJFunctorCudaGlobalsSoA : public FunctorCudaSoA<floatType> {
 public:
  /**
   * Constructor for only positions
   * @param size Number of particles
   * @posX x positions of the particles
   * @posY y positions of the particles
   * @posZ z positions of the particles
   */
  LJFunctorCudaGlobalsSoA(unsigned int size, floatType* posX, floatType* posY, floatType* posZ, floatType* globals)
      : _size(size),
        _posX(posX),
        _posY(posY),
        _posZ(posZ),
        _forceX(NULL),
        _forceY(NULL),
        _forceZ(NULL),
        _globals(globals) {}

  /**
   * Constructor for only positions
   * @param size Number of particles
   * @posX x positions of the particles
   * @posY y positions of the particles
   * @posZ z positions of the particles
   * @forceX x forces of the particles
   * @forceY y forces of the particles
   * @forceZ z forces of the particles
   */
  LJFunctorCudaGlobalsSoA(unsigned int size, floatType* posX, floatType* posY, floatType* posZ, floatType* forceX,
                          floatType* forceY, floatType* forceZ, floatType* globals)
      : _size(size),
        _posX(posX),
        _posY(posY),
        _posZ(posZ),
        _forceX(forceX),
        _forceY(forceY),
        _forceZ(forceZ),
        _globals(globals) {}

  /**
   * CopyConstructor
   * @param obj other object
   */
  LJFunctorCudaGlobalsSoA(const LJFunctorCudaGlobalsSoA& obj)
      : _size(obj._size),
        _posX(obj._posX),
        _posY(obj._posY),
        _posZ(obj._posZ),
        _forceX(obj._forceX),
        _forceY(obj._forceY),
        _forceZ(obj._forceZ),
        _globals(obj._globals) {}

  unsigned int _size;
  floatType* _posX;
  floatType* _posY;
  floatType* _posZ;
  floatType* _forceX;
  floatType* _forceY;
  floatType* _forceZ;
  floatType* _globals;
};

template <typename floatType>
class LJFunctorCudaGlobalsWrapper : public CudaWrapperInterface<floatType> {
 public:
  LJFunctorCudaGlobalsWrapper() { _num_threads = 64; }
  virtual ~LJFunctorCudaGlobalsWrapper() {}

  void setNumThreads(int num_threads) override { _num_threads = num_threads; }

  void loadConstants(FunctorCudaConstants<floatType>* constants) override;

  void SoAFunctorNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, cudaStream_t stream = 0) override;
  void SoAFunctorNoN3PairWrapper(FunctorCudaSoA<floatType>* cell1Base, FunctorCudaSoA<floatType>* cell2Base,
                                 cudaStream_t stream) override;

  void SoAFunctorN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, cudaStream_t stream = 0) override;
  void SoAFunctorN3PairWrapper(FunctorCudaSoA<floatType>* cell1Base, FunctorCudaSoA<floatType>* cell2Base,
                               cudaStream_t stream) override;

  void LinkedCellsTraversalNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, unsigned int reqThreads,
                                       unsigned int cids_size, unsigned int* cids, unsigned int cellSizes_size,
                                       size_t* cellSizes, cudaStream_t stream) override;

  void LinkedCellsTraversalN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, unsigned int reqThreads,
                                     unsigned int cids_size, unsigned int* cids, unsigned int cellSizes_size,
                                     size_t* cellSizes, cudaStream_t stream) override;

  void CellVerletTraversalNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, unsigned int ncells,
                                      unsigned int clusterSize, unsigned int others_size, unsigned int* other_ids,
                                      cudaStream_t stream) override;

  void CellVerletTraversalN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, unsigned int ncells, unsigned int clusterSize,
                                    unsigned int others_size, unsigned int* other_ids, cudaStream_t stream) override;

  void loadLinkedCellsOffsets(unsigned int offsets_size, int* offsets) override;

 private:
  int numRequiredBlocks(int n) { return ((n - 1) / _num_threads) + 1; }

  int _num_threads;
};

}  // namespace autopas
