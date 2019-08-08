/**
 * @file VerletClusterCellsTraversal.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <algorithm>
#include <vector>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/VerletClusterTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/CudaDeviceVector.h"
#include "autopas/utils/SoAView.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas {

/**
 * This Traversal is used to interact all clusters in VerletClusterCluster Container
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class VerletClusterCellsTraversal : public CellPairTraversal<ParticleCell>,
                                    public VerletClusterTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for the VerletClusterClusterTraversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  VerletClusterCellsTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>({1, 1, 1}),
        _functor(pairwiseFunctor),
        _neighborMatrixDim(nullptr),
        _clusterSize(nullptr) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClusterCellsTraversal; }

  bool isApplicable() const override {
    // TODO enable when functors use owned pointers correctly for soas and global calculation
    if (DataLayout == DataLayoutOption::soa) {
      return false;
    }
    if (DataLayout == DataLayoutOption::cuda) {
      int nDevices = 0;
#if defined(AUTOPAS_CUDA)
      cudaGetDeviceCount(&nDevices);
      if (not _functor->getCudaWrapper()) return false;
#endif
      return nDevices > 0;
    } else
      return true;
  }

  bool getUseNewton3() const override { return useNewton3; }

  DataLayoutOption getDataLayout() const override { return DataLayout; }

  std::tuple<TraversalOption, DataLayoutOption, bool> getSignature() override {
    return std::make_tuple(TraversalOption::verletClusterCellsTraversal, DataLayout, useNewton3);
  }

  void setVerletListPointer(unsigned int *clusterSize,
                            std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> *neighborCellIds,
                            size_t *neighborMatrixDim, utils::CudaDeviceVector<unsigned int> *neighborMatrix) override {
    _clusterSize = clusterSize;
    _neighborCellIds = neighborCellIds;
    _neighborMatrixDim = neighborMatrixDim;
    _neighborMatrix = neighborMatrix;
  }

  void rebuildVerlet(const std::array<unsigned long, 3> &dims, std::vector<ParticleCell> &cells,
                     std::vector<std::vector<std::array<double, 6>>> &boundingBoxes, int interactionCellRadius,
                     double distance) override {
    this->_cellsPerDimension = dims;

    const size_t cellsSize = cells.size();
    _neighborCellIds->clear();
    _neighborCellIds->resize(cellsSize, {});

    for (size_t i = 0; i < cellsSize; ++i) {
      auto pos = utils::ThreeDimensionalMapping::oneToThreeD(i, this->_cellsPerDimension);
      for (int x = -interactionCellRadius; x <= interactionCellRadius; ++x) {
        if (0 <= (pos[0] + x) and (pos[0] + x) < this->_cellsPerDimension[0]) {
          for (int y = -interactionCellRadius; y <= interactionCellRadius; ++y) {
            if (0 <= (pos[1] + y) and (pos[1] + y) < this->_cellsPerDimension[1]) {
              // add neighbors
              auto other = utils::ThreeDimensionalMapping::threeToOneD(pos[0] + x, pos[1] + y, (unsigned long)0,
                                                                       this->_cellsPerDimension);
              if (useNewton3 and other > i) {
                continue;
              }
              // own clusters
              for (size_t ownClusterId = 0; ownClusterId < boundingBoxes[i].size(); ++ownClusterId) {
                (*_neighborCellIds)[i].resize(boundingBoxes[i].size(), {});
                const std::array<double, 6> ownBox = boundingBoxes[i][ownClusterId];

                auto start = std::find_if(boundingBoxes[other].begin(), boundingBoxes[other].end(),
                                          [this, ownBox, distance](const std::array<double, 6> &otherbox) {
                                            return getMinDist(ownBox, otherbox) < distance;
                                          });
                auto end = std::find_if(start, boundingBoxes[other].end(),
                                        [this, ownBox, distance](const std::array<double, 6> &otherbox) {
                                          return getMinDist(ownBox, otherbox) > distance;
                                        });

                const size_t size = end - start;

                if (start != end) {
                  (*_neighborCellIds)[i][ownClusterId].reserve(size);
                  auto indexStart = start - boundingBoxes[other].begin();
                  if (other == i) {
                    for (size_t k = 0; k < size; ++k) {
                      if (useNewton3) {
                        if (indexStart + k > ownClusterId) {
                          (*_neighborCellIds)[i][ownClusterId].push_back(std::make_pair(other, indexStart + k));
                        }
                      } else {
                        if (indexStart + k != ownClusterId) {
                          (*_neighborCellIds)[i][ownClusterId].push_back(std::make_pair(other, indexStart + k));
                        }
                      }
                    }
                  } else {
                    for (size_t k = 0; k < size; ++k) {
                      (*_neighborCellIds)[i][ownClusterId].push_back(std::make_pair(other, indexStart + k));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (DataLayout == DataLayoutOption::cuda) {
      size_t neighborMatrixDim = 0;
      for (auto &cell : *_neighborCellIds) {
        for (auto &cluster : cell) {
          neighborMatrixDim = std::max(neighborMatrixDim, cluster.size());
        }
      }

      ++neighborMatrixDim;
      if (not useNewton3) {
        ++neighborMatrixDim;
      }
      *_neighborMatrixDim = neighborMatrixDim;

      std::vector<size_t> cellSizePartSums(cellsSize + 1, 0);
      for (size_t i = 0; i < cellsSize; ++i) {
        cellSizePartSums[i + 1] = boundingBoxes[i].size() + cellSizePartSums[i];
      }

      std::vector<unsigned int> neighborMatrix(cellSizePartSums.back() * neighborMatrixDim, UINT_MAX);

      for (size_t cell = 0; cell < cellsSize; ++cell) {
        for (size_t cluster = 0; cluster < (*_neighborCellIds)[cell].size(); ++cluster) {
          size_t i = 0;
          for (auto &neighbors : (*_neighborCellIds)[cell][cluster]) {
            neighborMatrix[(cellSizePartSums[cell] + cluster) * neighborMatrixDim + i] =
                cellSizePartSums[neighbors.first] + neighbors.second;
            ++i;
          }
          if (not useNewton3) {
            neighborMatrix[(cellSizePartSums[cell] + cluster) * neighborMatrixDim + i] =
                cellSizePartSums[cell] + cluster;
            ++i;
          }
        }
      }

#ifdef AUTOPAS_CUDA
      _neighborMatrix->copyHostToDevice(neighborMatrix.size(), neighborMatrix.data());
#endif
    }
  }

  void initTraversal() override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        for (size_t i = 0; i < (*this->_cells).size(); ++i) {
          _functor->SoALoader((*this->_cells)[i], (*this->_cells)[i]._particleSoABuffer);
        }

        return;
      }
      case DataLayoutOption::cuda: {
        size_t partSum = 0;
        for (size_t i = 0; i < (*this->_cells).size(); ++i) {
          _functor->SoALoader((*this->_cells)[i], _storageCell._particleSoABuffer, partSum);
          partSum += (*this->_cells)[i].numParticles();
        }
        _functor->deviceSoALoader(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
#ifdef AUTOPAS_CUDA
        utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
        return;
      }
    }
  }

  void endTraversal() override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < (*this->_cells).size(); ++i) {
          _functor->SoAExtractor((*this->_cells)[i], (*this->_cells)[i]._particleSoABuffer);
        }

        return;
      }
      case DataLayoutOption::cuda: {
        _functor->deviceSoAExtractor(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
#ifdef AUTOPAS_CUDA
        utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
        size_t partSum = 0;
        for (size_t i = 0; i < (*this->_cells).size(); ++i) {
          _functor->SoAExtractor((*this->_cells)[i], _storageCell._particleSoABuffer, partSum);
          partSum += (*this->_cells)[i].numParticles();
        }
        return;
      }
    }
  }

  void traverseParticlePairs() override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
        traverseCellPairsAoS(this->_cells);
        return;
      }
      case DataLayoutOption::soa: {
        traverseCellPairsSoA(this->_cells);
        return;
      }
      case DataLayoutOption::cuda: {
        traverseCellPairsGPU();
        return;
      }
    }
  }

 private:
  void traverseCellPairsAoS(std::vector<ParticleCell> *cells) {
    const auto clusterSize = *_clusterSize;

    // grid
    for (size_t i = 0; i < cells->size(); ++i) {
      // clusters
      for (size_t clusterId = 0; clusterId < (*_neighborCellIds)[i].size(); ++clusterId) {
        for (auto &neighbor : (*_neighborCellIds)[i][clusterId]) {
          // loop in cluster
          for (size_t ownPid = 0; ownPid < clusterSize; ++ownPid) {
            for (size_t otherPid = 0; otherPid < clusterSize; ++otherPid) {
              _functor->AoSFunctor((*cells)[i]._particles[clusterSize * clusterId + ownPid],
                                   (*cells)[neighbor.first]._particles[clusterSize * neighbor.second + otherPid],
                                   useNewton3);
            }
          }
        }
        // same cluster
        if (useNewton3) {
          for (size_t ownPid = 0; ownPid < clusterSize; ++ownPid) {
            for (size_t otherPid = ownPid + 1; otherPid < clusterSize; ++otherPid) {
              _functor->AoSFunctor((*cells)[i]._particles[clusterSize * clusterId + ownPid],
                                   (*cells)[i]._particles[clusterSize * clusterId + otherPid], useNewton3);
            }
          }
        } else {
          for (size_t ownPid = 0; ownPid < clusterSize; ++ownPid) {
            for (size_t otherPid = 0; otherPid < clusterSize; ++otherPid) {
              if (ownPid != otherPid) {
                _functor->AoSFunctor((*cells)[i]._particles[clusterSize * clusterId + ownPid],
                                     (*cells)[i]._particles[clusterSize * clusterId + otherPid], useNewton3);
              }
            }
          }
        }
      }
    }
  }

  void traverseCellPairsSoA(std::vector<ParticleCell> *cells) {
    auto clusterSize = *_clusterSize;
    // grid
    for (size_t i = 0; i < cells->size(); ++i) {
      // clusters
      for (size_t clusterId = 0; clusterId < (*_neighborCellIds)[i].size(); ++clusterId) {
        for (auto &neighbor : (*_neighborCellIds)[i][clusterId]) {
          const size_t c1start = clusterSize * clusterId;
          SoAView cluster1(&(*cells)[i]._particleSoABuffer, c1start, c1start + clusterSize);

          const size_t c2start = clusterSize * neighbor.second;
          SoAView cluster2(&(*cells)[neighbor.first]._particleSoABuffer, c2start, c2start + clusterSize);

          _functor->SoAFunctor(cluster1, cluster2, useNewton3);
        }
        // same cluster
        SoAView clusterSelf(&(*cells)[i]._particleSoABuffer, clusterId * clusterSize,
                            clusterId * clusterSize + clusterSize);
        _functor->SoAFunctor(clusterSelf, useNewton3);
      }
    }
  }

  void traverseCellPairsGPU() {
#ifdef AUTOPAS_CUDA
    if (!_functor->getCudaWrapper()) {
      _functor->CudaFunctor(_storageCell._particleSoABufferDevice, useNewton3);
      return;
    }

    auto cudaSoA = _functor->createFunctorCudaSoA(_storageCell._particleSoABufferDevice);

    if (useNewton3) {
      _functor->getCudaWrapper()->CellVerletTraversalN3Wrapper(
          cudaSoA.get(), _storageCell._particleSoABuffer.getNumParticles() / *_clusterSize, *_clusterSize,
          *_neighborMatrixDim, _neighborMatrix->get(), 0);
    } else {
      _functor->getCudaWrapper()->CellVerletTraversalNoN3Wrapper(
          cudaSoA.get(), _storageCell._particleSoABuffer.getNumParticles() / *_clusterSize, *_clusterSize,
          *_neighborMatrixDim, _neighborMatrix->get(), 0);
    }
    utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#else
    utils::ExceptionHandler::exception("VerletClusterCellsTraversal was compiled without Cuda support");
#endif
  }

  /**
   * Returns minimal distance between the two boxes
   * @param box1
   * @param box2
   * @return distance
   */
  inline double getMinDist(const std::array<double, 6> &box1, const std::array<double, 6> &box2) const {
    double sqrDist = 0;
    for (int i = 0; i < 3; ++i) {
      if (box2[i + 3] < box1[i]) {
        double d = box2[i + 3] - box1[i];
        sqrDist += d * d;
      } else if (box2[i] > box1[i + 3]) {
        double d = box2[i] - box1[i + 3];
        sqrDist += d * d;
      }
    }
    return sqrt(sqrDist);
  }

  /**
   * Pairwise functor used in this traversal
   */
  PairwiseFunctor *_functor;

  /**
   * SoA Storage cell containing SoAs and device Memory
   */
  ParticleCell _storageCell;

  // id of neighbor clusters of a clusters
  std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> *_neighborCellIds;

  size_t *_neighborMatrixDim;
  utils::CudaDeviceVector<unsigned int> *_neighborMatrix;

  unsigned int *_clusterSize;
};

}  // namespace autopas
