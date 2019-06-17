/**
 * @file VarVerletListAsBuild.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include "VerletNeighborListInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

class TraversalColorChangeObserver {
 public:
  virtual void receiveColorChange(unsigned long newColor) = 0;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C08TraversalColorChangeNotify : public C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3> {
 public:
  C08TraversalColorChangeNotify(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                TraversalColorChangeObserver *observer)
      : C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor),
        _observer(observer) {}

 protected:
  void notifyColorChange(unsigned long newColor) override { _observer->receiveColorChange(newColor); }

 private:
  TraversalColorChangeObserver *_observer;
};

template <class Particle>
class VerletNeighborListAsBuild : public VerletNeighborListInterface<Particle>, TraversalColorChangeObserver {
 private:
  /**
   * This functor can generate variable verlet lists using the typical pairwise
   * traversal.
   */
  class VarVerletListGeneratorFunctor
      : public autopas::Functor<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                typename VerletListHelpers<Particle>::SoAArraysType> {
    typedef typename VerletListHelpers<Particle>::VerletListParticleCellType ParticleCell;

   public:
    /**
     * Constructor
     * @param neighborList
     * @param cutoffskin
     */
    VarVerletListGeneratorFunctor(VerletNeighborListAsBuild &neighborList, double cutoffskin)
        : _list(neighborList), _cutoffskinsquared(cutoffskin * cutoffskin) {}

    bool isRelevantForTuning() override { return false; }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      double distsquare = ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffskinsquared) {
        _list.addPair(&i, &j);
      }
    }

    void SoAFunctor(SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa, bool newton3) override {}

    void SoAFunctor(SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa1,
                    SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa2, bool /*newton3*/) override {}

    void SoALoader(ParticleCell &cell, SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa,
                   size_t offset = 0) override {}

    void SoAExtractor(ParticleCell &cell, SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa,
                      size_t offset = 0) override {}

   private:
    VerletNeighborListAsBuild<Particle> &_list;
    double _cutoffskinsquared;
  };  // end functor

 public:
  using ThreadNeighborList = std::vector<std::pair<Particle *, Particle *>>;
  using ColorNeighborList = std::vector<ThreadNeighborList>;

  using SoAThreadNeighborList = std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>>;
  using SoAColorNeighborList = std::vector<SoAThreadNeighborList>;

  VerletNeighborListAsBuild() : _neighborList{}, _soaListIsValid(false), currentColor(0) {}

  std::vector<TraversalOption> getAllTraversals() const override {
    return std::vector<TraversalOption>{TraversalOption::varVerletTraversalAsBuild};
  };

  ContainerOption getContainerType() const override { return ContainerOption::varVerletListsAsBuild; }

  void buildNeighborList(LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                     typename VerletListHelpers<Particle>::SoAArraysType> &linkedCells,
                         bool useNewton3) override {
    _soaListIsValid = false;
    _baseLinkedCells = &linkedCells;

    // TODO: Maybe don't construct vectors for max number of threads
    unsigned int maxNumThreads = autopas_get_max_threads();
    for (int c = 0; c < 8; c++) {
      std::vector<ThreadNeighborList> &colorList = _neighborList[c];
      colorList.resize(maxNumThreads);
      for (unsigned int i = 0; i < maxNumThreads; i++) {
        // TODO: See if this call to clear takes a lot of performance because it is O(n)
        colorList[i].clear();
      }
    }

    VarVerletListGeneratorFunctor functor(*this, linkedCells.getCutoff());
    if (useNewton3) {
      auto traversal = C08TraversalColorChangeNotify<typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                                     VarVerletListGeneratorFunctor, DataLayoutOption::aos, true>(
          linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &functor, this);
      linkedCells.iteratePairwise(&functor, &traversal);
    } else {
      auto traversal = C08TraversalColorChangeNotify<typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                                     VarVerletListGeneratorFunctor, DataLayoutOption::aos, false>(
          linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &functor, this);
      linkedCells.iteratePairwise(&functor, &traversal);
    }
  }

  const auto &getInternalNeighborList() { return _neighborList; }

  const auto &getInternalSoANeighborList() { return _soaNeighborList; }

  void receiveColorChange(unsigned long newColor) override { currentColor = newColor; }

  void generateSoAFromAoS() override {
    std::unordered_map<Particle *, size_t> _aos2soaMap;
    _aos2soaMap.reserve(_baseLinkedCells->getNumParticles());
    size_t i = 0;
    for (auto iter = _baseLinkedCells->begin(); iter.isValid(); ++iter, ++i) {
      _aos2soaMap[&(*iter)] = i;
    }

    constexpr int numColors = 8;
    for (int color = 0; color < numColors; color++) {
      unsigned int numThreads = _neighborList[color].size();
      _soaNeighborList[color].resize(numThreads);
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel num_threads(numThreads)
#endif
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(static)
#endif
      for (unsigned int thread = 0; thread < numThreads; thread++) {
        auto &currentThreadList = _soaNeighborList[color][thread];
        currentThreadList.clear();
        currentThreadList.resize(_aos2soaMap.size());
        for (const auto &pair : _neighborList[color][thread]) {
          size_t indexFirst = _aos2soaMap[pair.first];
          size_t indexSecond = _aos2soaMap[pair.second];
          currentThreadList[indexFirst].push_back(indexSecond);
        }
      }
    }

    _soaListIsValid = true;
  }

  template <class TFunctor>
  auto *loadSoA(TFunctor *f) {
    _soa.clear();
    size_t offset = 0;
    for (auto &cell : _baseLinkedCells->getCells()) {
      f->SoALoader(cell, _soa, offset);
      offset += cell.numParticles();
    }

    return &_soa;
  }

  template <class TFunctor>
  void extractSoA(TFunctor *f) {
    size_t offset = 0;
    for (auto &cell : _baseLinkedCells->getCells()) {
      f->SoAExtractor(cell, _soa, offset);
      offset += cell.numParticles();
    }
  }

  bool isSoAListValid() const override { return _soaListIsValid; }

  long getNumberOfNeighborPairs() const override {
    long numPairs = 0;
    for (const auto &colorList : _neighborList) {
      for (const auto &threadList : colorList) {
        numPairs += threadList.size();
      }
    }
    return numPairs;
  }

 private:
  /*
   * Called from VarVerletListGeneratorFunctor
   */
  void addPair(Particle *first, Particle *second) {
    int currentThreadIndex = autopas_get_thread_num();
    _neighborList[currentColor][currentThreadIndex].push_back({first, second});
  }

 private:
  std::array<ColorNeighborList, 8> _neighborList;

  LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
              typename VerletListHelpers<Particle>::SoAArraysType> *_baseLinkedCells;

  std::array<SoAColorNeighborList, 8> _soaNeighborList;
  SoA<typename Particle::SoAArraysType> _soa;
  bool _soaListIsValid;

  int currentColor;
};

}  // namespace autopas
