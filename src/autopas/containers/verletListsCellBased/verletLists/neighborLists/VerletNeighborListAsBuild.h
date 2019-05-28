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

template<class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C08TraversalColorChangeNotify : public C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3> {
 public:
  C08TraversalColorChangeNotify(const std::array<unsigned long, 3> &dims,
                                PairwiseFunctor *pairwiseFunctor,
                                TraversalColorChangeObserver *observer)
      : C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor),
        _observer(observer) {}
 protected:
  void notifyColorChange(unsigned long newColor) override {
    _observer->receiveColorChange(newColor);
  }
 private:
  TraversalColorChangeObserver *_observer;
};

template<class Particle>
class VerletNeighborListAsBuild : public VerletNeighborListInterface<Particle>, TraversalColorChangeObserver {
 private:
  /**
   * This functor can generate variable verlet lists using the typical pairwise
   * traversal.
   */
  class VarVerletListGeneratorFunctor : public autopas::Functor<Particle,
                                                                typename VerletListHelpers<Particle>::VerletListParticleCellType,
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
                    SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa2,
                    bool /*newton3*/) override {}

    void SoALoader(ParticleCell &cell,
                   SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa,
                   size_t offset = 0) override {}

    void SoAExtractor(ParticleCell &cell,
                      SoA<typename VerletListHelpers<Particle>::SoAArraysType> &soa,
                      size_t offset = 0) override {}

   private:
    VerletNeighborListAsBuild<Particle> &_list;
    double _cutoffskinsquared;
  }; // end functor


 public:
  using ThreadNeighborList = std::vector<std::pair<Particle *, Particle *>>;
  using ColorNeighborList = std::vector<ThreadNeighborList>;

  VerletNeighborListAsBuild()
      : _neighborList{}, currentColor(0) {}

  std::vector<TraversalOption> getAllTraversals() override {
    return std::vector<TraversalOption>{TraversalOption::varVerletTraversalAsBuild};
  };

  void buildNeighborList(LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                     typename VerletListHelpers<Particle>::SoAArraysType> &linkedCells,
                         bool useNewton3) override {
    //TODO: Maybe don't construct vectors for max number of threads
    unsigned int maxNumThreads = autopas_get_max_threads();
    for (int c = 0; c < 8; c++) {
      std::vector<ThreadNeighborList> &colorList = _neighborList[c];
      colorList.resize(maxNumThreads);
      for (unsigned int i = 0; i < maxNumThreads; i++) {
        //TODO: See if this call to clear takes a lot of performance because it is O(n)
        colorList[i].clear();
      }
    }

    VarVerletListGeneratorFunctor functor(*this, linkedCells.getCutoff());
    if (useNewton3) {
      auto traversal =
          C08TraversalColorChangeNotify<typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                        VarVerletListGeneratorFunctor,
                                        DataLayoutOption::aos,
                                        true>(
              linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &functor, this);
      linkedCells.iteratePairwise(&functor, &traversal);
    } else {
      auto traversal =
          C08TraversalColorChangeNotify<typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                        VarVerletListGeneratorFunctor,
                                        DataLayoutOption::aos,
                                        false>(
              linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &functor, this);
      linkedCells.iteratePairwise(&functor, &traversal);
    }

  }

  auto getInternalNeighborList() {
    return _neighborList;
  }

  void receiveColorChange(unsigned long newColor) override {
    currentColor = newColor;
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

  int currentColor;
};

} // namespace autopas
