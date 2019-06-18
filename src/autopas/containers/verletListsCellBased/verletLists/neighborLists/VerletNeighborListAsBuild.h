/**
 * @file VerletNeighborListAsBuild.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include "VerletNeighborListInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * Observer interface that specifies handling color changes during a CBasedTraversal.
 */
class TraversalColorChangeObserver {
 public:
  /**
   * Gets called when the color changes during the observed traversal.
   * @param newColor The new color that the traversal handles now.
   */
  virtual void receiveColorChange(unsigned long newColor) = 0;
};

/**
 * @copydoc C08Traversal
 *
 * It furthermore notifies an observer when color changes during the traversal happen.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C08TraversalColorChangeNotify : public C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3> {
 public:
  /**
   * Constructor of the traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param observer The observer to notify when a color change happens during the traversal.
   */
  C08TraversalColorChangeNotify(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                TraversalColorChangeObserver *observer)
      : C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor),
        _observer(observer) {}

 protected:
  /**
   * Notifies the observer.
   * @param newColor The new color to notify the observer about.
   */
  void notifyColorChange(unsigned long newColor) override {
    if (_observer) _observer->receiveColorChange(newColor);
  }

 private:
  /**
   * The observer of the traversal.
   */
  TraversalColorChangeObserver *_observer;
};

/**
 * This class implements a neighbor list that remembers which thread added which particle pair and at which color
 * during the build with C08 from LinkedCells.
 * @tparam Particle The particle type the class uses.
 */
template <class Particle>
class VerletNeighborListAsBuild : public VerletNeighborListInterface<Particle>, TraversalColorChangeObserver {
 private:
  /**
   * This functor can generate or check variable verlet lists using the typical pairwise
   * traversal.
   *
   * @tparam callCheckInstead If false, generate a neighbor list. If true, check the current for validity.
   */
  template <bool callCheckInstead = false>
  class VarVerletListPairGeneratorFunctor
      : public autopas::Functor<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                typename VerletListHelpers<Particle>::SoAArraysType> {
    typedef typename VerletListHelpers<Particle>::VerletListParticleCellType ParticleCell;

   public:
    /**
     * Constructor of the functor.
     * @param neighborList The neighbor list to fill.
     * @param cutoffskin The cutoff skin to use.
     */
    VarVerletListPairGeneratorFunctor(VerletNeighborListAsBuild &neighborList, double cutoffskin)
        : _list(neighborList), _cutoffskinsquared(cutoffskin * cutoffskin) {}

    bool isRelevantForTuning() override { return false; }

    /**
     * Adds the given pair to the neighbor list.
     * @param i The first particle of the pair.
     * @param j The second particle of the pair.
     * @param newton3 Not used!
     */
    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      double distsquare = ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffskinsquared) {
        if (callCheckInstead) {
          _list.checkPair(&i, &j);
        } else {
          _list.addPair(&i, &j);
        }
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
    /**
     * The neighbor list to fill.
     */
    VerletNeighborListAsBuild<Particle> &_list;
    /**
     * The squared cutoff skin to determine if a pair should be added to the list.
     */
    double _cutoffskinsquared;
  };  // end functor

  /**
   * Starts the generate functor. _baseLinkedCells has to be set before!
   * @tparam useNewton3 If the functor should use newton 3.
   * @tparam callCheckInstead If false, start it in generate mode, if true, in check validity mode.
   * @param cutoff The cutoff to use for the particle pairs in the functor.
   */
  template <bool useNewton3, bool callCheckInstead = false>
  void startFunctor(double cutoff) {
    VarVerletListPairGeneratorFunctor<callCheckInstead> functor(*this, cutoff);
    auto traversal = C08TraversalColorChangeNotify<typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                                   VarVerletListPairGeneratorFunctor<callCheckInstead>,
                                                   DataLayoutOption::aos, useNewton3>(
        _baseLinkedCells->getCellBlock().getCellsPerDimensionWithHalo(), &functor, this);
    _baseLinkedCells->iteratePairwise(&functor, &traversal);
  }

 public:
  /**
   * This type represents the neighbor list that each thread has for each color.
   */
  using ThreadNeighborList = std::unordered_map<Particle *, std::vector<Particle *>>;
  /**
   * This type represents the thread lists for all colors.
   */
  using ColorNeighborList = std::vector<ThreadNeighborList>;

  /**
   * This type represents the SoA neighbor list that each thread has for each color.
   */
  using SoAThreadNeighborList = std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>>;
  /**
   * This type represents the SoA thread lists for all colors.
   */
  using SoAColorNeighborList = std::vector<SoAThreadNeighborList>;

  /**
   * Constructor for the VerletNeighborListAsBuild. Does only default initialization.
   */
  VerletNeighborListAsBuild() : _neighborList{}, _soaListIsValid(false), _currentColor(0) {}

  ContainerOption getContainerType() const override { return ContainerOption::varVerletListsAsBuild; }

  /**
   * @copydoc VerletNeighborListInterface::buildNeighborList()
   *
   * It executes C08 on the passed LinkedCells container and saves the resulting pairs in the neighbor list, remembering
   * the thread and current color for each pair.
   */
  void buildNeighborList(LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                     typename VerletListHelpers<Particle>::SoAArraysType> &linkedCells,
                         bool useNewton3) override {
    _soaListIsValid = false;
    _baseLinkedCells = &linkedCells;

    unsigned int maxNumThreads = autopas_get_max_threads();
    for (int c = 0; c < 8; c++) {
      std::vector<ThreadNeighborList> &colorList = _neighborList[c];
      colorList.resize(maxNumThreads);
      for (unsigned int i = 0; i < maxNumThreads; i++) {
        colorList[i].clear();
      }
    }

    if (useNewton3) {
      startFunctor<true>(linkedCells.getCutoff());
    } else {
      startFunctor<false>(linkedCells.getCutoff());
    }
  }

  bool checkNeighborListValidity(bool useNewton3, double cutoff) override {
    _allPairsPresent = true;
    if (_baseLinkedCells == nullptr) return false;

    constexpr bool callCheck = true;
    if (useNewton3) {
      startFunctor<true, callCheck>(cutoff);
    } else {
      startFunctor<false, callCheck>(cutoff);
    }

    return _allPairsPresent;
  }

  /**
   * Returns the internal AoS neighbor list. Should be used by traversals.
   *
   * The internal neighbor list is a vector of the neighbor lists of each color. Each of those neighbor lists is a
   * vector that contains one neighbor list for each thread. Each of those neighbor lists is a map from each particle to
   * a vector containing its neighbors.
   *
   * @return the internal AoS neighbor list.
   */
  const auto &getInternalNeighborList() { return _neighborList; }

  /**
   * Returns the internal SoA neighbor list. Should be used by traversals.
   *
   * The internal SoA neighbor list is a vector of the SoA neighbor lists of each color. Each of those SoA neighbor
   * lists is a vector that contains one SoA neighbor list for each thread. Each of those SoA neighbor lists is a vector
   * of vectors where the i-th vector contains the indices of all neighbors of particle i in the SoA.
   *
   * @return the internal SoA neighbor list.
   */
  const auto &getInternalSoANeighborList() { return _soaNeighborList; }

  void receiveColorChange(unsigned long newColor) override { _currentColor = newColor; }

  /**
   * @see getInternalSoANeighborList()
   */
  void generateSoAFromAoS() override {
    // Generate a map from pointer to particle index in the SoA. This works, because during loadSoA"()" the particles
    // are loaded in the same order.
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
          for (const auto &second : pair.second) {
            size_t indexSecond = _aos2soaMap[second];
            currentThreadList[indexFirst].push_back(indexSecond);
          }
        }
      }
    }

    _soaListIsValid = true;
  }

  /**
   * Loads the particle information in the SoA and returns a pointer to the filled SoA.
   * @tparam TFunctor The type of the functor to use for loading the particles.
   * @param f The functor to use for loading the particles.
   * @return A pointer to the SoA filled. Ownership is *not* passed.
   */
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
  /**
   * Extracts the particle information out of the SoA returned by loadSoA() before.
   * @tparam TFunctor The type of the functor to use for extracting the particles.
   * @param f The functor to use for extracting the particles.
   */
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
  /**
   * Called from VarVerletListGeneratorFunctor
   */
  void addPair(Particle *first, Particle *second) {
    int currentThreadIndex = autopas_get_thread_num();
    _neighborList[_currentColor][currentThreadIndex][first].push_back(second);
  }

  /**
   * Called from VarVerletListGeneratorFunctor
   */
  void checkPair(Particle *first, Particle *second) {
    int currentThreadIndex = autopas_get_thread_num();

    // Check all neighbor lists for the pair, but the one that the pair would be in if it was not moved first.
    auto &oldThreadNeighborList = _neighborList[_currentColor][currentThreadIndex];
    if (isPairInList(oldThreadNeighborList, first, second)) {
      for (int color = 0; color < 8; color++) {
        for (unsigned int thread = 0; thread < _neighborList[color].size(); thread++) {
          if (not isPairInList(_neighborList[_currentColor][currentThreadIndex], first, second)) {
            // this is thread safe, as _allPairsPresent is atomic
            _allPairsPresent = false;
            return;
          }
        }
      }
    } else {
      // this is thread safe, as _allPairsPresent is atomic
      _allPairsPresent = false;
    }
  }

  /**
   * Helper method for checkPair()
   * @return True, if the pair is not present, false otherwise.
   */
  bool isPairInList(ThreadNeighborList &currentNeighborList, Particle *first, Particle *second) {
    auto found = std::find(currentNeighborList[first].begin(), currentNeighborList[first].end(), second);
    if (found == currentNeighborList[first].end()) {
      return false;
    }
    return true;
  }

 private:
  /**
   * The internal AoS neighbor list. For format, see getInternalNeighborList().
   */
  std::array<ColorNeighborList, 8> _neighborList;

  /**
   * The LinkedCells object this neighbor list should use to build.
   */
  LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
              typename VerletListHelpers<Particle>::SoAArraysType> *_baseLinkedCells;

  /**
   * The internal SoA neighbor list. For format, see getInternalSoANeighborList().
   */
  std::array<SoAColorNeighborList, 8> _soaNeighborList;

  /**
   * The SoA used.
   */
  SoA<typename Particle::SoAArraysType> _soa;

  /**
   * If the SoA is valid, see isSoAListValid().
   */
  bool _soaListIsValid;

  /**
   * The current color in the traversal during the build of the neighbor list.
   */
  int _currentColor;

  /**
   * Used in checkNeighborListValidity(). Set to false in the pair generating functor.
   */
  std::atomic<bool> _allPairsPresent;
};

}  // namespace autopas
