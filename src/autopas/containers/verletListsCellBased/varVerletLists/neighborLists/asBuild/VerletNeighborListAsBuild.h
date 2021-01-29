/**
 * @file VerletNeighborListAsBuild.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

#include "AsBuildPairGeneratorFunctor.h"
#include "C08TraversalColorChangeNotify.h"
#include "autopas/containers/verletListsCellBased/varVerletLists/neighborLists/VerletNeighborListInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class implements a neighbor list that remembers which thread added which particle pair and at which color
 * during the build with C08 from LinkedCells.
 * @tparam Particle The particle type the class uses.
 */
template <class Particle>
class VerletNeighborListAsBuild : public VerletNeighborListInterface<Particle>, ColorChangeObserver {
  /**
   * Adds the generator functor for validation checks as friend so it can call checkPair().
   *  @param true mark that it is for validation checks.
   */
  friend class internal::AsBuildPairGeneratorFunctor<Particle, true>;
  /**
   * Adds the generator functor for adding pairs as friend so it can call addPair().
   *  @param false test mark that it is for adding pairs.
   */
  friend class internal::AsBuildPairGeneratorFunctor<Particle, false>;

 private:
  /**
   * Applies the generate functor. _baseLinkedCells has to be set before!
   * @tparam useNewton3 If the functor should use newton 3.
   * @tparam validationMode If false, start it in generate mode, if true, in check validity mode.
   * @param cutoff The cutoff to use for the particle pairs in the functor.
   */
  template <bool useNewton3, bool validationMode = false>
  void applyGeneratorFunctor(double cutoff) {
    internal::AsBuildPairGeneratorFunctor<Particle, validationMode> generatorFunctor(*this, cutoff);
    // Use SoA traversal for generation and AoS traversal for validation check.
    constexpr auto dataLayout = validationMode ? DataLayoutOption::aos : DataLayoutOption::soa;
    auto traversal = C08TraversalColorChangeNotify<FullParticleCell<Particle>,
                                                   internal::AsBuildPairGeneratorFunctor<Particle, validationMode>,
                                                   dataLayout, useNewton3>(
        _baseLinkedCells->getCellBlock().getCellsPerDimensionWithHalo(), &generatorFunctor,
        _baseLinkedCells->getInteractionLength(), _baseLinkedCells->getCellBlock().getCellLength(), this);
    _baseLinkedCells->iteratePairwise(&traversal);
  }

 public:
  /**
   * This type represents the neighbor list that each thread has for each color.
   */
  using AoSThreadNeighborList = std::unordered_map<Particle *, std::vector<Particle *>>;
  /**
   * This type represents the thread lists for all colors.
   */
  using AoSColorNeighborList = std::vector<AoSThreadNeighborList>;

  /**
   * This type represents the SoA neighbor list that each thread has for each color.
   */
  using SoAThreadNeighborList = std::vector<std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>>;
  /**
   * This type represents the SoA thread lists for all colors.
   */
  using SoAColorNeighborList = std::vector<SoAThreadNeighborList>;

  /**
   * Constructor for the VerletNeighborListAsBuild. Does only default initialization.
   */
  VerletNeighborListAsBuild() : _aosNeighborList{}, _soaListIsValid(false) {}

  /**
   * @copydoc VerletNeighborListInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::varVerletListsAsBuild; }

  /**
   * @copydoc VerletNeighborListInterface::buildAoSNeighborList()
   *
   * It executes C08 on the passed LinkedCells container and saves the resulting pairs in the neighbor list, remembering
   * the thread and current color for each pair.
   */
  void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3) override {
    _soaListIsValid = false;
    _baseLinkedCells = &linkedCells;

    auto maxNumThreads = autopas_get_max_threads();
    for (int color = 0; color < _numColors; color++) {
      _aosNeighborList[color].resize(maxNumThreads);
      for (auto &colorList : _aosNeighborList[color]) {
        colorList.clear();
      }
    }

    if (useNewton3) {
      applyGeneratorFunctor<true>(linkedCells.getInteractionLength());
    } else {
      applyGeneratorFunctor<false>(linkedCells.getInteractionLength());
    }
  }

  /**
   * @copydoc VerletNeighborListInterface::checkNeighborListValidity()
   */
  bool checkNeighborListValidity(bool useNewton3, double cutoff) override {
    _allPairsPresent = true;
    if (_baseLinkedCells == nullptr) return false;

    constexpr bool callCheck = true;
    if (useNewton3) {
      applyGeneratorFunctor<true, callCheck>(cutoff);
    } else {
      applyGeneratorFunctor<false, callCheck>(cutoff);
    }

    return _allPairsPresent;
  }

  /**
   * Returns the internal AoS neighbor list. Should be used by traversals.
   *
   * The internal neighbor list structure is an array of vectors for each color. Each of those vectors contains a
   * neighbor list for each thread. Each of those neighbor lists is a map from particle pointers to a vector containing
   * its neighbor pointers.
   *
   * Or in short:
   * _aosNeighborList
   * = std::array<AoSColorNeighborList, _numColors>
   * = std::array<std::vector<AoSThreadNeighborList>>, _numColors>
   * = std::array<std::vector<std::unordered_map<Particle *, std::vector<Particle *>>>, _numColors>
   *
   * @return the internal AoS neighbor list.
   */
  const auto &getAoSNeighborList() { return _aosNeighborList; }

  /**
   * Returns the internal SoA neighbor list. Should be used by traversals.
   *
   * The internal SoA neighbor list structure is an array of vectors for each color. Each of those vectors
   * contains one SoA neighbor list per thread. Each of those SoA neighbor lists is a vector of pairs mimicing a map.
   * Each pair contains an index in the SoA and a vector of the indices of all its neighbors in the SoA.
   *
   * Or in short:
   * _soaNeighborList
   * = std::array<SoAColorNeighborList, _numColors>
   * = std::array<std::vector<SoAThreadNeighborList>, _numColors>
   * = std::array<std::vector<std::pair<size_t, std::vector<size_t, autopas::AlignedAllocator<size_t>>>>, _numColors>
   *
   * @return the internal SoA neighbor list.
   */
  const auto &getSoANeighborList() { return _soaNeighborList; }

  /**
   * @copydoc ColorChangeObserver::receiveColorChange()
   */
  void receiveColorChange(uint64_t newColor) override { _currentColor = newColor; }

  /**
   * @see getSoANeighborList()
   */
  void generateSoAFromAoS() override {
    // Generate a map from pointer to particle index in the SoA. This works, because during loadSoA"()" the particles
    // are loaded in the same order.
    std::unordered_map<Particle *, size_t> _aos2soaMap;
    _aos2soaMap.reserve(_baseLinkedCells->getNumParticles());
    size_t i = 0;
    // needs to iterate also over dummies!
    for (auto iter = _baseLinkedCells->begin(IteratorBehavior::haloOwnedAndDummy); iter.isValid(); ++iter, ++i) {
      _aos2soaMap[&(*iter)] = i;
    }

    for (int color = 0; color < _numColors; color++) {
      unsigned int numThreads = _aosNeighborList[color].size();
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
        for (const auto &pair : _aosNeighborList[color][thread]) {
          size_t indexFirst = _aos2soaMap[pair.first];
          std::vector<size_t, AlignedAllocator<size_t>> neighbors;
          neighbors.reserve(pair.second.size());
          for (const auto &second : pair.second) {
            size_t indexSecond = _aos2soaMap[second];
            neighbors.push_back(indexSecond);
          }
          currentThreadList.push_back({indexFirst, std::move(neighbors)});
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

  /**
   * @copydoc VerletNeighborListInterface::isSoAListValid()
   */
  bool isSoAListValid() const override { return _soaListIsValid; }

  /**
   * @copydoc VerletNeighborListInterface::getNumberOfNeighborPairs()
   */
  long getNumberOfNeighborPairs() const override {
    long numPairs = 0;
    for (const auto &colorList : _aosNeighborList) {
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
    _aosNeighborList[_currentColor][currentThreadIndex][first].push_back(second);
  }

  /**
   * Called from VarVerletListGeneratorFunctor
   */
  void checkPair(Particle *first, Particle *second) {
    int currentThreadIndex = autopas_get_thread_num();

    // Check all neighbor lists for the pair, but the one that the pair would be in if it was not moved first.
    auto &oldThreadNeighborList = _aosNeighborList[_currentColor][currentThreadIndex];
    if (isPairInList(oldThreadNeighborList, first, second)) {
      for (int color = 0; color < _numColors; color++) {
        for (unsigned int thread = 0; thread < _aosNeighborList[color].size(); thread++) {
          if (not isPairInList(_aosNeighborList[_currentColor][currentThreadIndex], first, second)) {
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
   * @return True, if the pair is present, false otherwise.
   */
  bool isPairInList(AoSThreadNeighborList &currentNeighborList, Particle *first, Particle *second) {
    auto iteratorFound = std::find(currentNeighborList[first].begin(), currentNeighborList[first].end(), second);

    return iteratorFound != currentNeighborList[first].end();
  }

 private:
  /**
   * Number of colors used for the domain coloring during parallelization.
   */
  constexpr static size_t _numColors = 8;

  /**
   * The internal AoS neighbor list. For format, see getAoSNeighborList().
   */
  std::array<AoSColorNeighborList, _numColors> _aosNeighborList;

  /**
   * The LinkedCells object this neighbor list should use to build.
   */
  LinkedCells<Particle> *_baseLinkedCells;

  /**
   * The internal SoA neighbor list. For format, see getSoANeighborList().
   */
  std::array<SoAColorNeighborList, _numColors> _soaNeighborList;

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
  int _currentColor{0};

  /**
   * Used in checkNeighborListValidity(). Set to false in the pair generating functor.
   */
  std::atomic<bool> _allPairsPresent;
};

}  // namespace autopas
