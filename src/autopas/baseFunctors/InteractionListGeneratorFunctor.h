/**
 * @file InteractionListGeneratorFunctor.h
 * @author seckler
 * @date 27.04.18
 */

#pragma once

#include <atomic>

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {

/**
 * This functor generates lists of particles within interactionLength of each other, which could be used within
 * user simulators to replace their contact detection passes, although they should consider the warning below.
 *
 * @warning Performing contact detection with this functor and then applying forces externally to AutoPas is not
 * recommended, as it does not take advantage of AutoPas's full algorithm library and will also be at best inefficient
 * or at worse unsupported by any future GPU extension of AutoPas. We provide this functor primarily for codes which
 * perform a separate contact detection and force calculation passes, to easily experience some benefit of AutoPas. For
 * the full capabilities of AutoPas, we recommend writing a functor class that directly applies relevant interactions.
 * See applicationLibrary for examples.
 *
 * @details After applying AutoPas's computeInteractions function with this functor, a
 * std::unordered_map<Particle_T *, std::vector<Particle_T *>> is filled, mapping from each particle pointer to a vector
 * of pointers to all particles within interactionLength of the first.
 *
 * @note This functor is also used internally for some Verlet List generation.
 *
 * @tparam Particle_T The type of Particle class used.
 * @tparam isInternal Should be set true if used internally within AutoPas. Makes the functor irrelevant for tuning.
 */
template <class Particle_T, bool isInternal = false>
class InteractionListGeneratorFunctor
    : public PairwiseFunctor<Particle_T, InteractionListGeneratorFunctor<Particle_T>> {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = Particle_T::SoAArraysType;

  /**
   * Neighbor list AoS style.
   */
  using NeighborListAoSType = std::unordered_map<Particle_T *, std::vector<Particle_T *>>;

  /**
   * Constructor
   * @param neighborListsAoS
   * @param interactionLength The distance between particles within which particle pairs get added to the neighbor
   * lists.
   * @param gatherNewton3Lists If false, for a particle pair i, j in contact, **both** particle j will be in particle
   * i's list **and** particle i will be in particle j's list. If true, for a particle pair i, j in contact, **either**
   * particle j will be in particle i's list **or** particle i will be in particle j's list. The latter can be used more
   * easily to reduce calculations by applying forces to both particles in the pair, but care should be taken to avoid
   * race conditions. If true, we make no guarantee which of particle i or j will be in the other's list. If true, this
   * functor will only allow functor calls with newton3 enabled.
   */
  InteractionListGeneratorFunctor(NeighborListAoSType &neighborListsAoS, double interactionLength,
                                  bool gatherNewton3Lists)
      : PairwiseFunctor<Particle_T, InteractionListGeneratorFunctor>(interactionLength),
        _neighborListsAoS(neighborListsAoS),
        _interactionLengthSquared(interactionLength * interactionLength),
        _gatherNewton3Lists(gatherNewton3Lists) {}

  std::string getName() override { return "InteractionListGeneratorFunctor"; }

  /**
   * Whether particle is relevant for tuning. True, unless functor used internally.
   * @return
   */
  bool isRelevantForTuning() override { return not isInternal; }

  /**
   * Whether InteractionListGeneratorFunctor allows non-newton3. Is always allowed.
   * @return
   */
  bool allowsNewton3() override { return true; }

  /**
   * Whether InteractionListGeneratorFunctor allows non-newton3. This is not allowed if not gathering N3 lists. (Not
   * a fundamental issue but messy implementation and not really needed).
   * @return
   */
  bool allowsNonNewton3() override { return not _gatherNewton3Lists; }

  /**
   * AoSFunctor for interaction list generation.
   * @param i
   * @param j
   * @param newton3 Whether AutoPas is using N3 or not for list generation. This is regardless of whether the user
   * requests N3 lists or not. If this and gatherNewton3Lists match, the handling is trivial. If newton3=true and
   * gatherNewton3Lists=false, we just add each particle to each other's list. We do not allow newton3=false with
   * gatherNewton3Lists=true - there are no fundamental issues with this but messy to implement and probably not too
   * needed.
   */
  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
    using namespace autopas::utils::ArrayMath::literals;

    [[unlikely]] if (_gatherNewton3Lists and not newton3) {
      utils::ExceptionHandler::exception(
          "InteractionListGeneratorFunctor should not be used with newton3=false and gatherNewton3Lists=true.");
    }

    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto dist = i.getR() - j.getR();

    double distsquare = utils::ArrayMath::dot(dist, dist);
    if (distsquare < _interactionLengthSquared) {
      // Assuming this functor is used like any other functor, this is thread safe: _neighborListsAoS is an
      // unordered_map, meaning we can push_back to particle i's list with the same thread safety as for writing to its
      // force buffer in e.g. a LJ functor.

      // This is only thread-safe if all keys are inserted prior to applying this functor, as a rehash might lead to
      // dangling references, but at() protects against this: if the key exists, we can push_back safely; if a key
      // doesn't, an exception is thrown anyway.

      // - If newton3=false & gatherNewton3Lists=false, we only need to add the i->j interaction to i's list as the j->i
      // interaction is handled in another call.
      // - If newton3=true & gatherNewton3Lists=false, we need to add the i->j interaction to j's list and the j->i
      // interaction to j's list.
      // - If newton3=true & gatherNewton3Lists=true, we only need to add the i->j interaction to i's list as we don't
      // want the j->i interaction in the lists. (Could also be the other way around)

      _neighborListsAoS.at(&i).push_back(&j);
      if (newton3 and not _gatherNewton3Lists) {
        _neighborListsAoS.at(&j).push_back(&i);
      }
    }
  }

  /**
   * SoAFunctor for verlet list generation. (single cell version)
   * @param soa the soa
   * @param newton3 Whether AutoPas is using N3 or not for list generation. For this function, this is ignored and
   * N3 is always used. This is regardless of whether the user requests N3 lists or not. If this and gatherNewton3Lists
   * match, the handling is trivial. If newton3=true and gatherNewton3Lists=false, we just add each particle to each
   * other's list.
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool /*newton3*/) override {
    if (soa.size() == 0) return;

    auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

    size_t numPart = soa.size();
    for (unsigned int i = 0; i < numPart; ++i) {
      for (unsigned int j = i + 1; j < numPart; ++j) {
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _interactionLengthSquared) {
          // This SoAFunctorSingle implementation always used newton3 in practice, regardless of the option.
          _neighborListsAoS.at(ptrptr[i]).push_back(ptrptr[j]);
          if (not _gatherNewton3Lists) {
            _neighborListsAoS.at(ptrptr[j]).push_back(ptrptr[i]);
          }
        }
      }
    }
  }

  /**
   * SoAFunctor for the verlet list generation. (two cell version)
   * @param soa1 soa of first cell
   * @param soa2 soa of second cell
   * @param newton3 Whether AutoPas is using N3 or not for list generation. This is regardless of whether the user
   * requests N3 lists or not. If this and gatherNewton3Lists match, the handling is trivial. If newton3=true and
   * gatherNewton3Lists=false, we just add each particle to each other's list. We do not allow newton3=false with
   * gatherNewton3Lists=true - there are no fundamental issues with this but messy to implement and probably not too
   * needed.
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3) override {
    [[unlikely]] if (_gatherNewton3Lists and not newton3) {
      utils::ExceptionHandler::exception(
          "InteractionListGeneratorFunctor should not be used with newton3=false and gatherNewton3Lists=true.");
    }

    if (soa1.size() == 0 or soa2.size() == 0) return;

    auto **const __restrict ptr1ptr = soa1.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();

    auto **const __restrict ptr2ptr = soa2.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();

    size_t numPart1 = soa1.size();
    for (unsigned int i = 0; i < numPart1; ++i) {
      size_t numPart2 = soa2.size();

      for (unsigned int j = 0; j < numPart2; ++j) {
        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _interactionLengthSquared) {
          _neighborListsAoS.at(ptr1ptr[i]).push_back(ptr2ptr[j]);
          if (newton3 and not _gatherNewton3Lists) {
            _neighborListsAoS.at(ptr2ptr[j]).push_back(ptr1ptr[i]);
          }
        }
      }
    }
  }

  /**
   * SoAFunctorVerlet for interaction list generation.
   *
   * Generates interaction list entries for the particle at index indexFirst in soa and every potential neighbor
   * given by the Verlet list.
   *
   * @param soa the soa
   * @param indexFirst index of the particle whose neighbors are given in verletList
   * @param verletList indices of the potential neighbors of indexFirst
   * @param newton3 Whether AutoPas is using N3 or not for list generation. This is regardless of whether the user
   * requests N3 lists or not. If this and gatherNewton3Lists match, the handling is trivial. If newton3=true and
   * gatherNewton3Lists=false, we just add each particle to each other's list. We do not allow newton3=false with
   * gatherNewton3Lists=true - there are no fundamental issues with this but messy to implement and probably not too
   * needed.
   */
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, AlignedAllocator<size_t>> &verletList, bool newton3) override {
    [[unlikely]] if (_gatherNewton3Lists and not newton3) {
      utils::ExceptionHandler::exception(
          "InteractionListGeneratorFunctor should not be used with newton3=false and gatherNewton3Lists=true.");
    }

    if (soa.size() == 0 or verletList.empty()) return;

    auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto &firstList = _neighborListsAoS.at(ptrptr[indexFirst]);
    const double xFirst = xptr[indexFirst];
    const double yFirst = yptr[indexFirst];
    const double zFirst = zptr[indexFirst];

    for (const size_t j : verletList) {
      const double drx = xFirst - xptr[j];
      const double dry = yFirst - yptr[j];
      const double drz = zFirst - zptr[j];

      const double drx2 = drx * drx;
      const double dry2 = dry * dry;
      const double drz2 = drz * drz;

      const double dr2 = drx2 + dry2 + drz2;

      if (dr2 < _interactionLengthSquared) {
        firstList.push_back(ptrptr[j]);
        if (newton3 and not _gatherNewton3Lists) {
          _neighborListsAoS.at(ptrptr[j]).push_back(ptrptr[indexFirst]);
        }
      }
    }
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static std::array<typename Particle_T::AttributeNames, 4> getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 4>{
        Particle_T::AttributeNames::ptr, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
        Particle_T::AttributeNames::posZ};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static std::array<typename Particle_T::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 0>{/*Nothing*/};
  }

 private:
  NeighborListAoSType &_neighborListsAoS;
  double _interactionLengthSquared;

  /**
   * If true, each particle pair will appear in only one particle's list. Else, each particle pair will appear in each
   * particle's list.
   */
  bool _gatherNewton3Lists{false};
};

}  // namespace autopas
