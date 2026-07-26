/**
 * @file NeighborIdentificationFunctor.h
 * @author S. Newcome
 * @date 22.07.2026
 */

#pragma once

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"

namespace autopas {

/**
 * This functor generates lists of particles within interactionLength of each other, which could be used within
 * user simulators to replace their neighbor identification/contact detection/cutoff check passes, although they should
 * consider the warning below.
 *
 * In some fields this can be called "cutoff checking" or "broad-phase contact detection" but we will refer to this
 * in this file as neighbor identification.
 *
 * @warning To take fully advantage of AutoPas, we recommend, instead of using this functor and apply forces externally
 * to AutoPas, using a functor which  integrates neighbor identification and force calculation in one call over using
 * this functor (see @ref LJFunctor.h as an example). We provide this functor primarily for codes which already
 * perform separate neighbor identification and force calculation passes, to easily experience some benefit of AutoPas
 * with minimal code refactors.
 *
 * @details After applying AutoPas's computeInteractions function with this functor, a
 * std::unordered_map<Particle_T *, std::vector<Particle_T *>> is filled, mapping from each particle pointer to a vector
 * of pointers to all particles within interactionLength of the first.
 *
 * @tparam Particle_T The type of Particle class used.
 */
template <class Particle_T>
class NeighborIdentificationFunctor : public PairwiseFunctor<Particle_T, NeighborIdentificationFunctor<Particle_T>> {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Neighbor list AoS style.
   */
  using NeighborListAoSType = std::unordered_map<Particle_T *, std::vector<Particle_T *>>;

  /**
   * Constructor
   * @param neighborListsAoS Reference to the neighbor list map. When used with AutoPas::computeInteractions, the map
   * will be overridden.
   * @param interactionLength The distance between particles within which particle pairs get added to the neighbor
   * lists.
   * @param gatherNewton3Lists If false, for a particle pair i, j that are neighbors, **both** particle j will be in
   * particle i's list **and** particle i will be in particle j's list. If true, for a particle pair i, j that are
   * neighbors, **either** particle j will be in particle i's list **or** particle i will be in particle j's list. The
   * latter can be used more easily to reduce calculations by applying forces to both particles in the pair, but care
   * should be taken to avoid race conditions. If true, we make no guarantee which of particle i or j will be in the
   * other's list. If true, this functor will only allow functor calls with newton3 enabled.
   */
  NeighborIdentificationFunctor(NeighborListAoSType &neighborListsAoS, double interactionLength,
                                bool gatherNewton3Lists)
      : PairwiseFunctor<Particle_T, NeighborIdentificationFunctor>(interactionLength),
        _neighborListsAoS(neighborListsAoS),
        _interactionLengthSquared(interactionLength * interactionLength),
        _gatherNewton3Lists(gatherNewton3Lists) {}

  std::string getName() override { return "NeighborIdentificationFunctor"; }

  /**
   * Initializes the neighbor list map for the given range of particles: clears any previous contents and inserts an
   * empty neighbor list for every particle. This must happen before the functor is applied, and it must also happen
   * after any particle rearrangement in memory. AutoPas's computeInteractions() calls
   * this automatically when this functor (or a child of it) is used, so the handling of the list is safe when used
   * with computeInteractions.
   *
   * @tparam ParticleIterator_T Type of the particle iterator
   * @param particlesBegin Iterator to the first particle. Iterated until no longer valid.
   */
  template <class ParticleIterator_T>
  void initializeNeighborList(ParticleIterator_T particlesBegin) {
    _neighborListsAoS.clear();
    for (auto iter = particlesBegin; iter.isValid(); ++iter) {
      _neighborListsAoS[&(*iter)];
    }
  }

  /**
   * Whether particle is relevant for tuning.
   * @return true
   */
  bool isRelevantForTuning() override { return true; }

  /**
   * Whether NeighborIdentificationFunctor allows Newton3. Is always allowed.
   * @return true
   */
  bool allowsNewton3() override { return true; }

  /**
   * Whether NeighborIdentificationFunctor allows non-newton3. This is not allowed if
gathering N3 lists.
   * @return true if gatherNewton3Lists is false, otherwise false.
   */
  bool allowsNonNewton3() override { return not _gatherNewton3Lists; }

  /**
   * AoSFunctor for interaction list generation.
   * @param i First particle
   * @param j Second particle
   * @param newton3 Whether AutoPas is using N3 or not for list generation.
   */
  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) override {
    using namespace autopas::utils::ArrayMath::literals;

    if (_gatherNewton3Lists and not newton3) [[unlikely]] {
      utils::ExceptionHandler::exception(
          "NeighborIdentificationFunctor should not be used with newton3=false and gatherNewton3Lists=true.");
    }

    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto dist = i.getR() - j.getR();

    double distSquared = utils::ArrayMath::dot(dist, dist);
    if (distSquared < _interactionLengthSquared) {
      // This is only thread-safe if all keys are inserted prior to applying this functor, as a rehash might lead to
      // dangling references, but at() protects against this: if the key exists, we can push_back safely; if a key
      // doesn't, an exception is thrown anyway.
      _neighborListsAoS.at(&i).push_back(&j);
      if (newton3 and not _gatherNewton3Lists) {
        _neighborListsAoS.at(&j).push_back(&i);
      }
    }
  }

  /**
   * SoAFunctor for verlet list generation. (single cell version)
   * @param soa the soa
   * @param newton3 Whether newton3 is used.
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool /*newton3*/) override {
    if (soa.size() == 0) return;

    auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    size_t numPart = soa.size();
    for (size_t i = 0; i < numPart; ++i) {
      if (ownedStatePtr[i] == OwnershipState::dummy) continue;
      auto &iList = _neighborListsAoS.at(ptrptr[i]);
      for (size_t j = i + 1; j < numPart; ++j) {
        if (ownedStatePtr[j] == OwnershipState::dummy) continue;
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _interactionLengthSquared) {
          iList.push_back(ptrptr[j]);
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
   * @param newton3 Whether newton3 is used.
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3) override {
    if (_gatherNewton3Lists and not newton3) [[unlikely]] {
      utils::ExceptionHandler::exception(
          "NeighborIdentificationFunctor should not be used with newton3=false and gatherNewton3Lists=true.");
    }

    if (soa1.size() == 0 or soa2.size() == 0) return;

    auto **const __restrict ptr1ptr = soa1.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedState1Ptr = soa1.template begin<Particle_T::AttributeNames::ownershipState>();

    auto **const __restrict ptr2ptr = soa2.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedState2Ptr = soa2.template begin<Particle_T::AttributeNames::ownershipState>();

    const size_t numPart1 = soa1.size();
    for (size_t i = 0; i < numPart1; ++i) {
      auto &iList = _neighborListsAoS.at(ptr1ptr[i]);
      if (ownedState1Ptr[i] == OwnershipState::dummy) continue;
      const size_t numPart2 = soa2.size();

      for (size_t j = 0; j < numPart2; ++j) {
        if (ownedState2Ptr[j] == OwnershipState::dummy) continue;
        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _interactionLengthSquared) {
          iList.push_back(ptr2ptr[j]);
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
   * Generates interaction list entries for the particle at index indexFirst in soa and
every potential neighbor
   * given by the Verlet list.
   *
   * @param soa the soa
   * @param indexFirst index of the particle whose neighbors are given in verletList
   * @param verletList indices of the potential neighbors of indexFirst
   * @param newton3 Whether newton3 is used.
   */
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, std::span<const size_t> verletList,
                        bool newton3) override {
    if (_gatherNewton3Lists and not newton3) [[unlikely]] {
      utils::ExceptionHandler::exception(
          "NeighborIdentificationFunctor should not be used with newton3=false and gatherNewton3Lists=true.");
    }

    if (soa.size() == 0 or verletList.empty()) return;

    auto **const __restrict ptrptr = soa.template begin<Particle_T::AttributeNames::ptr>();
    const double *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const double *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const double *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    if (ownedStatePtr[indexFirst] == OwnershipState::dummy) return;

    auto &firstList = _neighborListsAoS.at(ptrptr[indexFirst]);
    const double xFirst = xptr[indexFirst];
    const double yFirst = yptr[indexFirst];
    const double zFirst = zptr[indexFirst];

    for (const size_t j : verletList) {
      if (ownedStatePtr[j] == OwnershipState::dummy) continue;
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
  constexpr static std::array<typename Particle_T::AttributeNames, 5> getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 5>{
        Particle_T::AttributeNames::ptr, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
        Particle_T::AttributeNames::posZ, Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static std::array<typename Particle_T::AttributeNames, 5> getNeededAttr(std::false_type) {
    return getNeededAttr();
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
   * If true, each particle pair will appear in only one particle's list. Else, each
particle pair will appear in each
   * particle's list.
   */
  bool _gatherNewton3Lists{false};
};

}  // namespace autopas