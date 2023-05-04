/**
 * @file DynamicVLGeneratorFunctor.h
 * @author Luis Gall
 * @date 04.05.2023
 *
 * oriented on
 * @file VLCAllCellsGeneratorFunctor.h
 * @author tirgendetwas
 * @date 05.12.2020
 */

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/containers/verletListsCellBased/verletLists/NewVerletListHelpers.h"

namespace autopas {

template <class Particle>
class DynamicVLGeneratorFunctor : public Functor<Particle, DynamicVLGeneratorFunctor<Particle>> {
  using NeighborListsType = typename NewVerletListHelpers<Particle>::DynamicNeighborListsType;
  using SoAArraysType = typename Particle::SoAArraysType;

 public:

  DynamicVLGeneratorFunctor(NeighborListsType &neighborLists,
                           double interactionLength)
      : Functor<Particle, DynamicVLGeneratorFunctor<Particle>>(0.),
        _neighborLists(neighborLists),
        _interactionLengthSquared(interactionLength * interactionLength) {}

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override {
    utils::ExceptionHandler::exception(
        "DynamicVLGeneratorFunctor::allowsNewton3() is not implemented, because it should not be called.");
    return true;
  }

  bool allowsNonNewton3() override {
    utils::ExceptionHandler::exception(
        "DynamicVLGeneratorFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
    return true;
  }

  void AoSFunctor(Particle &i, Particle &j, bool /*newton3*/) override {
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto dist = utils::ArrayMath::sub(i.getR(), j.getR());

    double distsquare = utils::ArrayMath::dot(dist, dist);
    if (distsquare < _interactionLengthSquared) {
      // this is thread safe, only if particle i is accessed by only one
      // thread at a time. which is ensured, as particle i resides in a
      // specific cell and each cell is only accessed by one thread at a time
      // (ensured by traversals)
      // also the list is not allowed to be resized!

      _neighborLists.at(&i).first.push_back(&j);
      // no newton3 here, as AoSFunctor(j,i) will also be called if newton3 is disabled.
    }
  }

  /**
     * SoAFunctor for verlet list generation. (single cell version)
     * @param soa the soa
     * @param newton3 whether to use newton 3
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
    if (soa.getNumberOfParticles() == 0) return;

    auto **const __restrict ptrptr = soa.template begin<Particle::AttributeNames::ptr>();
    const double *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const double *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const double *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    size_t numPart = soa.getNumberOfParticles();
    for (unsigned int i = 0; i < numPart; ++i) {
      auto &currentList = _neighborLists.at(ptrptr[i]).first;

      for (unsigned int j = i + 1; j < numPart; ++j) {
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _interactionLengthSquared) {
          currentList.push_back(ptrptr[j]);
          if (not newton3) {
            // we need this here, as SoAFunctorSingle will only be called once for both newton3=true and false.
            _neighborLists.at(ptrptr[j]).first.push_back(ptrptr[i]);
          }
        }
      }
    }
  }

  /**
     * SoAFunctor for the verlet list generation. (two cell version)
     * @param soa1 soa of first cell
     * @param soa2 soa of second cell
     * @note newton3 is ignored here, as for newton3=false SoAFunctorPair(soa2, soa1) will also be called.
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool /*newton3*/) override {
    if (soa1.getNumberOfParticles() == 0 || soa2.getNumberOfParticles() == 0) return;

    auto **const __restrict ptr1ptr = soa1.template begin<Particle::AttributeNames::ptr>();
    const double *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const double *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const double *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();

    auto **const __restrict ptr2ptr = soa2.template begin<Particle::AttributeNames::ptr>();
    const double *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const double *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const double *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    size_t numPart1 = soa1.getNumberOfParticles();
    for (unsigned int i = 0; i < numPart1; ++i) {
      auto &currentList = _neighborLists.at(ptr1ptr[i]).first;

      size_t numPart2 = soa2.getNumberOfParticles();

      for (unsigned int j = 0; j < numPart2; ++j) {
        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _interactionLengthSquared) {
          currentList.push_back(ptr2ptr[j]);
        }
      }
    }
  }

  /**
     * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static std::array<typename Particle::AttributeNames, 4> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 4>{
        Particle::AttributeNames::ptr, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ};
  }

  /**
     * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 0>{/*Nothing*/};
  }

 private:
  NeighborListsType &_neighborLists;
  double _interactionLengthSquared;

};

}